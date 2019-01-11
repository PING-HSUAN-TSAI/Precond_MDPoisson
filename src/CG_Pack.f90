      subroutine CG(x_out,x_in,f,l,iterNum,tol)
      use constants
      use MD2D_Grid
      use CG_Var
      use State_Var
      use ERR_Var
      use Legendre
      implicit none
      
!     Solve Ax=f where A is SPD and is invoked by the routine ax()
!
!     Output:  x - vector of length n
!
!     Input:   f - vector of length n
!
!     Work arrays:   r,w,p,z  - vectors of length n
!
!     User-provided routine ax(w,z,n) returns  w := Az,
!
!     User-provided routine solveM(z,r,n) ) returns  z := M^-1 r,
!
!     User-provided array wght() is used to scale inner-products
!     of the form (p,z) = p'*wght*z     

      integer :: iterNum, i,j,k,method,dummy1,dummy2
      integer :: l, Nx, Ny
      real(kind=8):: tol
      real(kind=8):: x_coord, y_coord
      real(kind=8):: F_init, v_init
      real(kind=8):: alpha 
      real(kind=8):: rsold, rsnew, pAp
      real(kind=8)::  x_in(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8):: x_out(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8)::     f(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8)::    Ap(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8)::    Ax(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8)::     p(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8)::     r(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8):: glsc2

      Nx = PolyDegN_DM(1,1,l); Ny = PolyDegN_DM(2,1,l)

!      tol = 1E-8
      
!     initialize Ax matrix
      call axhelm2(Nx,Ny,x_in,Ax,l)

      
!     initial residue and p
      call add3s2(r,f,Ax,1.0,-1.0,Nx,Ny,l)
      call copy(p,r,Nx,Ny,l)

      call chk_amax('inr',r,Nx,Ny,l)
      call chk_amax('inp',p,Nx,Ny,l)

      rsold = glsc2(Nx,Ny,r,r,l)
      
!     start iteration
      do k=1,3*iterNum

         call axhelm2(Nx,Ny,p,Ap,l)

         pAp = glsc2(Nx,Ny,p,Ap,l)
         alpha = rsold / pAp

         call add2s2(x_in,p,alpha,Nx,Ny,l)
         call add2s2(r,Ap,-alpha,Nx,Ny,l)
      
         rsnew = glsc2(Nx,Ny,r,r,l)

         if (abs(rsnew) .lt. tol) then
      
            write(10,fmt="(i5,A36,es24.15)")k," iterations to reach the tolerance: ", tol
            exit
         endif
      
         call add3s2(p,r,p,1.0,(rsnew/rsold),Nx,Ny,l)
      
         rsold = rsnew
      
      enddo! iteration

      call copy(x_out,x_in,Nx,Ny,l)
      
      write(10,9999) k,tol
      
9999 format(' ',' ', 'CG   : iteration#',i5,1p3e12.4)
         return
      end subroutine CG
!----------------------------------------------------------------------
      subroutine axhelm2(LD1,LD2,u,Au,l)
      use MD2D_Grid    ! MD2D_Grid.f90
      use State_Var    ! State_Var.f90
      use Legendre     ! Legendre.f90
      use CG_Var       ! CG_Var.f90
      implicit none
      
      integer :: i, j, LD1,LD2, NDFix
      real(kind=8) :: u(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: Au(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: dudx(0:LD1,0:LD2), dudy(0:LD1,0:LD2)
      real(kind=8) :: tmp1(0:LD1,0:LD2,1:TotNum_DM), tmp2(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: vEdge(0:LD1,1:4,1:TotNum_DM), dvdnEdge(0:LD1,1:4,1:TotNum_DM) 
      real(kind=8) :: BvEdge(0:LD1,1:4,1:TotNum_DM) 
      real(kind=8) :: BITF1(0:LD1,1:4,1:TotNum_DM), BITF2(0:LD1,1:4,1:TotNum_DM)
      integer:: l
      
      !write(*,*)'l in cg=',l
      !real(kind=8):: Bq_loc, x_coord, y_coord
      !real(kind=8):: nx, ny, alpha, beta, b_loc, g_ext
      
      !open(180,file="differential.text")
      
      !do DDK=1,TotNum_DM
      !ND1=PolyDegN_DM(1,DDK,l)
      !ND2=PolyDegN_DM(2,DDK,l)
      !write(*,*)LD1,LD2
      !do j=0,LD2
      !do i=0,LD1
      !write(*,*)i,j,u(i,j,DDK),Au(i,j,DDK)
      !enddo
      !enddo
      !enddo
      
      Au=0.0d0; dudx=0.0d0; dudy=0.0d0
      tmp1=0.0d0; tmp2=0.0d0
      vEdge=0.0d0; dvdnEdge=0.0d0; BvEdge=0.0d0
      BITF1=0.0d0; BITF2=0.0d0

      write(*,*)'u',maxval(u)
      
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l)
         ND2=PolyDegN_DM(2,DDK,l)
      !   compute gradient x
      !open(181,file="u.text",position="append",action="write")
!      do j=0,ND2
!      do i=0,ND1
!      write(*,*)i,j,u(i,j,DDK)
!!      !write(181,fmt="(2I3,1x,1E11.5e3)")i,j,u(i,j,DDK)
!      enddo
!      enddo
      !close(181)
         dudx(0:ND1,0:ND2) = Matmul(Diff_xi1(0:ND1,0:ND1,ND1), &
                                           u(0:ND1,0:ND2,DDK))
      
      
         dudy(0:ND1,0:ND2) = &
                    Matmul(u(0:ND1,0:ND2,DDK),Diff_xi2(0:ND2,0:ND2,ND2))
!      write(*,*)'dudx,dudy',maxval(abs(dudx(0:ND1,0:ND2))),maxval(abs(dudy(0:ND1,0:ND2)))      
      !write(*,*)"cg_dqdx,cg_dqdy"
      
!      do j=0,ND2
!      do i=0,ND1
!      write(*,*) i,j,dudx(i,j),dudy(i,j)
      !write(180,fmt="(2I3,1x,3E11.5e3)") i,j,cg_dqdx(i,j),cg_dqdy(i,j)
!      enddo
!      enddo
      !close(180)
         tmp1(0:ND1,0:ND2,DDK) = dxi1_dx1(0:ND1,0:ND2,DDK,l) * &
                                     dudx(0:ND1,0:ND2) &
                               + dxi2_dx1(0:ND1,0:ND2,DDK,l) * &
                                     dudy(0:ND1,0:ND2)
      
         tmp2(0:ND1,0:ND2,DDK) = dxi1_dx2(0:ND1,0:ND2,DDK,l) * &
                                     dudx(0:ND1,0:ND2) &
                               + dxi2_dx2(0:ND1,0:ND2,DDK,l) * &
                                     dudy(0:ND1,0:ND2)
      
!      write(*,*)'tmp1,tmp2',maxval(abs(tmp1(0:ND1,0:ND2,DDK))),maxval(abs(tmp2(0:ND1,0:ND2,DDK)))
      enddo !DDK  

      ! assign u and dx_cgdn on the edge number
      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,l);ND2=PolyDegN_DM(2,DDK,l)
      
         do Edge_Num=1,4
            select case (Edge_Num)
               case (2,4)
      
               ND=ND2; NDFix=0; if (Edge_Num .eq. 2) NDFix=ND1
      
               vEdge(0:ND,Edge_Num,DDK) = u(NDFix,0:ND,DDK)
               dvdnEdge(0:ND,Edge_Num,DDK) &
                  = Norvec_x1(0:ND,Edge_Num,DDK,l) * tmp1(NDFix,0:ND,DDK) &
                  + Norvec_x2(0:ND,Edge_Num,DDK,l) * tmp2(NDFix,0:ND,DDK)
               BvEdge(0:ND,Edge_Num,DDK) &
                  =  bnd_alpha(0:ND,Edge_Num,DDK,l) * vEdge(0:ND,Edge_Num,DDK) &
                  +   bnd_beta(0:ND,Edge_Num,DDK,l) * bEdge(0:ND,Edge_Num,DDK,l) &
                  *   dvdnEdge(0:ND,Edge_Num,DDK)

               case(1,3)
      
               ND=ND1; NDFix=0; if (Edge_Num .eq. 3) NDFix=ND2
      
               vEdge(0:ND,Edge_Num,DDK) = u(0:ND,NDFix,DDK)
               dvdnEdge(0:ND,Edge_Num,DDK) &
                  = Norvec_x1(0:ND,Edge_Num,DDK,l) * tmp1(0:ND,NDFix,DDK) &
                  + Norvec_x2(0:ND,Edge_Num,DDK,l) * tmp2(0:ND,NDFix,DDK)
      
               BvEdge(0:ND,Edge_Num,DDK) &
                  =  bnd_alpha(0:ND,Edge_Num,DDK,l) * vEdge(0:ND,Edge_Num,DDK) &
                  +   bnd_beta(0:ND,Edge_Num,DDK,l) * bEdge(0:ND,Edge_Num,DDK,l) &
                  *   dvdnEdge(0:ND,Edge_Num,DDK)
      
            end select !Edge_Num
      
         enddo !Edge_Num
      
      enddo ! DDK
      
      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,l);ND2=PolyDegN_DM(2,DDK,l)
      
         do Edge_Num=1,4
            select case (Edge_Num)
               case (2,4)
      
               ND=ND2; NDFix=0; if (Edge_Num .eq. 2) NDFix=ND1
      
               select case (BC_Type(Edge_Num,DDK))
      
                  case(0) ! internal BC
                  DDK_Connect   = DM_Connect(1,Edge_Num,DDK)
                  Edge_Connect  = DM_Connect(2,Edge_Num,DDK)
                  Patch_Type    = DM_Connect(3,Edge_Num,DDK)
                  select case (Patch_Type)
                     case(1)
                     do j=0,ND
      
                        BITF1(j,Edge_Num,DDK) = vEdge(j,Edge_Num,DDK) &
                                              - vEdge(j,Edge_Connect,DDK_Connect)
      
                        BITF2(j,Edge_Num,DDK) = a(DDK) * dvdnEdge(j,Edge_Num,DDK) &
                                              + a(DDK_Connect) &
                                              * dvdnEdge(j,Edge_Connect,DDK_Connect)
      
                     enddo
      
                     case (-1) ! Reverse Patching
                     do j=0,ND
      
                        BITF1(j,Edge_Num,DDK) = vEdge(j,Edge_Num,DDK) &
                                              - vEdge(ND-j,Edge_Connect,DDK_Connect)
      
                        BITF2(j,Edge_Num,DDK) = a(DDK) * dvdnEdge(j,Edge_Num,DDK) &
                                              + a(DDK_Connect) &
                                              * dvdnEdge(ND-j,Edge_Connect,DDK_Connect)
      
                     enddo
      
                  end select ! Patch_type
      
               end select ! bc_type
      
            case(1,3)
      
               ND=ND1; NDFix=0; if (Edge_Num .eq. 3) NDFix=ND2
      
               select case (BC_Type(Edge_Num,DDK))
                  case (0)
                  DDK_Connect  = DM_Connect(1,Edge_Num,DDK)
                  Edge_Connect = DM_Connect(2,Edge_Num,DDK)
                  Patch_Type   = DM_Connect(3,Edge_Num,DDK)
      !           write(*,10000)'P',DDK,Edge_Num,DDK_Connect,Edge_Connect,Patch_Type
      !   10000 format(1A3,1x,4i5)
                  select case (Patch_Type)
                     case(1)
                     do i=0,ND
                        BITF1(i,Edge_Num,DDK) = vEdge(i,Edge_Num,DDK) &
                                              - vEdge(i,Edge_Connect,DDK_Connect)
      
                        BITF2(i,Edge_Num,DDK) = a(DDK) * dvdnEdge(i,Edge_Num,DDK) &
                                              + a(DDK_Connect) &
                                              * dvdnEdge(i,Edge_Connect,DDK_Connect)
      
                     enddo
                     case (-1) ! Reverse Patching
                     do i=0,ND
      
                        BITF1(i,Edge_Num,DDK) = vEdge(i,Edge_Num,DDK) &
                                              - vEdge(ND-i,Edge_Connect,DDK_Connect)
      
                        BITF2(i,Edge_Num,DDK) = a(DDK) * dvdnEdge(i,Edge_Num,DDK) &
                                              + a(DDK_Connect) &
                                              * dvdnEdge(ND-i,Edge_Connect,DDK_Connect)
      
                     enddo ! i
                  end select ! Patch_Type
               end select ! BC_Type
            end select !Edge_Num
          enddo !Edge_Num
      enddo ! DDK
      
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l)
         ND2=PolyDegN_DM(2,DDK,l)
      
!     compute divergence b*F  (div dot bF)
         dudx(0:ND1,0:ND2) = (dxi1_dx1(0:ND1,0:ND2,DDK,l) *  tmp1(0:ND1,0:ND2,DDK)  &
                           +  dxi1_dx2(0:ND1,0:ND2,DDK,l) *  tmp2(0:ND1,0:ND2,DDK)) &
                           *   Jacobin(0:ND1,0:ND2,DDK,l) * b_pts(0:ND1,0:ND2,DDK,l)
      
         dudy(0:ND1,0:ND2) = (dxi2_dx1(0:ND1,0:ND2,DDK,l) *  tmp1(0:ND1,0:ND2,DDK)  &
                           +  dxi2_dx2(0:ND1,0:ND2,DDK,l) *  tmp2(0:ND1,0:ND2,DDK)) &
                           *   Jacobin(0:ND1,0:ND2,DDK,l) * b_pts(0:ND1,0:ND2,DDK,l)
      
         Au(0:ND1,0:ND2,DDK) = a_pts(0:ND1,0:ND2,DDK,l) &
                             * ( - Matmul( Diff_xi1(0:ND1,0:ND1,ND1), &
                                            dudx(0:ND1,0:ND2)) &
                                 - Matmul( dudy(0:ND1,0:ND2), &
                                          Diff_xi2(0:ND2,0:ND2,ND2)) )
      
      ! add penalty of boundary condition
      
      ! side 4
         Edge_Num=4
         select case (BC_Type(Edge_Num,DDK))
            case(1)
            do j = 0,ND2
               Au(0:ND1,j,DDK)  =   Au(0:ND1,j,DDK) &
                                + (    tauD(0:ND1,j,Edge_Num,DDK,l) &
                                - tau_tilde(0:ND1,j,Edge_Num,DDK,l)) &
                                * BvEdge(j,Edge_Num,DDK)
               Au(0,0:ND2,DDK)  =     Au(0,0:ND2,DDK) &
                                + (     tauD2(0:ND2,j,Edge_Num,DDK,l) &
                                *   BvEdge(j,Edge_Num,DDK) )
          
            enddo
      
            case(2,3)
            do j = 0,ND2
               Au(0:ND1,j,DDK) =   Au(0:ND1,j,DDK) &
                               +    tauND1(0:ND1,j,Edge_Num,DDK,l) &
                               * BvEdge(j,Edge_Num,DDK)
      
            enddo
            case(0)
      
!     add interface penalty : first set of dirichlet
            do j=0,ND2
               Au(0:ND1,j,DDK) =  Au(0:ND1,j,DDK) &
                               +     tau1(0:ND1,j,Edge_Num,DDK,l)  &
                               *  BITF1(j,Edge_Num,DDK)
               Au(0,0:ND2,DDK) =  Au(0,0:ND2,DDK) &
                               +     tau3(j,0:ND2,Edge_Num,DDK,l) &
                               *  BITF1(j,Edge_Num,DDK)
            enddo

            Au(0,0:ND2,DDK) =   Au(0,0:ND2,DDK) &
                            +      tau2(0:ND2,Edge_Num,DDK,l) &
                            *   BITF1(0:ND2,Edge_Num,DDK)
!     add interface penalty : second set of dirichlet
            Au(0,0:ND2,DDK) =      Au(0,0:ND2,DDK) &
                            +   Sigma_tild(0:ND2,Edge_Num,DDK,l) &
                            *      BITF1(0:ND2,Edge_Num,DDK)
!     add interface penalty : neumann
            Au(0,0:ND2,DDK) =  Au(0,0:ND2,DDK) &
                            + SigmaHat(0:ND2,Edge_Num,DDK,l) &
                            *  BITF2(0:ND2,Edge_Num,DDK)
      
         end select
      
      ! side 2
         Edge_Num=2
         select case (BC_Type(Edge_Num,DDK))
            case(1)
            do j=0,ND2
               Au(0:ND1,j,DDK) =     Au(0:ND1,j,DDK) &
                               + (      tauD(0:ND1,j,Edge_Num,DDK,l) &
                               -   tau_tilde(0:ND1,j,Edge_Num,DDK,l)) &
                               *   BvEdge(j,Edge_Num,DDK)
               Au(ND1,0:ND2,DDK) =     Au(ND1,0:ND2,DDK) &
                                 + (     tauD2(0:ND2,j,Edge_Num,DDK,l) &
                                 *   BvEdge(j,Edge_Num,DDK) )
            enddo

            case(2,3)
            do j=0,ND2
               Au(0:ND1,j,DDK) =    Au(0:ND1,j,DDK) &
                               +     tauND1(0:ND1,j,Edge_Num,DDK,l) &
                               *  BvEdge(j,Edge_Num,DDK)
      
            enddo
            case(0)
      
!     add interface penalty : first set of dirichlet
            do j=0,ND2
               Au(0:ND1,j,DDK) =  Au(0:ND1,j,DDK) &
                               +     tau1(0:ND1,j,Edge_Num,DDK,l)  &
                               *  BITF1(j,Edge_Num,DDK)
      
               Au(ND1,0:ND2,DDK) =   Au(ND1,0:ND2,DDK) &
                                 + tau3(j,0:ND2,Edge_Num,DDK,l) &
                                 *  BITF1(j,Edge_Num,DDK)
            enddo
               Au(ND1,0:ND2,DDK) = Au(ND1,0:ND2,DDK) &
                                 +    tau2(0:ND2,Edge_Num,DDK,l) &
                                 * BITF1(0:ND2,Edge_Num,DDK)
!     add interface penalty : second set of dirichlet
               Au(ND1,0:ND2,DDK) =     Au(ND1,0:ND2,DDK) &
                                 +Sigma_tild(0:ND2,Edge_Num,DDK,l) &
                                 *     BITF1(0:ND2,Edge_Num,DDK)
!     add interface penalty : neumann
               Au(ND1,0:ND2,DDK) =  Au(ND1,0:ND2,DDK) &
                                 + SigmaHat(0:ND2,Edge_Num,DDK,l) &
                                 *  BITF2(0:ND2,Edge_Num,DDK)
      
         end select
      
      ! side 1
         Edge_Num=1
         select case (BC_Type(Edge_Num,DDK))
            case(1)
            do i=0,ND1
               Au(i,0:ND2,DDK) =   Au(i,0:ND2,DDK) &
                               + (    tauD(i,0:ND2,Edge_Num,DDK,l) &
                               - tau_tilde(i,0:ND2,Edge_Num,DDK,l)) &
                               * BvEdge(i,Edge_Num,DDK)
      
               Au(0:ND1,0,DDK) =     Au(0:ND1,0,DDK) &
                               + (     tauD2(0:ND1,i,Edge_Num,DDK,l) &
                               *   BvEdge(i,Edge_Num,DDK) )
            enddo
            case(2,3)
            do i=0,ND1
               Au(i,0:ND2,DDK) =   Au(i,0:ND2,DDK) &
                               +    tauND1(i,0:ND2,Edge_Num,DDK,l) &
                               * BvEdge(i,Edge_Num,DDK)
      
            enddo
            case(0)
      
!     add interface penalty : first set of dirichlet
            do i=0,ND1
               Au(i,0:ND2,DDK) = Au(i,0:ND2,DDK) &
                               +    tau1(i,0:ND2,Edge_Num,DDK,l)  &
                               * BITF1(i,Edge_Num,DDK)
               Au(0:ND1,0,DDK) = Au(0:ND1,0,DDK) &
                               +    tau3(i,0:ND1,Edge_Num,DDK,l) &
                               * BITF1(i,Edge_Num,DDK)
            enddo

            Au(0:ND1,0,DDK) = Au(0:ND1,0,DDK) &
                            +    tau2(0:ND1,Edge_Num,DDK,l) &
                            * BITF1(0:ND1,Edge_Num,DDK)
!     add interface penalty : second set of dirichlet
            Au(0:ND1,0,DDK) =    Au(0:ND1,0,DDK) &
                            + Sigma_tild(0:ND1,Edge_Num,DDK,l) &
                            *    BITF1(0:ND1,Edge_Num,DDK)
!     add interface penalty : neumann
            Au(0:ND1,0,DDK) =  Au(0:ND1,0,DDK) &
                            + SigmaHat(0:ND1,Edge_Num,DDK,l) &
                            *  BITF2(0:ND1,Edge_Num,DDK)
         end select
      
      ! side 3
         Edge_Num=3
         select case (BC_Type(Edge_Num,DDK))
            case(1)
            do i=0,ND1
               Au(i,0:ND2,DDK) =   Au(i,0:ND2,DDK) &
                               + (    tauD(i,0:ND2,Edge_Num,DDK,l) &
                               - tau_tilde(i,0:ND2,Edge_Num,DDK,l)) &
                               * BvEdge(i,Edge_Num,DDK)
      
               Au(0:ND1,ND2,DDK) =     Au(0:ND1,ND2,DDK) &
                                 + (     tauD2(0:ND1,i,Edge_Num,DDK,l) &
                                 *   BvEdge(i,Edge_Num,DDK) )
      
            enddo
      
            case(2,3)
            do i=0,ND1
               Au(i,0:ND2,DDK) =   Au(i,0:ND2,DDK) &
                               +    tauND1(i,0:ND2,Edge_Num,DDK,l) &
                               * BvEdge(i,Edge_Num,DDK)
      
            enddo
            case(0)
      
!     add interface penalty : first set of dirichlet
            do i=0,ND1
               Au(i,0:ND2,DDK) = Au(i,0:ND2,DDK) &
                               +    tau1(i,0:ND2,Edge_Num,DDK,l)  &
                               * BITF1(i,Edge_Num,DDK)
               Au(0:ND1,ND2,DDK) = Au(0:ND1,ND2,DDK) &
                                 +    tau3(i,0:ND1,Edge_Num,DDK,l) &
                                 * BITF1(i,Edge_Num,DDK)
      
            enddo
            Au(0:ND1,ND2,DDK) = Au(0:ND1,ND2,DDK) &
                              +    tau2(0:ND1,Edge_Num,DDK,l) &
                              * BITF1(0:ND1,Edge_Num,DDK)
!     add interface penalty : second set of dirichlet
            Au(0:ND1,ND2,DDK) =    Au(0:ND1,ND2,DDK) &
                              + Sigma_tild(0:ND1,Edge_Num,DDK,l) &
                              *    BITF1(0:ND1,Edge_Num,DDK)
!     add interface penalty : neumann
            Au(0:ND1,ND2,DDK) =  Au(0:ND1,ND2,DDK) &
                              + SigmaHat(0:ND1,Edge_Num,DDK,l) &
                              *  BITF2(0:ND1,Edge_Num,DDK)
         end select
      
      enddo !DDK
      
!     multiply mass matrix
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l)
         ND2=PolyDegN_DM(2,DDK,l)
         do j=0,ND2
            do i=0,ND1
               Au(i,j,DDK) =    Au(i,j,DDK) &
                           * LGLWeights(j,ND2) &
                           * LGLWeights(i,ND1)
      
            enddo
         enddo
      enddo
      
      return
      
      end subroutine axhelm2
!---------------------------------------------------------------------------------------------
      subroutine inner_product(LD1,LD2,vector1,vector2,tmp,l)

      use MD2D_Grid    ! MD2D_Grid.f90
      use CG_Var       ! CG_Var.f90

      implicit none
      
      integer:: i,j, LD1, LD2, l
      real(kind=8) :: vector1(0:LD1,0:LD2,1:TotNum_DM) 
      real(kind=8) :: vector2(0:LD1,0:LD2,1:TotNum_DM) 
      real(kind=8) :: tmp 
      
      tmp=0.0
      
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l)
         ND2=PolyDegN_DM(2,DDK,l)
         do j=0,ND2
            do i=0,ND1
            tmp = tmp + (vector1(i,j,DDK) * vector2(i,j,DDK))
            enddo
         enddo
      enddo

      return
      end subroutine
!-----------------------------------------------------------------------------------------------
      
