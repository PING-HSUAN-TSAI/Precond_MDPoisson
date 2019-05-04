!------------------------------------------------------------- 
      subroutine Init_Material_Parameters(mgl)
      use State_Var
      use MD2D_Grid
      implicit none

!     Declare local arguments
      integer:: lid, count
      integer:: i,j
      integer:: l,mgl
      
      lid=80
      open(lid,file='Material.in',form='formatted', status='unknown')

         read(lid,*) !'============================================'
         do DDK=1,TotNum_DM
            read(lid,1000) count, a(DDK), b(DDK)
         enddo
         read(lid,*) !'============================================'

      close(lid)
      
!     Assign grid point values of a and b 
      do l = mgl,1,-1
         do DDK = 1, TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
            a_pts(0:ND1,0:ND2,DDK,l)=a(DDK)
            b_pts(0:ND1,0:ND2,DDK,l)=b(DDK)
      
         enddo !DDK
      enddo !l
      
      
      write(10,*)'Complete Initializing Material Parameters'
      
1000 format(i6,2f10.4)
      
      end subroutine Init_Material_Parameters
!------------------------------------------------------------- 
      subroutine Init_BC_Variables(mgl)
      use MD2D_Grid
      use State_Var
      implicit none
      integer:: N_max,mgl
        
!      N_max=maxval(PolyDegN_Max(1:2));
      N_max = PolyDegN_DM(1,1,mgl)
      call alloc_mem_BC_Var(N_max,TotNum_DM,mgl)
        
      call Init_BC_Type
      
      call Init_Boundary_alpha_beta(mgl)
        
      call Init_Penalty_Parameters(mgl)
      
      end subroutine Init_BC_Variables
!------------------------------------------------------------- 
      Subroutine Init_BC_Type
      use MD2D_Grid
      use State_Var
      implicit none
      integer:: lid, Num_DM_BC, tmp_int, ierr
        
! subroutine begin
      lid=81
      open(lid, file='BC.in', form='formatted', status='old')
         read(lid,*) !'==============================================='
         read(lid,*) Num_DM_BC
         read(lid,*) strength
         if (Num_DM_BC .ne. TotNum_DM) then 
            write(*,*)'Message from State_Pack.f90'
            write(*,*)'Inconsistent specificition of Total Number of BC'
            write(*,*)'Stop!'
            close(lid)
            stop
         endif
      
         read(lid,*) !'===============================================' 
         do DDK =1,TotNum_DM 
            read(lid,*) tmp_int, BC_type(1,DDK),BC_type(2,DDK),&
                                 BC_type(3,DDK),BC_type(4,DDK)  
            read(lid,*) !'============================================='
         enddo 
        
      close(lid)
      return 
        
      end subroutine Init_BC_Type
!-------------------------------------------------------------
      subroutine init_boundary_alpha_beta(mgl)
      use MD2D_Grid
      use Legendre
      use State_Var
      use INPUT
      implicit none

      integer:: i, j, NDFix
      real(kind=8):: OmegaSelf, OmegaSurr
      real(kind=8):: c_bar_default
      integer:: l,mgl
      
      c_sigma(1) = param(34)
      c_sigma(2) = param(35)
      write(10,*)'Interface Parameters:',c_sigma(1),c_sigma(2)

      do l=mgl,1,-1
         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
!     assign c_bar = 2 as default
!     so c_bar=2       for edge points and 
!        c_bar=2*2=4   for vertex points
            c_bar_default=1.0d0
            c_bar(0:ND2,2,DDK,l) = c_bar_default;
            c_bar(0:ND2,4,DDK,l) = c_bar_default;
            c_bar(0:ND1,1,DDK,l) = c_bar_default;
            c_bar(0:ND1,3,DDK,l) = c_bar_default;
      
            c_bar(0:ND2:ND2,2,DDK,l)=c_bar(0:ND2:ND2,2,DDK,l)*2.d0;
            c_bar(0:ND2:ND2,4,DDK,l)=c_bar(0:ND2:ND2,4,DDK,l)*2.d0;
            c_bar(0:ND1:ND1,1,DDK,l)=c_bar(0:ND1:ND1,1,DDK,l)*2.d0;
            c_bar(0:ND1:ND1,3,DDK,l)=c_bar(0:ND1:ND1,3,DDK,l)*2.d0;

!     Assign b grid value on the edges        
            bEdge(0:ND1,1,DDK,l)=b_pts(0:ND1,0,DDK,l)
            bEdge(0:ND2,2,DDK,l)=b_pts(ND1,0:ND2,DDK,l)
            bEdge(0:ND1,3,DDK,l)=b_pts(0:ND1,ND2,DDK,l)
            bEdge(0:ND2,4,DDK,l)=b_pts(0,0:ND2,DDK,l)
      
         enddo ! DDK
      enddo ! l
      
      do l = mgl,1,-1
         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
            do Edge_Num=1,4 

!     Assign beta and alpha
               select case (Edge_Num)
                  case (2,4)  ! Edge 2 and 4 are curves xi= \pm 1
                     ND=ND2; OmegaSelf=LGLWeights(0,ND1)
                     NDFix=0; if (Edge_Num .eq. 2 ) NDFix = ND1
      
                     do j=0,ND
!     compute J * (\nabla \xi \cdot \nabla \eta)
                        J_xi_eta(j,Edge_Num,DDK,l) = Jacobin(NDFix,j,DDK,l) &
                         * ( dxi2_dx1(NDFix,j,DDK,l) * dxi1_dx1(NDFix,j,DDK,l) &
                         +   dxi2_dx2(NDFix,j,DDK,l) * dxi1_dx2(NDFix,j,DDK,l) )
!     compute normal vector \cdot \nabla \xi
                        Nor_dot_gradxi(j,Edge_Num,DDK,l) = &
                           NorVec_x1(j,Edge_Num,DDK,l) * dxi1_dx1(NDFix,j,DDK,l) &
                         + NorVec_x2(j,Edge_Num,DDK,l) * dxi1_dx2(NDFix,j,DDK,l)
!     compute normal vector \cdot \nabla \eta
                        Nor_dot_gradeta(j,Edge_Num,DDK,l) = &
                           NorVec_x1(j,Edge_Num,DDK,l) * dxi2_dx1(NDFix,j,DDK,l) &
                         + NorVec_x2(j,Edge_Num,DDK,l) * dxi2_dx2(NDFix,j,DDK,l)
                     enddo
      
                  case (1,3)  ! Edge 1 and 3 are curves eta= \pm 1
                     ND=ND1; OmegaSelf=LGLWeights(0,ND2)
                     NDFix=0; if (Edge_Num .eq. 3 ) NDFix = ND2
      
                     do i=0,ND
      
!     compute J * (\nabla \xi \cdot \nabla \eta)
                        J_eta_xi(i,Edge_Num,DDK,l) = Jacobin(i,NDFix,DDK,l) &
                        * ( dxi2_dx1(i,NDFix,DDK,l) * dxi1_dx1(i,NDFix,DDK,l) &
                        +   dxi2_dx2(i,NDFix,DDK,l) * dxi1_dx2(i,NDFix,DDK,l) )
      
!     compute normal vector \cdot \nabla \xi
                        Nor_dot_gradxi(i,Edge_Num,DDK,l) = &
                           NorVec_x1(i,Edge_Num,DDK,l) * dxi1_dx1(i,NDFix,DDK,l) &
                         + NorVec_x2(i,Edge_Num,DDK,l) * dxi1_dx2(i,NDFix,DDK,l)
      
!compute normal vector \cdot \nabla \eta
                        Nor_dot_gradeta(i,Edge_Num,DDK,l) = &
                           NorVec_x1(i,Edge_Num,DDK,l) * dxi2_dx1(i,NDFix,DDK,l) &
                         + NorVec_x2(i,Edge_Num,DDK,l) * dxi2_dx2(i,NDFix,DDK,l)
                     enddo
               end select! Edge_Num
      
               select case (BC_Type(Edge_Num,DDK))
                  case (1) ! Dirichlet BC 
                     bnd_alpha(0:ND,Edge_Num,DDK,l) = 1.0
                     bnd_beta(0:ND,Edge_Num,DDK,l) = 0.0

                  case (2) ! Neumann BC 
                     bnd_alpha(0:ND,Edge_Num,DDK,l) = 0.0
                     bnd_beta(0:ND,Edge_Num,DDK,l) = 1.0
      
                  case (3) ! Robin BC 
                     bnd_alpha(0:ND,Edge_Num,DDK,l) = 1.d0
                     bnd_beta(0:ND,Edge_Num,DDK,l) = 1.d0
      
               end select ! BC_Type
            enddo ! Edge_Num
         enddo ! DDK
      enddo ! l
      end subroutine 
!------------------------------------------------------------- 
      subroutine Init_Penalty_Parameters(mgl)
      use MD2D_Grid
      use Legendre
      use State_Var
      use INPUT
      implicit none
      
      integer:: i, j, NDFix, l,mgl
      real(kind=8):: omega, denum, OmegaSurr,OmegaSelf
      
      c_tau = param(33)
      write(10,*)'Dirichlet BC Parameter:',c_tau
      
      do l = mgl,1,-1
         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         
            do Edge_Num=1,4
               select case (Edge_Num)
                  case (2,4)
                     ND=ND2; omega = LGLWeights(0,ND1)
                     e_first(0,DDK,l)=1
                     e_end(ND1,DDK,l)=1
      
                  case (1,3)
                     ND=ND1; omega = LGLWeights(0,ND2)
                     e_first(0,DDK,l)=1
                     e_end(ND2,DDK,l)=1
               end select
      
               select case (BC_Type(Edge_Num,DDK))
!----------------------------------------------------------------------
                  case(1) ! Dirichlet
                     if (Edge_Num .eq. 3) then
                        do i=0,ND
!     Penalty for orthogonal part
                           tauD(i,0:ND2,Edge_Num,DDK,l) = &
                              - Diff_xi2(0:ND2,ND2,ND2) &
                              *      a_pts(i,ND2,DDK,l) &
                              *  JacNorVec(i,Edge_Num,DDK,l) &
                              *  NorVec_mg(i,Edge_Num,DDK,l) &
                              /  LGLWeights(0:ND2,ND2)
!     Penalty for nonorthogonal part
                           tauD2(0:ND1,i,Edge_Num,DDK,l) = &
                              - Diff_xi1(i,0:ND1,ND1) &
                              *      a_pts(i,ND2,DDK,l) &
                              *   J_eta_xi(i,Edge_Num,DDK,l) &
                              * LGLWeights(i,ND1) &
                              / (LGLWeights(0:ND1,ND1) &
                              *   omega )
!     Penalty for ensuring positive definite
                           tau_tilde(i,0:ND2,Edge_Num,DDK,l) = &
                              - c_tau * e_end(0:ND2,DDK,l) &
                              *   JacNorVec(i,Edge_Num,DDK,l) &
                              *       a_pts(i,ND2,DDK,l) &
                              *   NorVec_mg(i,Edge_Num,DDK,l) &
                              / (omega**2 * c_bar(i,Edge_Num,DDK,l))
                        enddo
!                     do j=0,ND2
!                        do i=0,ND1
!                        write(10,*)l,i,j,DDK,tauD(i,j,Edge_Num,DDK,l)
!                        enddo
!                     enddo
                     endif
      
                     if (Edge_Num .eq. 1) then
                        do i=0,ND
!     Penalty for orthogonal part
                           tauD(i,0:ND2,Edge_Num,DDK,l) = &
                              Diff_xi2(0:ND2,0,ND2) &
                           *      a_pts(i,0,DDK,l) &
                           *  JacNorVec(i,Edge_Num,DDK,l) &
                           *  NorVec_mg(i,Edge_Num,DDK,l) &
                           / LGLWeights(0:ND2,ND2)
!     Penalty for nonorthogonal part
                           tauD2(0:ND1,i,Edge_Num,DDK,l) = &
                              Diff_xi1(i,0:ND1,ND1) &
                           *       a_pts(i,0,DDK,l) &
                           *    J_eta_xi(i,Edge_Num,DDK,l)  &
                           *  LGLWeights(i,ND1) &
                           / (LGLWeights(0:ND1,ND1) &
                           *   omega )
!     Penalty for ensuring positive definite
                           tau_tilde(i,0:ND2,Edge_Num,DDK,l) = &
                              - c_tau * e_first(0:ND2,DDK,l) &
                           *     JacNorVec(i,Edge_Num,DDK,l) &
                           *         a_pts(i,0,DDK,l) &
                           *     NorVec_mg(i,Edge_Num,DDK,l) &
                           / (omega**2 * c_bar(i,Edge_Num,DDK,l))
                        enddo
                     endif
      
                     if (Edge_Num .eq. 2) then
                        do j=0,ND
!     Penalty for orthogonal part
                           tauD(0:ND1,j,Edge_Num,DDK,l) = &
                              - Diff_xi1(ND1,0:ND1,ND1) &
                           *      a_pts(ND1,j,DDK,l) &
                           *  JacNorVec(j,Edge_Num,DDK,l) &
                           *  NorVec_mg(j,Edge_Num,DDK,l) &
                           / LGLWeights(0:ND1,ND1)
!     Penalty for nonorthogonal part
                           tauD2(0:ND2,j,Edge_Num,DDK,l) = &
                              - Diff_xi2(0:ND2,j,ND2) &
                           *      a_pts(ND1,j,DDK,l) &
                           *   J_xi_eta(j,Edge_Num,DDK,l) &
                           * LGLWeights(j,ND2) &
                           / (LGLWeights(0:ND2,ND2) &
                           *   omega )
!     Penatly for ensuring positive definite
                           tau_tilde(0:ND1,j,Edge_Num,DDK,l) = &
                              - c_tau * e_end(0:ND1,DDK,l) &
                           *   JacNorVec(j,Edge_Num,DDK,l) &
                           *     a_pts(ND1,j,DDK,l) &
                           *   NorVec_mg(j,Edge_Num,DDK,l) &
                           / (omega**2 * c_bar(j,Edge_Num,DDK,l))
                        enddo
                     endif

                     if (Edge_Num .eq. 4) then
                        do j=0,ND
!     Penalty for orthogonal part
                           tauD(0:ND1,j,Edge_Num,DDK,l) = &
                              Diff_xi1(0,0:ND1,ND1) &
                           *      a_pts(0,j,DDK,l) &
                           *  JacNorVec(j,Edge_Num,DDK,l) &
                           *  NorVec_mg(j,Edge_Num,DDK,l) &
                           / LGLWeights(0:ND1,ND1)
!     Penalty for nonorthogonal part
                           tauD2(0:ND2,j,Edge_Num,DDK,l) = &
                              Diff_xi2(0:ND2,j,ND2) &
                           *      a_pts(0,j,DDK,l) &
                           *   J_xi_eta(j,Edge_Num,DDK,l) &
                           * LGLWeights(j,ND2) &
                           / (LGLWeights(0:ND2,ND2) &
                           *   omega )
!     Penalty for ensuring positive definite
                           tau_tilde(0:ND1,j,Edge_Num,DDK,l) = &
                              - c_tau * e_first(0:ND1,DDK,l) &
                           *         a_pts(0,j,DDK,l) &
                           *     JacNorVec(j,Edge_Num,DDK,l) &
                           *     NorVec_mg(j,Edge_Num,DDK,l) &
                           / (omega**2 * c_bar(j,Edge_Num,DDK,l))
                        enddo
                     endif
!-------------------------------------------------------------------------
                  case (2,3) ! Neumann and Robin
                     if (Edge_Num .eq. 3) then
                        do i=0,ND
                           tauND1(i,0:ND2,Edge_Num,DDK,l) = &
                              e_end(0:ND2,DDK,l) &
                           * JacNorVec(i,Edge_Num,DDK,l) &
                           *     a_pts(i,ND2,DDK,l) &
                           / omega
                        enddo
                     endif
      
                     if (Edge_Num .eq. 1) then
                        do i=0,ND
                           tauND1(i,0:ND2,Edge_Num,DDK,l) = &
                              e_first(0:ND2,DDK,l) &
                           * JacNorVec(i,Edge_Num,DDK,l) &
                           *     a_pts(i,0,DDK,l) &
                           /  omega
                        enddo
                     endif
      
                     if (Edge_Num .eq. 2) then
                        do j=0,ND
                           tauND1(0:ND1,j,Edge_Num,DDK,l) = &
                              e_end(0:ND1,DDK,l) &
                           * JacNorVec(j,Edge_Num,DDK,l) &
                           *     a_pts(ND1,j,DDK,l) &
                           / omega
                        enddo
                     endif

                     if (Edge_Num .eq. 4) then
                        do j=0,ND
                           tauND1(0:ND1,j,Edge_Num,DDK,l) = &
                              e_first(0:ND1,DDK,l) &
                           * JacNorVec(j,Edge_Num,DDK,l) &
                           *     a_pts(0,j,DDK,l) &
                           / omega
                        enddo
                     endif
!=============================================================================
                  case(0) ! Interface penalty
                     DDK_Connect  =  DM_Connect(1,Edge_Num,DDK)
                     Edge_Connect =  DM_Connect(2,Edge_Num,DDK)
                     ND1_Connect  = PolyDegN_DM(1,DDK_Connect,l)
                     ND2_Connect  = PolyDegN_DM(2,DDK_Connect,l)
      
                     select case(Edge_Connect)
                        case(2,4)
                           NDFix=0; if (Edge_Connect .eq. 2 ) NDFix = ND1_Connect
      
                           OmegaSurr=LGLWeights(0,ND1_Connect)

                           do j=0,ND2_Connect
                              SigmaSurr(j,Edge_Connect,DDK_Connect,l) = &
                              - (  Jacobin(NDFix,j,DDK_Connect,l) &
                              *  NorVec_mg(j,Edge_Connect,DDK_Connect,l) &
                              * LGLWeights(j,ND2_Connect)) &
                              / 2
      
                              SigmaSurr2(j,Edge_Connect,DDK_Connect,l) = &
                                 (  c_sigma(l) * LGLWeights(j,ND2_Connect) &
                                 *   Jacobin(NDFix,j,DDK_Connect,l) &
                                 *     a_pts(NDFix,j,DDK_Connect,l) &
                                 * NorVec_mg(j,Edge_Connect,DDK_Connect,l)**2) &
                                 / (4 * c_bar(j,Edge_Connect,DDK_Connect,l) &
                                 * OmegaSurr )
                           enddo
                        case(1,3)
                           NDFix=0; if (Edge_Connect .eq. 3 ) NDFix = ND2_Connect
      
                           OmegaSurr=LGLWeights(0,ND2_Connect)
                           
                           do i =0,ND1_Connect
                              SigmaSurr(i,Edge_Connect,DDK_Connect,l) = &
                              - (  Jacobin(i,NDFix,DDK_Connect,l) &
                              *  NorVec_mg(i,Edge_Connect,DDK_Connect,l) &
                              * LGLWeights(i,ND1_Connect) ) &
                              / 2
      
                              SigmaSurr2(i,Edge_Connect,DDK_Connect,l) = &
                              ( c_sigma(l) * LGLWeights(i,ND1_Connect) &
                              *   Jacobin(i,NDFix,DDK_Connect,l) &
                              *     a_pts(i,NDFix,DDK_Connect,l) &
                              * NorVec_mg(i,Edge_Connect,DDK_Connect,l)**2 ) &
                              / (4 * c_bar(i,Edge_Connect,DDK_Connect,l) &
                              * OmegaSurr )
                           enddo
                     end select !Edge_Connect

                     select case (Edge_Num)
                        case (2,4)
                           NDFix=0; if (Edge_Num .eq. 2 ) NDFix = ND1
                           OmegaSelf=LGLWeights(0,ND1)
      
                           do j=0,ND2
                              SigmaHat(j,Edge_Num,DDK,l) = &
                              (   Jacobin(NDFix,j,DDK,l) &
                              * NorVec_mg(j,Edge_Num,DDK,l) ) &
                              / (2 * OmegaSelf )

                              Sigma_tild(j,Edge_Num,DDK,l) = &
                              (c_sigma(l) * Jacobin(NDFix,j,DDK,l) &
                              *  NorVec_mg(j,Edge_Num,DDK,l)**2) &
                              *      a_pts(NDFix,j,DDK,l) &
                              / (4 * c_bar(j,Edge_Num,DDK,l) &
                              * OmegaSelf**2 )
      
                              tau1(0:ND1,j,Edge_Num,DDK,l) = &
                                 Diff_xi1(NDFix,0:ND1,ND1) &
                              *  Nor_dot_gradxi(j,Edge_Num,DDK,l) &
                              *           a_pts(NDFix,j,DDK,l) &
                              / ( LGLWeights(j,ND2) &
                              *   LGLWeights(0:ND1,ND1) )
      
                              tau2(j,Edge_Num,DDK,l) = &
                                 1 / (LGLWeights(NDFix,ND1) * LGLWeights(j,ND2))
      
                              tau3(j,0:ND2,Edge_Num,DDK,l) = &
                                 Diff_xi2(0:ND2,j,ND2) &
                              *  Nor_dot_gradeta(j,Edge_Num,DDK,l) &
                              *            a_pts(NDFix,j,DDK,l) &
                              / ( LGLWeights(0:ND2,ND2) &
                              *   LGLWeights(NDFix,ND1) )
                           enddo
!     pathch direction
                           if (DM_Connect(3,Edge_Num,DDK) .eq. 1) then
      
                              do j=0,ND2
                                 tau1(0:ND1,j,Edge_Num,DDK,l) = &
                                    SigmaSurr(j,Edge_Connect,DDK_Connect,l) &
                                 *      tau1(0:ND1,j,Edge_Num,DDK,l)
      
                                 tau2(j,Edge_Num,DDK,l) = &
                                    SigmaSurr2(j,Edge_Connect,DDK_Connect,l) &
                                 *       tau2(j,Edge_Num,DDK,l)

                                 tau3(j,0:ND2,Edge_Num,DDK,l) = &
                                    SigmaSurr(j,Edge_Connect,DDK_Connect,l) &
                                 *      tau3(j,0:ND2,Edge_Num,DDK,l)
                              enddo
                           else

                              do j=0,ND2
                                 tau1(0:ND1,j,Edge_Num,DDK,l) = &
                                    SigmaSurr(ND2-j,Edge_Connect,DDK_Connect,l) &
                                 *      tau1(0:ND1,j,Edge_Num,DDK,l)
      
                                 tau2(j,Edge_Num,DDK,l) = &
                                    SigmaSurr2(ND2-j,Edge_Connect,DDK_Connect,l) &
                                 *       tau2(j,Edge_Num,DDK,l)
      
                                 tau3(j,0:ND2,Edge_Num,DDK,l) = &
                                    SigmaSurr(ND2-j,Edge_Connect,DDK_Connect,l) &
                                 *      tau3(j,0:ND2,Edge_Num,DDK,l)
                              enddo
                           endif !Patch_direction

                        case(1,3)
                           NDFix=0; if (Edge_Num .eq. 3 ) NDFix = ND2
      
                           OmegaSelf=LGLWeights(0,ND2)
      
                           do i =0,ND1
                              SigmaHat(i,Edge_Num,DDK,l) = &
                              (   Jacobin(i,NDFix,DDK,l) &
                              * NorVec_mg(i,Edge_Num,DDK,l) ) &
                              / (2  * OmegaSelf )
      
                              Sigma_tild(i,Edge_Num,DDK,l) = &
                              (   c_sigma(l) * Jacobin(i,NDFix,DDK,l) &
                              *  NorVec_mg(i,Edge_Num,DDK,l)**2 ) &
                              *      a_pts(i,NDFix,DDK,l) &
                              / (4 * c_bar(i,Edge_Num,DDK,l) &
                              * OmegaSelf**2  )
      
                              tau1(i,0:ND2,Edge_Num,DDK,l) = &
                                 Diff_xi2(0:ND2,NDFix,ND2) &
                              * Nor_dot_gradeta(i,Edge_Num,DDK,l) &
                              *           a_pts(i,NDFix,DDK,l) &
                              / ( LGLWeights(i,ND1) &
                              * LGLWeights(0:ND2,ND2) )
      
                              tau2(i,Edge_Num,DDK,l) = &
                                 1 / (LGLWeights(NDFix,ND2) * LGLWeights(i,ND1))
      
                              tau3(i,0:ND1,Edge_Num,DDK,l) = &
                                 Diff_xi1(i,0:ND1,ND1) &
                              * Nor_dot_gradxi(i,Edge_Num,DDK,l) &
                              *      a_pts(i,NDFix,DDK,l) &
                              / ( LGLWeights(0:ND1,ND1) &
                              *  LGLWeights(NDFix,ND2))
                           enddo
!     pathch direction
                           if (DM_Connect(3,Edge_Num,DDK) .eq. 1) then
                              do i=0,ND1
                                 tau1(i,0:ND2,Edge_Num,DDK,l) = &
                                    SigmaSurr(i,Edge_Connect,DDK_Connect,l) &
                                 *      tau1(i,0:ND2,Edge_Num,DDK,l)
      
      
                                 tau2(i,Edge_Num,DDK,l) = &
                                    SigmaSurr2(i,Edge_Connect,DDK_Connect,l) &
                                 *       tau2(i,Edge_Num,DDK,l)
      
                                 tau3(i,0:ND1,Edge_Num,DDK,l) = &
                                    SigmaSurr(i,Edge_Connect,DDK_Connect,l) &
                                 *      tau3(i,0:ND1,Edge_Num,DDK,l)
                              enddo
                           else
                              do i=0,ND1
                                 tau1(i,0:ND2,Edge_Num,DDK,l) = &
                                    SigmaSurr(ND1-i,Edge_Connect,DDK_Connect,l) &
                                 *      tau1(i,0:ND2,Edge_Num,DDK,l)
      
                                 tau2(i,Edge_Num,DDK,l) = &
                                    SigmaSurr2(ND1-i,Edge_Connect,DDK_Connect,l) &
                                 *      tau2(i,Edge_Num,DDK,l)
      
                                 tau3(i,0:ND1,Edge_Num,DDK,l) = &
                                    SigmaSurr(ND1-i,Edge_Connect,DDK_Connect,l) &
                                 *      tau3(i,0:ND1,Edge_Num,DDK,l)
                              enddo
                           endif ! Patch_direction
                     end select ! Edge_Num
                  end select ! BC_Type
            enddo ! Edge_Num
         enddo ! DDK
      enddo ! l
      
      return
      
      end subroutine Init_Penalty_Parameters
!-------------------------------------------------------------
      subroutine Initial_Field(mgl)
      use constants
      use State_Var
      use MD2D_Grid
      implicit none
      
      integer:: demo_case,remainder
      real(kind=8):: time_eval
      real(kind=8):: v_init, dFdx_init, dFdy_init, F_init
      
      integer:: i,j,mgl 
      real(kind=8):: x_coord, y_coord

      do DDK=1,TotNum_DM 
         ND1=PolyDegN_DM(1,DDK,mgl); ND2=PolyDegN_DM(2,DDK,mgl)
         do j=0,ND2
            do i=0,ND1
               x_coord=x1(i,j,DDK,mgl); y_coord=x2(i,j,DDK,mgl)
               v(i,j,DDK)    = v_init(x_coord,y_coord,DDK)
               rhs(i,j,DDK)  = F_init(x_coord,y_coord,DDK) ! compute flux
               dFdx(i,j,DDK) = dFdx_init(x_coord,y_coord,DDK)
               dFdy(i,j,DDK) = dFdy_init(x_coord,y_coord,DDK)
            enddo ! i
         enddo ! j
      enddo !DDK

      open(1238,file='rhs_withoutmask.text')
      do DDK = 1, TotNum_DM
         ND1=PolyDegN_DM(1,DDK,mgl); ND2=PolyDegN_DM(2,DDK,mgl)
            do j = 0, ND2
               do i = 0, ND1
                  write(1238,*)(DDK-1)*(ND1+1)*(ND2+1)+j*(ND1+1)+i, &
                                rhs(i,j,DDK)
               enddo
            enddo
      enddo
      close(1238)

      call outpost(rhs,mgl,'rwo')
      return
      
      end subroutine Initial_Field
!------------------------------------------------------------------------------
      subroutine FluxComp(mgl)
      use constants
      use Legendre
      use MD2D_Grid
      use State_Var
      implicit none

      integer :: i, j, mgl

!     Compute PBC
      call Compute_PBC(mgl)
      
!     add BC
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,mgl); ND2=PolyDegN_DM(2,DDK,mgl)

         rhs(0:ND1,0:ND2,DDK) = rhs(0:ND1,0:ND2,DDK) * Jacobin(0:ND1,0:ND2,DDK,mgl)
      
!     Add second set penalty boundary condition
!     side 4
         Edge_Num=4
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               do j = 0,ND2
                  rhs(0:ND1,j,DDK) =       rhs(0:ND1,j,DDK) &
                                   + (    tauD(0:ND1,j,Edge_Num,DDK,mgl) &
                                   - tau_tilde(0:ND1,j,Edge_Num,DDK,mgl)) &
                                   *         G(j,Edge_Num,DDK)
                  rhs(0,0:ND2,DDK) =     rhs(0,0:ND2,DDK) &
                                   + ( tauD2(0:ND2,j,Edge_Num,DDK,mgl) &
                                   *       G(j,Edge_Num,DDK) )
               enddo
            case(2,3)
               do j = 0,ND2
                  rhs(0:ND1,j,DDK) =    rhs(0:ND1,j,DDK) &
                                   + tauND1(0:ND1,j,Edge_Num,DDK,mgl) &
                                     *    G(j,Edge_Num,DDK)
               enddo
         end select
!     side 2
         Edge_Num=2
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               do j=0,ND2
                  rhs(0:ND1,j,DDK) = rhs(0:ND1,j,DDK) &
                                   + (    tauD(0:ND1,j,Edge_Num,DDK,mgl) &
                                   - tau_tilde(0:ND1,j,Edge_Num,DDK,mgl)) &
                                   *         G(j,Edge_Num,DDK)
                  rhs(ND1,0:ND2,DDK) = rhs(ND1,0:ND2,DDK) &
                                     + (   tauD2(0:ND2,j,Edge_Num,DDK,mgl) &
                                     *         G(j,Edge_Num,DDK) )
               enddo
            case(2,3)
               do j=0,ND2
                  rhs(0:ND1,j,DDK) = rhs(0:ND1,j,DDK) &
                                   + tauND1(0:ND1,j,Edge_Num,DDK,mgl) &
                                   *      G(j,Edge_Num,DDK)
               enddo
         end select
!     side 1
         Edge_Num=1
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               do i=0,ND1
                  rhs(i,0:ND2,DDK) = rhs(i,0:ND2,DDK) &
                                   + (    tauD(i,0:ND2,Edge_Num,DDK,mgl) &
                                   - tau_tilde(i,0:ND2,Edge_Num,DDK,mgl)) &
                                   *         G(i,Edge_Num,DDK)
                  rhs(0:ND1,0,DDK) = rhs(0:ND1,0,DDK) &
                                   + ( tauD2(0:ND1,i,Edge_Num,DDK,mgl) &
                                   *       G(i,Edge_Num,DDK) )
               enddo
            case(2,3)
               do i=0,ND1
                  rhs(i,0:ND2,DDK) = rhs(i,0:ND2,DDK) &
                                   + tauND1(i,0:ND2,Edge_Num,DDK,mgl) &
                                   *      G(i,Edge_Num,DDK)
               enddo
         end select
!     side 3
         Edge_Num=3
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               do i=0,ND1
                  rhs(i,0:ND2,DDK) = rhs(i,0:ND2,DDK) &
                                   + (    tauD(i,0:ND2,Edge_Num,DDK,mgl) &
                                   - tau_tilde(i,0:ND2,Edge_Num,DDK,mgl)) &
                                   *         G(i,Edge_Num,DDK)
                  rhs(0:ND1,ND2,DDK) = rhs(0:ND1,ND2,DDK) &
                                     + ( tauD2(0:ND1,i,Edge_Num,DDK,mgl)&
                                     *       G(i,Edge_Num,DDK) )
               enddo
            case(2,3)
               do i=0,ND1
                  rhs(i,0:ND2,DDK) = rhs(i,0:ND2,DDK) &
                                   + tauND1(i,0:ND2,Edge_Num,DDK,mgl) &
                                   *      G(i,Edge_Num,DDK)
               enddo
         end select
      enddo !DDK
      
!     Multiply mass matrix in y direction
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,mgl); ND2=PolyDegN_DM(2,DDK,mgl)
         do j=0,ND2
            do i=0,ND1
               rhs(i,j,DDK) = rhs(i,j,DDK) * LGLWeights(j,ND2) * LGLWeights(i,ND1)
            enddo
         enddo
      enddo

!     Check gloabl rhs vector times Mass matrix
      open(1238,file='rhs.text')
      do DDK = 1, TotNum_DM
         ND1=PolyDegN_DM(1,DDK,mgl); ND2=PolyDegN_DM(2,DDK,mgl)
            do j = 0, ND2
               do i = 0, ND1
                  write(1238,*)(DDK-1)*(ND1+1)*(ND2+1)+j*(ND1+1)+i, &
                                rhs(i,j,DDK)
               enddo
            enddo
      enddo
      close(1238)
      call outpost(rhs,mgl,'rhs')
      
      return
      
      end subroutine FluxComp
!--------------------------------------------------------------------------------------
      subroutine compute_PBC(mgl)
      use MD2D_Grid
      use State_Var
      implicit none

      integer:: i,j, NDFix, mgl
      real(kind=8):: x_coord, y_coord
      real(kind=8):: nx, ny, alpha, beta, g_ext, b_loc

!     Comput boundary condition for each element and each edge
      
      do DDK=1,TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,mgl); ND2 = PolyDegN_DM(2,DDK,mgl)
         do Edge_Num = 1,4
            select case (Edge_Num)
               case (2,4)
                  ND=ND2; NDFix=0; if (Edge_Num .eq. 2 ) NDFix = ND1
              
                  select case (BC_Type(Edge_Num,DDK))
                     case(1,2,3)
                        do j=0,ND  !(ND=ND2)
                           x_coord= x1(NDFix,j,DDK,mgl);
                           y_coord= x2(NDFix,j,DDK,mgl);
                           nx     = NorVec_x1(j,Edge_Num,DDK,mgl)
                           ny     = NorVec_x2(j,Edge_Num,DDK,mgl)
                           alpha  = bnd_alpha(j,Edge_Num,DDK,mgl)
                           beta   =  bnd_beta(j,Edge_Num,DDK,mgl)
                           G(j,Edge_Num,DDK) = g_ext(x_coord,y_coord,&
                                             nx,ny,alpha,beta,b_loc,Edge_Num)
                        enddo
                  end select ! BC_Type
               case (1,3)
                  ND=ND1; NDFix=0; if (Edge_Num .eq. 3 ) NDFix = ND2
              
                  select case (BC_Type(Edge_Num,DDK)) 
                     case(1,2,3)     
                        do i=0,ND ! ND=ND1
                           x_coord= x1(i,NDFix,DDK,mgl);
                           y_coord= x2(i,NDFix,DDK,mgl);
                           nx     = NorVec_x1(i,Edge_Num,DDK,mgl)
                           ny     = NorVec_x2(i,Edge_Num,DDK,mgl)
                           alpha  = bnd_alpha(i,Edge_Num,DDK,mgl)
                           beta   =  bnd_beta(i,Edge_Num,DDK,mgl)
                           G(i,Edge_Num,DDK) = g_ext(x_coord,y_coord,&
                                             nx,ny,alpha,beta,b_loc,Edge_Num)
                        enddo !i
                  end select ! BC_Type
            end select !Edge_Num
         enddo ! Edge_Num
      enddo ! DDK
      
      end subroutine 
