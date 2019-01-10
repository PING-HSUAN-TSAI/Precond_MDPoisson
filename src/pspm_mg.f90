      subroutine hsmg_setup
      use MD2D_Grid
      use MultiGrid_Var
      implicit none

      call hsmg_setup_mg_nx !! set nx values for each level of multigrid
      call hsmg_setup_semhat
      call hsmg_setup_intp
!      call hsmg_do_wt
      call hsmg_setup_fast
      return
      end subroutine
!----------------------------------------------------------------------
      subroutine hsmg_setup_mg_nx
      use MD2D_Grid
      use Multigrid_Var
      implicit none

      integer :: i,l,n

      mg_lmax = 3

      if (mg_lmax == 3) then
      
         mg_nx(1) = 1
         mg_ny(1) = 1
         mg_nz(1) = 1

!      mg_nx(2) = 2
!      mg_ny(2) = 2
!      mg_nz(2) = 2
         mg_nx(2) = PolyDegN_DM(1,1,1)/2 
         mg_ny(2) = PolyDegN_DM(1,1,1)/2
         mg_nz(2) = PolyDegN_DM(1,1,1)/2
   
         mg_nx(3) = PolyDegN_DM(1,1,1) 
         mg_ny(3) = PolyDegN_DM(1,1,1)
         mg_nz(3) = PolyDegN_DM(1,1,1)

      else if (mg_lmax == 2) then

         mg_nx(1) = PolyDegN_DM(1,1,1)/2 
         mg_ny(1) = PolyDegN_DM(1,1,1)/2
         mg_nz(1) = PolyDegN_DM(1,1,1)/2
   
         mg_nx(2) = PolyDegN_DM(1,1,1) 
         mg_ny(2) = PolyDegN_DM(1,1,1)
         mg_nz(2) = PolyDegN_DM(1,1,1)

      endif

      write(10,*) 'mg_nx:',(mg_nx(i),i=1,mg_lmax)
      write(10,*) 'mg_ny:',(mg_ny(i),i=1,mg_lmax)
      write(10,*) 'mg_nz:',(mg_nz(i),i=1,mg_lmax)

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine hsmg_setup_semhat
      use Legendre
      use MD2D_Grid
      use Multigrid_Var
      implicit none 

      integer:: n, l  

      do l=1,mg_lmax
         n = mg_nx(l)     ! polynomial order

         mg_zh(0:n,l) = LGLCoord(0:n,n)
!         write(*,*)mg_zh(0:n,l)
         mg_nh(l)=n+1

      enddo

      return
      end subroutine
!----------------------------------------------------------------------
      subroutine hsmg_setup_intp
      use Multigrid_Var
      implicit none
      integer:: l,nf,nc
      integer:: i,j

      do l=1,mg_lmax-1

         nf=mg_nh(l+1)
         nc=mg_nh(l)
!         write(*,*)nf,nc
!         write(*,*)mg_zh(0:nf-1,l+1)

!        Standard multigrid coarse-to-fine interpolation

         call hsmg_setup_intpm(mg_jh(0:nf-1,0:nc-1,l),mg_zh(:,l+1),mg_zh(:,l),nf,nc)

         mg_jht(0:nc-1,0:nf-1,l) = transpose(mg_jh(:,:,l))

!        Fine-to-coarse interpolation for variable-coefficient operators
!         call hsmg_setup_intpm(mg_jhfc(1,1,1,l),mg_zh(1,1,l),mg_zh(1,1,l+1),nc,nf)
!         call transpose(mg_jhfct(1,l),nf,mg_jhfc(1,l),nc)

      enddo
      end
!!----------------------------------------------------------------------
      subroutine hsmg_setup_intpm(jh,zf,zc,nf,nc)
      use Legendre
      use MD2D_Grid
      use Multigrid_Var
      implicit none

      integer:: level, i, j, ND1c, ND2c
      integer:: nf,nc
      real(kind=8):: jh(0:nf-1,0:nc-1),zf(0:nf-1),zc(0:nc-1)!,zc(0:nc-1,1,1)
      real(kind=8):: w_mg(0:nc-1,1:2)

      w_mg = 0.d0

      do i=0,nf-1
         call fd_weights_full(zf(i),zc,nc-1,1,w_mg(0:nc-1,1:2))
         do j=0,nc-1
            jh(i,j)=w_mg(j,1)
!            write(*,*)i,j,jh(i,j)
         enddo
      enddo
      return
      end
!!----------------------------------------------------------------------
      subroutine hsmg_tnsr(v,nv,u,nu,A,At,LD1,LD2)
      use MD2D_Grid
!     computes
!     v = [A (x) A] u      or
!     v = [A (x) A (x) A] u
      integer:: nv,nu, LD1, LD2, i, j
      real(kind=8):: v(0:nv-1,0:nv-1,1:TotNum_DM)
      real(kind=8):: u(0:nu-1,0:nu-1,1:TotNum_DM)
      real(kind=8):: tmp(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8):: A(0:LD1,0:LD2)
      real(kind=8):: At(0:LD1,0:LD2)

      do DDK = 1, TotNum_DM
         tmp(0:nu-1,0:nv-1,DDK) = &
            Matmul(u(0:nu-1,0:nu-1,DDK), At(0:nu-1,0:nv-1))

         v(0:nv-1,0:nv-1,DDK) = &
            Matmul(A(0:nv-1,0:nu-1), tmp(0:nu-1,0:nv-1,DDK))
      enddo
      return
      end
!----------------------------------------------------------------------
      subroutine hsmg_intp(uf,uc,l,LD1,LD2) ! l is coarse level
      use MD2D_Grid
      use Multigrid_Var
      implicit none
      real(kind=8):: uf(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8):: uc(0:LD1,0:LD2,1:TotNum_DM)
      integer:: l,LD1,LD2
      call hsmg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(:,:,l),mg_jht(:,:,l),LD1,LD2)
      return
      end
!----------------------------------------------------------------------
      subroutine hsmg_do_wt
      use MD2D_Grid
      use Multigrid_Var
      use State_Var
      implicit none
      integer :: i, j, n, level 

      level = 1
      open(unit=30,file='scale_c')

      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
!     Side 4
         Edge_Num=4
         select case (BC_Type(Edge_Num,DDK))
            case(0)

            do j=0,ND2
               n=0+1+j*(ND2+1)
               scale_c(n,DDK) = scale_c(n,DDK) + 1
            enddo

         end select
      ! side 2
         Edge_Num=2
         select case (BC_Type(Edge_Num,DDK))
            case(0)

            do j=0,ND2
               n=ND1+1+j*(ND2+1)
               scale_c(n,DDK) = scale_c(n,DDK) + 1
            enddo
         end select
      ! side 1
         Edge_Num=1
         select case (BC_Type(Edge_Num,DDK))
            case(0)

            do i=0,ND1
               n=i+1+0*(ND2+1)
               scale_c(n,DDK) = scale_c(n,DDK) + 1
            enddo
         end select
         Edge_Num=3
         select case (BC_Type(Edge_Num,DDK))
            case(0)

            do i=0,ND1
               n=i+1+ND2*(ND2+1)
               scale_c(n,DDK) = scale_c(n,DDK) + 1
            enddo
         end select

         do j=0,ND2
            do i=0,ND1
               n=i+1+j*(ND2+1)
               write(30,*)i,j,n,scale_c(n,DDK)
            enddo
         enddo

      enddo
      return
      end
!!----------------------------------------------------------------------
      subroutine hsmg_setup_fast
      use Multigrid_Var
      use MD2D_Grid
      implicit none
      integer :: il, i, j, l

      call hsmg_setup_fast1d_a
      call hsmg_setup_fast1d_b
      call generalev

      do l=2,mg_lmax
         il = 1
   
         open(unit=49,file='DDiagonal.text')
         do DDK = 1, TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
            do j =0,ND2
               do i = 0,ND1
                  DDiagonal(il,DDK,l) = 1.0 / ( lamxx(i,DDK,l) + lamyy(j,DDK,l))
                  write(49,*)lamxx(i,DDK,l),lamyy(j,DDK,l)
                  write(49,*)i,j,DDK,il,DDiagonal(il,DDK,l)
                  il = il +1
               enddo
            enddo
         enddo
         close(49)
      enddo
      end
!!----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d_a
      use Legendre
      use MD2D_Grid
      use Multigrid_Var
      use State_Var
      implicit none
      integer :: i,j,l
      integer :: ND1p, ND2p
      real(kind=8) :: dudx(0:PolyDegN_DM(1,1,1),0:PolyDegN_DM(2,1,1),1:TotNum_DM) 
      real(kind=8) :: dudy(0:PolyDegN_DM(1,1,1),0:PolyDegN_DM(2,1,1),1:TotNum_DM) 

      l = 1

!     Initialize variables for Lx and Ly operator


!     Construct 1d Lx, Ly operator without adding penalties
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1

         dudx(0:ND1,0:ND1,DDK) = a_pts(0:ND1,0:ND2,DDK,l)  &
                               * (dxi1_dx1(0:ND1,0:ND2,DDK,l)) &
                               *  Diff_xi1(0:ND1,0:ND1,ND1)

         dudy(0:ND2,0:ND2,DDK) = a_pts(0:ND1,0:ND2,DDK,l)  &
                               * (dxi2_dx2(0:ND1,0:ND2,DDK,l)) &
                               *  Diff_xi1(0:ND1,0:ND1,ND1)
!
         Lxx(0:ND1,0:ND1,DDK) = - Matmul(Diff_xi1(0:ND1,0:ND1,ND1), &
                                             dudx(0:ND1,0:ND1,DDK) )

         Lyy(0:ND2,0:ND2,DDK) = - Matmul(Diff_xi1(0:ND1,0:ND1,ND1), &
                                             dudy(0:ND2,0:ND2,DDK) )
      enddo

      write(10,*)'L operator without penalty'
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         do j=0,ND2
            do i=0,ND1
               write(10,*)i,j,DDK,Lxx(i,j,DDK),Lyy(i,j,DDK)
            enddo
         enddo
      enddo
      write(10,*)'end'

!     Construct 1d penalty
      call hsmg_setup_penalty1d

!     Start adding penalty for every DDK
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
!     side 4
         Edge_Num=4
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
               do i=0,ND1
                  write(10,*)i,tau(i,Edge_Num,DDK,l),tautil(i,Edge_Num,DDK,l)
               enddo
               Lxx(0:ND1,0,DDK) = Lxx(0:ND1,0,DDK) &
                                + (   tau(0:ND1,Edge_Num,DDK,l)   &
                                -  tautil(0:ND1,Edge_Num,DDK,l) )
            case(2,3)
               do j = 0,ND2
                  Lxx(0:ND1,j,DDK) = Lxx(0:ND1,j,DDK) &
                                   + tauND1(0:ND1,j,Edge_Num,DDK,l) !&
               enddo
            case(0)
               write(10,*)'Edge:',Edge_Num,'interface penalty'

!     add interface penalty : first set of dirichlet
               Lxx(0:ND1,0,DDK) = Lxx(0:ND1,0,DDK) &
                                + sigmatil_1(0:ND1,Edge_Num,DDK,l)
               Lxx(0    ,0,DDK) = Lxx(0    ,0,DDK) &
                                + sigmatil_2(Edge_Num,DDK,l) !&
!     add interface penalty : second set of dirichlet
               Lxx(0    ,0,DDK) = Lxx(0    ,0,DDK) &
                                + Sigmabar(Edge_Num,DDK,l)  !&
!     add interface penalty : neumann
               Lxx(0,0:ND2,DDK) = Lxx(0,0:ND2,DDK) &
                                - ( SigmaHat1d(Edge_Num,DDK,l) &
                                *     Diff_xi1(0,0:ND1,ND1) ) &
                                / dx1_dxi1(0,0,DDK,l) 

               write(10,*)'Edge:',Edge_Num,'interface penalty'
               do i=0,ND1
                  write(10,*)i,sigmatil_1(i,Edge_Num,DDK,l)
               enddo
               write(10,*)sigmatil_2(Edge_Num,DDK,l),Sigmabar(Edge_Num,DDK,l)&
                  ,SigmaHat1d(Edge_Num,DDK,l)
         end select

         write(10,*)'Lx operator after adding left penalty'
         do j=0,ND1
            do i =0,ND1
               write(10,*)i,j,Lxx(i,j,DDK)
            enddo
         enddo

!     side 2
         Edge_Num=2
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
               do i=0,ND1
                  write(10,*)i,tau(i,Edge_Num,DDK,l),tautil(i,Edge_Num,DDK,l)
               enddo

               Lxx(0:ND1,ND2,DDK) = Lxx(0:ND1,ND2,DDK) &
                                  + (  tau(0:ND1,Edge_Num,DDK,l) &
                                  - tautil(0:ND1,Edge_Num,DDK,l) )
            case(2,3)
               do j=0,ND2
                  Lxx(0:ND1,j,DDK) = Lxx(0:ND1,j,DDK) &
                                   + tauND1(0:ND1,j,Edge_Num,DDK,l) !&
               enddo
            case(0)
!     add interface penalty : first set of dirichlet
               Lxx(0:ND1,ND1,DDK) = Lxx(0:ND1,ND1,DDK) &
                                  + sigmatil_1(0:ND1,Edge_Num,DDK,l) !
               Lxx(  ND1,ND1,DDK) = Lxx(  ND1,ND1,DDK) &
                                  + sigmatil_2(Edge_Num,DDK,l)  !&
!     add interface penalty : second set of dirichlet
               Lxx(  ND1,ND1,DDK) = Lxx(  ND1,ND1,DDK) &
                                  + Sigmabar(Edge_Num,DDK,l)  !&
!     add interface penalty : neumann
               Lxx(ND1,0:ND2,DDK) = Lxx(ND1,0:ND2,DDK) &
                                    + ( SigmaHat1d(Edge_Num,DDK,l) &
                                    *     Diff_xi1(ND1,0:ND1,ND1) ) &
                                    / dx1_dxi1(ND1,0,DDK,l) 
               write(10,*)'Edge:',Edge_Num,'interface penalty'
               do i=0,ND1
                  write(10,*)i,sigmatil_1(i,Edge_Num,DDK,l)
               enddo
               write(10,*)sigmatil_2(Edge_Num,DDK,l),Sigmabar(Edge_Num,DDK,l)&
                  ,SigmaHat1d(Edge_Num,DDK,l)
         end select

         write(10,*)'Lx operator after adding penalty'
         do j=0,ND1
            do i =0,ND1
               write(10,*)i,j,Lxx(i,j,DDK)
            enddo
         enddo
!     side 1
         Edge_Num=1
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
               do j=0,ND2
                  write(10,*)i,tau(j,Edge_Num,DDK,l),tautil(j,Edge_Num,DDK,l)
               enddo

               Lyy(0:ND1,0,DDK) = Lyy(0:ND1,0,DDK) &
                                + (  tau(0:ND2,Edge_Num,DDK,l) &
                                - tautil(0:ND2,Edge_Num,DDK,l) )
            case(2,3)
               do i=0,ND1
                  Lyy(i,0:ND2,DDK) = Lyy(i,0:ND2,DDK) &
                                  +    tauND1(i,0:ND2,Edge_Num,DDK,l) !&
               enddo
            case(0)
!     Add interface penalty : first set of dirichlet
               Lyy(0:ND1,0,DDK) = Lyy(0:ND1,0,DDK) &
                                + sigmatil_1(0:ND2,Edge_Num,DDK,l)
               Lyy(0    ,0,DDK) = Lyy(0    ,0,DDK) &
                                + sigmatil_2(Edge_Num,DDK,l) !&
!     Add interface penalty : second set of dirichlet
               Lyy(0    ,0,DDK) = Lyy(0    ,0,DDK) &
                                + Sigmabar(Edge_Num,DDK,l)  !&
!     Add interface penalty : neumann
               Lyy(0,0:ND2,DDK) = Lyy(0,0:ND2,DDK) &
                                - ( SigmaHat1d(Edge_Num,DDK,l) &
                                *     Diff_xi1(0,0:ND2,ND2) ) &
                                / dx2_dxi2(0,0,DDK,l) 
         end select
         write(10,*)'Ly operator after adding left penalty'
         do j=0,ND1
            do i =0,ND1
               write(10,*)i,j,Lyy(i,j,DDK)
            enddo
         enddo
!     side 3
         Edge_Num=3
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
               do j=0,ND2
                  write(10,*)i,tau(j,Edge_Num,DDK,l),tautil(j,Edge_Num,DDK,l)
               enddo

               Lyy(0:ND1,ND2,DDK) = Lyy(0:ND1,ND2,DDK) &
                                  + (  tau(0:ND2,Edge_Num,DDK,l)  &
                                  - tautil(0:ND2,Edge_Num,DDK,l) )  
            case(2,3)
               do i=0,ND1
                  Lyy(i,0:ND2,DDK) = Lyy(i,0:ND2,DDK) &
                                   + tauND1(i,0:ND2,Edge_Num,DDK,l) !&
               enddo
            case(0)
!     Add interface penalty : first set of dirichlet
               Lyy(0:ND1,ND2,DDK) = Lyy(0:ND1,ND2,DDK) &
                                  + sigmatil_1(0:ND2,Edge_Num,DDK,l)
               Lyy(  ND1,ND2,DDK) = Lyy(  ND1,ND2,DDK) &
                                  + sigmatil_2(Edge_Num,DDK,l) !&
!     Add interface penalty : second set of dirichlet
               Lyy(  ND1,ND2,DDK) = Lyy(  ND1,ND2,DDK) &
                                  + Sigmabar(Edge_Num,DDK,l) !&
!     Add interface penalty : neumann
               Lyy(ND1,0:ND2,DDK) = Lyy(ND1,0:ND2,DDK) &
                                  + ( SigmaHat1d(Edge_Num,DDK,l) &
                                  *     Diff_xi1(ND2,0:ND2,ND2) ) &
                                  / dx2_dxi2(ND2,ND2,DDK,l)
         end select
         write(10,*)'Ly operator after adding penalty'
         do j=0,ND1
            do i =0,ND1
               write(10,*)i,j,Lyy(i,j,DDK)
            enddo
         enddo
      enddo !DDK

      write(10,*)'Lx,Ly after Mass'
      do DDK=1,TotNum_DM
         write(10,*)'DDK:',DDK
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do j = 0, ND2
            do i = 0, ND1
               Lxx(i,j,DDK) = Lxx(i,j,DDK) * LGLWeights(i,ND1)
               Lyy(i,j,DDK) = Lyy(i,j,DDK) * LGLWeights(i,ND2)
               write(10,*)i,j,Lxx(i,j,DDK),Lyy(i,j,DDK)
            enddo
         enddo
      enddo
      end
!!----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d_b
      use Multigrid_Var
      use Legendre
      use MD2D_Grid
      implicit none
      integer :: i,j,l
      
      do l=2,mg_lmax
      
!     Construct Matrix Bx and By which will be used in solving the generalized eigenvalue problem Ax = lambda B x
         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
            do i = 0, ND1
               Bxx(i,i,DDK,l) = LGLWeights(i,ND1)
            enddo
            do j= 0, ND2
               Byy(j,j,DDK,l) = LGLWeights(j,ND2)
            enddo
         enddo
      enddo

      end
!!----------------------------------------------------------------------
      subroutine generalev
      use Multigrid_Var
      use Legendre
      use MD2D_Grid
      implicit none
      integer :: i,j,lbw,info,l

      do l=2,mg_lmax

         do DDK = 1, TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l)
            lbw =  (ND1+1) * (ND1+1)
            call dsygv(1,'V','U',ND1+1,Lxx(0:ND1,0:ND1,DDK,l),&
                  ND1+1,Bxx(0:ND1,0:ND1,DDK,l),ND1+1,lamxx(0:ND1,DDK,l)&
                  ,bwxx(0:4*(ND1+1)*(ND2+1),DDK,l),lbw,info)
         enddo
         write(10,*)'Complete diagonalize Lx operator'

!     Storing the eigenvectors that are not normalized yet
         do DDK = 1, TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l)
            Sxx_t(0:ND1,0:ND1,DDK,l) = transpose(Lxx(0:ND1,0:ND1,DDK,l))
            Sxx(0:ND1,0:ND1,DDK,l) = Lxx(0:ND1,0:ND1,DDK,l)
         enddo
!     Start Diagonalize Ly operator
         do DDK = 1, TotNum_DM
            ND2 = PolyDegN_DM(2,DDK,l)
            lbw = 4 * (ND2+1) * (ND2+1)
            call dsygv(1,'V','U',ND2+1,Lyy(0:ND2,0:ND2,DDK,l),&
                  ND2+1,Byy(0:ND2,0:ND2,DDK,l),ND2+1,lamyy(0:ND2,DDK,l)&
                  ,bwyy(0:4*(ND1+1)*(ND2+1),DDK,l),lbw,info)
         enddo
         write(10,*)'Complete diagonalize Ly operator'

!     Storing the eigenvectors that are not normalized yet
         do DDK = 1, TotNum_DM
            ND2 = PolyDegN_DM(2,DDK,l)
            Syy_t(0:ND2,0:ND2,DDK,l) = transpose(Lyy(0:ND2,0:ND2,DDK,l))
            Syy(0:ND2,0:ND2,DDK,l) = Lyy(0:ND2,0:ND2,DDK,l)
         enddo

      enddo

!      write(*,*)'Output the normalizd eigenvectors of L_x and L_y'
!      do DDK = 1, TotNum_DM
!         ND1 = PolyDegN_DM(1,DDK,l)
!         do j=0,ND1
!            do i =0,ND1
!               write(10,*)i,j,Sxx(i,j,DDK),Syy(i,j,DDK)
!            enddo
!         enddo
!      enddo
      end
!!----------------------------------------------------------------------
      subroutine hsmg_setup_penalty1d
      use Legendre
      use MD2D_Grid
      use Multigrid_Var
      use State_Var
      implicit none
      integer :: i,j,l
      integer :: ND1p, ND2p, NDFix
      real(kind=8):: omega, denum, OmegaSurr,OmegaSelf

      l = 1

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
!======================================================================
               case(1) ! Dirichlet
                  if (Edge_Num .eq. 3) then
!     Penalty for orthogonal part
                     tau(0:ND2,Edge_Num,DDK,l) = &
                        - Diff_xi2(0:ND2,ND2,ND2) &
                        *    a_pts(ND1,ND2,DDK,l) &
                        / ( LGLWeights(0:ND2,ND2) &
                        * dx2_dxi2(ND1,0:ND2,DDK,l))
!     Penalty for ensuring positive definite
                     tautil(0:ND2,Edge_Num,DDK,l) = &
                        - c_tau * e_end(0:ND2,DDK,l) &
                        *       a_pts(ND1,ND2,DDK,l) &
                        / (omega**2 * dx2_dxi2(ND1,0:ND2,DDK,l))

                     do j=0,ND2
                        write(10,*)DDK,Edge_Num,j,tau(j,Edge_Num,DDK,l),tautil(j,Edge_Num,DDK,l)
                     enddo
                  endif

                  if (Edge_Num .eq. 1) then
!     Penalty for orthogonal part
                     tau(0:ND2,Edge_Num,DDK,l) = &
                        Diff_xi2(0:ND2,0,ND2) &
                        *      a_pts(0,0,DDK,l) &
                        / ( LGLWeights(0:ND2,ND2) &
                        * dx2_dxi2(0,0:ND2,DDK,l))
!     Penalty for ensuring positive definite
                     tautil(0:ND2,Edge_Num,DDK,l) = &
                        - c_tau * e_first(0:ND2,DDK,l) &
                        *         a_pts(0,0,DDK,l) &
                        / (omega**2 * dx2_dxi2(0,0:ND2,DDK,l))

                        do j=0,ND2
                           write(10,*)DDK,Edge_Num,j,tau(j,Edge_Num,DDK,l),tautil(j,Edge_Num,DDK,l)
                        enddo
                  endif

                  if (Edge_Num .eq. 2) then
!     Penalty for orthogonal part
                     tau(0:ND1,Edge_Num,DDK,l) = &
                        -  Diff_xi1(ND1,0:ND1,ND1) &
                        *      a_pts(ND1,ND2,DDK,l) &
                        / ( LGLWeights(0:ND1,ND1) &
                        *   dx1_dxi1(0:ND1,ND2,DDK,l))
!     Penatly for ensuring positive definite
                     tautil(0:ND1,Edge_Num,DDK,l) = &
                        - c_tau * e_end(0:ND1,DDK,l) &
                        *     a_pts(ND1,ND2,DDK,l) &
                        / (omega**2 * dx1_dxi1(0:ND1,ND2,DDK,l))

                        do j=0,ND2
                           write(10,*)DDK,Edge_Num,j,tau(j,Edge_Num,DDK,l),tautil(j,Edge_Num,DDK,l)
                        enddo
                  endif

                  if (Edge_Num .eq. 4) then
!     Penalty for orthogonal part
                     tau(0:ND1,Edge_Num,DDK,l) = &
                        Diff_xi1(0,0:ND1,ND1) &
                        *      a_pts(0,0,DDK,l) &
                        / ( LGLWeights(0:ND1,ND1) &
                        *   dx1_dxi1(0:ND1,0,DDK,l))
!     Penalty for ensuring positive definite
                     tautil(0:ND1,Edge_Num,DDK,l) = &
                        - c_tau * e_first(0:ND1,DDK,l) &
                        *         a_pts(0,0,DDK,l) &
                        / (omega**2 * dx1_dxi1(0:ND1,0,DDK,l))

                     do j=0,ND2
                        write(10,*)DDK,Edge_Num,j,tau(j,Edge_Num,DDK,l),tautil(j,Edge_Num,DDK,l)
                     enddo
                  endif
!!=========================================================================
!                 case (2,3) ! Neumann and Robin
!                     if (Edge_Num .eq. 3) then
!                        do i=0,ND
!                           tauND1(i,0:ND2,Edge_Num,DDK,l) = &
!                              e_end(0:ND2,DDK,l) &
!                           * JacNorVec(i,Edge_Num,DDK,l) &
!                           *     a_pts(i,ND2,DDK,l) &
!                           / omega
!                        enddo
!                     endif
!
!                     if (Edge_Num .eq. 1) then
!                        do i=0,ND
!                           tauND1(i,0:ND2,Edge_Num,DDK,l) = &
!                              e_first(0:ND2,DDK,l) &
!                           * JacNorVec(i,Edge_Num,DDK,l) &
!                           *     a_pts(i,0,DDK,l) &
!                           /  omega
!                        enddo
!                     endif
!
!                     if (Edge_Num .eq. 2) then
!                        do j=0,ND
!                           tauND1(0:ND1,j,Edge_Num,DDK,l) = &
!                              e_end(0:ND1,DDK,l) &
!                           * JacNorVec(j,Edge_Num,DDK,l) &
!                           *     a_pts(ND1,j,DDK,l) &
!                           / omega
!                        enddo
!                     endif
!
!                     if (Edge_Num .eq. 4) then
!                        do j=0,ND
!                           tauND1(0:ND1,j,Edge_Num,DDK,l) = &
!                              e_first(0:ND1,DDK,l) &
!                           * JacNorVec(j,Edge_Num,DDK,l) &
!                           *     a_pts(0,j,DDK,l) &
!                           / omega
!                        enddo
!                     endif
!!=============================================================================
               case(0) ! Interface penalty
                  DDK_Connect  =  DM_Connect(1,Edge_Num,DDK)
                  Edge_Connect =  DM_Connect(2,Edge_Num,DDK)
                  ND1_Connect  = PolyDegN_DM(1,DDK_Connect,l)
                  ND2_Connect  = PolyDegN_DM(2,DDK_Connect,l)

                  select case(Edge_Connect)
                     case(2,4)
                        NDFix=0; if (Edge_Connect .eq. 2 ) NDFix = ND1_Connect
                        OmegaSurr=LGLWeights(0,ND1_Connect)

                        SigmaSurr2_1d(Edge_Connect,DDK_Connect,l) = &
                           ( c_sigma(l) * a_pts(NDFix,0,DDK_Connect,l) )&
                           / (4 * dx1_dxi1(NDFix,0,DDK_Connect,l) * OmegaSurr )
                        write(10,*)'NDFix',NDFix
                        write(10,*)'SigmaSurr2_1d:','a_pts',a_pts(NDFix,0,DDK_Connect,l),&
                        c_sigma(l),OmegaSurr,SigmaSurr2_1d(Edge_Connect,DDK_connect,l)

                     case(1,3)
                        NDFix=0; if (Edge_Connect .eq. 3 ) NDFix = ND2_Connect
                        OmegaSurr=LGLWeights(0,ND2_Connect)

                        SigmaSurr2_1d(Edge_Connect,DDK_Connect,l) = &
                           ( c_sigma(l) * a_pts(0,NDFix,DDK_Connect,l) )&
                           / (4 * dx2_dxi2(0,NDFix,DDK_Connect,l) * OmegaSurr )

                  end select !Edge_Connect

                  select case (Edge_Num)
                     case (2,4)
                        NDFix=0; if (Edge_Num .eq. 2 ) NDFix = ND1
                        OmegaSelf=LGLWeights(0,ND1)

                        SigmaHat1d(Edge_Num,DDK,l) = 1 / ( 2 * OmegaSelf )

                        Sigmabar(Edge_Num,DDK,l) = &
                           ( c_sigma(l) * a_pts(NDFix,0,DDK,l) )&
                           / (4 * dx1_dxi1(NDFix,0,DDK,l) * OmegaSelf**2 )

                        sigmatil_1(0:ND1,Edge_Num,DDK,l) = &
                           Diff_xi1(NDFix,0:ND1,ND1) &
                           *  NorVec_x1(NDFix,Edge_Connect,DDK_Connect,l) &
                           *           a_pts(NDFix,0,DDK,l) &
                           / ( 2 * dx1_dxi1(NDFix,0,DDK,l) &
                           *   LGLWeights(0:ND1,ND1) )

                        sigmatil_2(Edge_Num,DDK,l) = 1 / LGLWeights(NDFix,ND1) 

!     pathch direction
                        if (DM_Connect(3,Edge_Num,DDK) .eq. 1) then

                           sigmatil_2(Edge_Num,DDK,l) = &
                              SigmaSurr2_1d(Edge_Connect,DDK_Connect,l) &
                              *       sigmatil_2(Edge_Num,DDK,l)
                        endif !Patch_direction

                        write(10,*)DDK,Edge_Num
                        write(10,*)SigmaHat1d(Edge_Num,DDK,l),Sigmabar(Edge_Num,DDK,l),&
                           sigmatil_1(0:ND1,Edge_Num,DDK,l),sigmatil_2(Edge_Num,DDK,l),SigmaSurr2_1d(Edge_Connect,DDK_Connect,l)

                     case(1,3)
                        NDFix=0; if (Edge_Num .eq. 3 ) NDFix = ND2
                        OmegaSelf=LGLWeights(0,ND2)

                        SigmaHat1d(Edge_Num,DDK,l) = 1 / ( 2 * OmegaSelf )

                        Sigmabar(Edge_Num,DDK,l) = &
                           (   c_sigma(l) *  a_pts(0,NDFix,DDK,l) ) &
                           / (4 * dx2_dxi2(0,NDFix,DDK,l) &
                           * OmegaSelf**2  )

                        sigmatil_1(0:ND2,Edge_Num,DDK,l) = &
                           Diff_xi2(0:ND2,NDFix,ND2) &
                           *  NorVec_x2(NDFix,Edge_Connect,DDK_Connect,l) &
                           *           a_pts(0,NDFix,DDK,l) &
                           / ( 2* dx2_dxi2(0,NDFix,DDK,l) &
                           * LGLWeights(0:ND2,ND2) )

                        sigmatil_2(Edge_Num,DDK,l) = 1 / LGLWeights(NDFix,ND2) 
!     pathch direction
                        if (DM_Connect(3,Edge_Num,DDK) .eq. 1) then

                           sigmatil_2(Edge_Num,DDK,l) = &
                              SigmaSurr2_1d(Edge_Connect,DDK_Connect,l) &
                              *       sigmatil_2(Edge_Num,DDK,l)
                        endif ! Patch_direction

                        write(10,*)DDK,Edge_Num
                        write(10,*)SigmaHat1d(Edge_Num,DDK,l),Sigmabar(Edge_Num,DDK,l),&
                           sigmatil_1(0:ND2,Edge_Num,DDK,l),sigmatil_2(Edge_Num,DDK,l)
                     end select ! Edge_Num
                  end select ! BC_Type
            enddo ! Edge_Num
         enddo ! DDK
!
      return
!
      end

      !================================================================
!      subroutine hsmg_schwarz(LD1,LD2,e,r,smoothpar,l)
!      use Legendre
!      use MD2D_Grid
!      use State_Var
!      use CG_Var
!      use Multigrid_Var
!      implicit none
!
!      integer :: i, j, LD1,LD2, vcycle, Nx, Ny
!      integer :: m_vcycle, m_smooth, TotN, TotNc
!      integer :: k, n, ND1p, ND2p,l, ND1c, ND2c, ind_JS
!      real(kind=8) :: ri(0:LD1,0:LD2,1:TotNum_DM)
!     real(kind=8) :: r(0:LD1,0:LD2,1:TotNum_DM)
!      real(kind=8) :: e(0:LD1,0:LD2,1:TotNum_DM)
!      real(kind=8) :: x_vcp(0:LD1,0:LD2,1:TotNum_DM)
!      real(kind=8):: smoothpar, err_sm, r_sm, err_vc, max_ef,glamax
!!     temporary variable for fdm
!      real(kind=8) :: tmp(0:LD1,0:LD2,1:TotNum_DM),tmp1(0:LD1,0:LD2,1:TotNum_DM)
!      real(kind=8) :: tmp2(0:LD1,0:LD2,1:TotNum_DM),tmp3(0:LD1,0:LD2,1:TotNum_DM)
!      real(kind=8) :: tmp4(0:LD1,0:LD2,1:TotNum_DM),tmp5(0:LD1,0:LD2,1:TotNum_DM)
!
!      Nx = PolyDegN_DM(1,1,l); Ny = PolyDegN_DM(2,1,l)
!      write(10,*)'hsmg_schwarz:',Nx,Ny
!
!      call copy(r,ri,Nx,Ny,l)
!
!      do DDK = 1, TotNum_DM
!         ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
!
!         tmp2(0:ND1,0:ND2,DDK) = &
!            Matmul(r(0:ND1,0:ND2,DDK), Syy(0:ND2,0:ND2,DDK,l))
!
!         tmp3(0:ND1,0:ND2,DDK) = &
!            Matmul(Sxx_t(0:ND1,0:ND1,DDK,l), tmp2(0:ND1,0:ND2,DDK))
!      enddo
!
!      do DDK = 1, TotNum_DM
!         ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
!         do j=0,ND2
!            do i=0,ND1
!               n = (DDK-1)*(ND1+1)*(ND2+1) + i + j*(ND2+1) +1
!!               n = i + j*(ND2+1) +1
!               tmp4(i,j,DDK) = DDiagonal(n,DDK,l) * tmp3(i,j,DDK)
!            enddo
!         enddo
!      enddo
!
!      do DDK = 1, TotNum_DM
!         ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
!         tmp5(0:ND1,0:ND2,DDK) = &
!            Matmul(tmp4(0:ND1,0:ND2,DDK), Syy_t(0:ND2,0:ND2,DDK,l))
!
!            e(0:ND1,0:ND2,DDK) = &
!            Matmul(Sxx(0:ND1,0:ND1,DDK,l), tmp5(0:ND1,0:ND2,DDK))
!      enddo
!
!
!            do DDK = 1, TotNum_DM
!               ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
!               do j=0,ND2
!                  do i=0,ND1
!                     e(i,j,DDK,l) = e(i,j,DDK,l) * smoothpar
!                  enddo
!               enddo
!            enddo
!      return
!      call add2s2(x_vcp,e(0:Nx,0:Ny,1:TotNum_DM,l),smoothpar,Nx,Ny,l)
!
!      call axhelm2(Nx,Ny,x_vcp,Ax_vc,l)
!
!      call add2s2(r,Ax_vc,-1.0,Nx,Ny,l)
!      r_sm = glamax(r,Nx,Ny,l)
!
!      end
