      module Multigrid_Var

      implicit none
      
!     Allocate variables for ML operator
      real(kind=8), save, allocatable :: I_U(:,:,:,:)
      real(kind=8), save, allocatable :: ML(:,:,:,:)
      
!     Smoothing variables
      real(kind=8), save, allocatable :: r_smooth(:,:,:)
      real(kind=8), save, allocatable :: rc_smooth(:,:,:)
      real(kind=8), save, allocatable :: b_smooth(:,:,:)
      real(kind=8), save, allocatable :: x_vc(:,:,:)
      real(kind=8), save, allocatable :: Ax_vc(:,:,:)
      real(kind=8), save, allocatable :: M(:) ! for Jacobi smoothing
      real(kind=8), save, allocatable :: Mr(:,:,:)
      
!     Interpolate matrix
      real(kind=8), save, allocatable :: xo(:,:), yo(:,:)
      real(kind=8), save, allocatable :: xi(:,:), yi(:,:)
      real(kind=8), save, allocatable :: Ihx(:,:,:)
      real(kind=8), save, allocatable :: Ihy(:,:,:)
      real(kind=8), save, allocatable :: Ihx_transpose(:,:,:)
      real(kind=8), save, allocatable :: Ihy_transpose(:,:,:)
      real(kind=8), save, allocatable :: wx(:,:,:)
      real(kind=8), save, allocatable :: wy(:,:,:)
      
      real(kind=8), save, allocatable :: test(:,:,:)
      
      !-coarse correction variables
      real(kind=8), save, allocatable :: xc_in(:,:,:)
      real(kind=8), save, allocatable :: ec(:,:,:)
      real(kind=8), save, allocatable :: ef(:,:,:)
      
      !-Lx, Ly operator
      real(kind=8), save, allocatable :: Lx(:,:,:)
      real(kind=8), save, allocatable :: Ly(:,:,:)
      real(kind=8), save, allocatable :: Bx(:,:,:)
      real(kind=8), save, allocatable :: By(:,:,:)
      real(kind=8), save, allocatable :: dudx_tmp(:,:,:)
      real(kind=8), save, allocatable :: dudy_tmp(:,:,:)
      real(kind=8), save, allocatable :: lamx(:,:)
      real(kind=8), save, allocatable :: lamy(:,:)
      real(kind=8), save, allocatable :: Lxx(:,:,:,:)
      real(kind=8), save, allocatable :: Lyy(:,:,:,:)
      real(kind=8), save, allocatable :: Bxx(:,:,:,:)
      real(kind=8), save, allocatable :: Byy(:,:,:,:)
      real(kind=8), save, allocatable :: lamxx(:,:,:)
      real(kind=8), save, allocatable :: lamyy(:,:,:)
      real(kind=8), save, allocatable :: bwxx(:,:,:)
      real(kind=8), save, allocatable :: bwyy(:,:,:)
      real(kind=8), save, allocatable :: Sxx_t(:,:,:,:)
      real(kind=8), save, allocatable :: Sxx(:,:,:,:)
      real(kind=8), save, allocatable :: Syy_t(:,:,:,:)
      real(kind=8), save, allocatable :: Syy(:,:,:,:)
      
      real(kind=8), save, allocatable :: z_ov(:,:,:)
      real(kind=8), save, allocatable :: z_ov_sol(:,:,:)
      real(kind=8), save, allocatable :: bwx(:,:)
      real(kind=8), save, allocatable :: bwy(:,:)
      real(kind=8), save, allocatable :: Diagonal(:,:)
      real(kind=8), save, allocatable :: DDiagonal(:,:,:)
      real(kind=8), save, allocatable :: Sx_t(:,:,:)
      real(kind=8), save, allocatable :: Sx(:,:,:)
      real(kind=8), save, allocatable :: Sy_t(:,:,:)
      real(kind=8), save, allocatable :: Sy(:,:,:)
      real(kind=8), save, allocatable :: Sx_t_norm(:,:,:)
      real(kind=8), save, allocatable :: Sx_norm(:,:,:)
      real(kind=8), save, allocatable :: Sy_t_norm(:,:,:)
      real(kind=8), save, allocatable :: Sy_norm(:,:,:)
      real(kind=8), save, allocatable :: BSx(:,:,:)
      real(kind=8), save, allocatable :: BSy(:,:,:)
      real(kind=8), save, allocatable :: SxBSx(:,:,:)
      real(kind=8), save, allocatable :: SyBSy(:,:,:)
      real(kind=8), save, allocatable :: scale_c(:,:)
      
      !-Wrapper variables
      real(kind=8), save, allocatable :: b_wrapper(:,:,:)
      real(kind=8), save, allocatable :: x_precond(:,:,:)
      real(kind=8), save, allocatable :: Ax_precond(:,:,:)
      real(kind=8), save, allocatable :: r_wrapper(:,:,:)
      real(kind=8), save, allocatable :: r_tmp(:,:,:)
      real(kind=8), save, allocatable :: pc_proAp(:)
      real(kind=8), save, allocatable :: Proj_AP(:,:,:,:)
      real(kind=8), save, allocatable :: Proj_P(:,:,:,:)
      real(kind=8), save, allocatable :: p_cond(:,:,:)
      real(kind=8), save, allocatable :: w(:,:,:)
      real(kind=8), save, allocatable :: Pap_pcond(:,:,:)
      real(kind=8):: pAp_wrapper
      real(kind=8):: pp
      real(kind=8):: alpha_wrapper
      integer :: Nk


      !mghs
      real(kind=8), save, allocatable :: mg_jh(:,:,:)
      real(kind=8), save, allocatable :: mg_jhfc(:,:,:)
      real(kind=8), save, allocatable :: mg_jht(:,:,:)
      real(kind=8), save, allocatable :: mg_jhfct(:,:,:)
      real(kind=8), save, allocatable :: mg_zh(:,:)
      integer :: mg_lmax, mg_nx(1:3), mg_ny(1:3), mg_nz(1:3)
      integer :: mg_nh(1:3)

      contains
!----------------------------------------------------------------------      
      subroutine alloc_mem_wrapper_var(DegMax,DegNcMax,TotNumDomain,step)

      implicit none
      integer:: DegMax(2)
      integer:: TotNumDomain
      integer:: ND1, ND2, DegNcMax
      integer:: ND1p, ND2p
      integer:: ierr,step
      
      ND1 = DegMax(1); ND2 = DegMax(2);
      ND1p = ND1 + 1; ND2p = ND2 + 1;
      allocate(   b_wrapper(0:ND1,0:ND2,1:TotNumDomain), &
                  x_precond(0:ND1,0:ND2,1:TotNumDomain), &
                 Ax_precond(0:ND1,0:ND2,1:TotNumDomain), &
                  r_wrapper(0:ND1,0:ND2,1:TotNumDomain), &
                      r_tmp(0:ND1,0:ND2,1:TotNumDomain), &
                    Proj_AP(0:ND1,0:ND2,1:TotNumDomain,1:step), &
                     Proj_P(0:ND1,0:ND2,1:TotNumDomain,1:step), &
                     p_cond(0:ND1,0:ND2,1:TotNumDomain), &
                          w(0:ND1,0:ND2,1:TotNumDomain), &
                  Pap_pcond(0:ND1,0:ND2,1:TotNumDomain), &
                   pc_proAp(1:step), &
                         ef(0:ND1,0:ND2,1:TotNumDomain), &
                       x_vc(0:ND1,0:ND2,1:TotNumDomain), &
                      Ax_vc(0:ND1,0:ND2,1:TotNumDomain), &
                       z_ov(0:ND1,0:ND2,1:TotNumDomain), &
                    z_ov_sol(0:ND1,0:ND2,1:TotNumDomain), &
                     stat = ierr)

      allocate(   xc_in(0:DegNcMax,0:DegNcMax,1:TotNumDomain), &
                     ec(0:DegNcMax,0:DegNcMax,1:TotNumDomain), &
              rc_smooth(0:DegNcMax,0:DegNcMax,1:TotNumDomain), &
                     stat = ierr)

      if (ierr .ne. 0 ) then
         write(*,*)'Cannot allocate memory for Wrapper variables'
         write(*,*)'Abort!'
         stop
      endif
   

      b_wrapper = 0.d0; p_cond =0.d0; w =0.d0 
      x_precond =0.d0; Ax_precond =0.d0
      r_wrapper =0.d0; r_tmp=0.d0
      Proj_AP =0.d0; Proj_P =0.d0
      Pap_pcond = 0.d0; pc_proAp = 0.d0
      x_vc =0.d0; Ax_vc =0.d0
      z_ov = 0.d0; z_ov_sol =0.d0
      xc_in = 0.d0
      ec = 0.d0; ef = 0.d0
      rc_smooth = 0.d0
      
      end subroutine alloc_mem_wrapper_var
      
!----------------------------------------------------------------------      
      subroutine alloc_mem_MLoperator_var(DegMax,TotNumDomain)

      implicit none
      integer:: DegMax(2)
      integer:: TotNumDomain
      integer:: ND1, ND2
      integer:: ND1p, ND2p
      integer:: ierr,step

      ND1 = DegMax(1); ND2 = DegMax(2);
      ND1p = ND1 + 1; ND2p = ND2 + 1;
      allocate( I_U(0:ND1,0:ND2,ND1p*ND2p,1:TotNumDomain), &
                 ML(0:ND1,0:ND2,ND1p*ND2p,1:TotNumDomain), &
               stat=ierr)

      if (ierr .ne. 0 ) then
         write(*,*)'Cannot allocate memory for ML operator variables'
         write(*,*)'Abort!'
         stop
      endif

      I_U=0.d0; ML = 0.d0


      end subroutine 
      subroutine alloc_mem_jacobismooth_var(DegMax,TotNumDomain)

      implicit none
      integer:: DegMax(2)
      integer:: TotNumDomain
      integer:: ND1, ND2
      integer:: ND1p, ND2p
      integer:: ierr,step

      ND1 = DegMax(1); ND2 = DegMax(2);
      ND1p = ND1 + 1; ND2p = ND2 + 1;
      allocate(  r_smooth(0:ND1,0:ND2,1:TotNumDomain), &
                  rc_smooth(0:ND1,0:ND2,1:TotNumDomain), &
                  b_smooth(0:ND1,0:ND2,1:TotNumDomain), &
                     xc_in(0:ND1,0:ND2,1:TotNumDomain), &
                        ec(0:ND1,0:ND2,1:TotNumDomain), &
                        ef(0:ND1,0:ND2,1:TotNumDomain), &
                      x_vc(0:ND1,0:ND2,1:TotNumDomain), &
                     Ax_vc(0:ND1,0:ND2,1:TotNumDomain), &
                        Mr(0:ND1,0:ND2,1:TotNumDomain), &
                         M(1:(ND1p*ND2p*TotNumDomain)), &
               stat=ierr)

      if (ierr .ne. 0 ) then
         write(*,*)'Cannot allocate memory for Jacobi smoothing & 
         operator variables'
         write(*,*)'Abort!'
         stop
      endif

      r_smooth=0.d0; b_smooth=0.d0
      rc_smooth=0.d0
      x_vc =0.d0; Ax_vc =0.d0
      M = 0.d0; Mr =0.d0
      xc_in = 0.d0
      ec = 0.d0; ef = 0.d0


      end subroutine
!----------------------------------------------------------------------      
      subroutine alloc_mem_Interpolatematrix_var(DegMax,TotNumDomain)

      implicit none
      integer:: DegMax(2)
      integer:: TotNumDomain
      integer:: ND1, ND2
      integer:: ND1p, ND2p
      integer:: ierr,step

      ND1 = DegMax(1); ND2 = DegMax(2);
      ND1p = ND1 + 1; ND2p = ND2 + 1;

      allocate( Ihx(0:ND1,0:ND1,1:TotNumDomain), &
                Ihy(0:ND2,0:ND2,1:TotNumDomain), &
                Ihx_transpose(0:ND1,0:ND1,1:TotNumDomain), &
                Ihy_transpose(0:ND2,0:ND2,1:TotNumDomain), &
                mg_jh(0:ND2,0:ND2,1:3), &
                mg_jht(0:ND1,0:ND1,1:3), &
                mg_jhfc(0:ND2,0:ND2,1:3), &
                mg_jhfct(0:ND2,0:ND2,1:3), &
                mg_zh(0:ND1,1:3), &
                xo(0:ND1,1:TotNumDomain), &
                xi(0:ND1,1:TotNumDomain), &
                yo(0:ND2,1:TotNumDomain), &
                yi(0:ND2,1:TotNumDomain), &
                wx(0:ND1,1:2,1:TotNumDomain), &
                wy(0:ND2,1:2,1:TotNumDomain), &
                stat=ierr)

      if (ierr .ne. 0 ) then
         write(*,*)'Cannot allocate memory for Interpolate matrix &
         variables'
         write(*,*)'Abort!'
         stop
      endif

      Ihx = 0.d0; Ihy = 0.d0;
      Ihx_transpose = 0.d0; Ihy_transpose = 0.d0
      xo = 0.d0; xi = 0.d0;
      yo = 0.d0; yi = 0.d0;
      wx = 0.d0; wy = 0.d0

      mg_jh=0.d0; mg_jht=0.d0
      mg_jhfc=0.d0; mg_jhfct=0.d0
      mg_zh=0.d0
      end subroutine
!----------------------------------------------------------------------      
      subroutine alloc_mem_ovlapsmooth_var(DegMax,TotNumDomain)

      implicit none
      integer:: DegMax(2)
      integer:: TotNumDomain
      integer:: ND1, ND2
      integer:: ND1p, ND2p
      integer:: ierr,step

      ND1 = DegMax(1); ND2 = DegMax(2);
      ND1p = ND1 + 1; ND2p = ND2 + 1;
      allocate(  r_smooth(0:ND1,0:ND2,1:TotNumDomain), &
                rc_smooth(0:ND1,0:ND2,1:TotNumDomain), &
                 b_smooth(0:ND1,0:ND2,1:TotNumDomain), &
                    xc_in(0:ND1,0:ND2,1:TotNumDomain), &
                       ec(0:ND1,0:ND2,1:TotNumDomain), &
                       ef(0:ND1,0:ND2,1:TotNumDomain), &
                     x_vc(0:ND1,0:ND2,1:TotNumDomain), &
                    Ax_vc(0:ND1,0:ND2,1:TotNumDomain), &
                     z_ov(0:ND1,0:ND2,1:TotNumDomain), &
                 z_ov_sol(0:ND1,0:ND2,1:TotNumDomain), &
               stat=ierr)

      if (ierr .ne. 0 ) then
         write(*,*)'Cannot allocate memory for overlapping variables'
         write(*,*)'Abort!'
         stop
      endif

      r_smooth=0.d0; rc_smooth=0.d0
      b_smooth=0.d0
      x_vc =0.d0; Ax_vc =0.d0
      z_ov = 0.d0; z_ov_sol =0.d0
      xc_in = 0.d0
      ec = 0.d0; ef = 0.d0


      end subroutine
!----------------------------------------------------------------------      
      subroutine alloc_mem_Lx_Ly_var(DegMax,TotNumDomain,mg_len)

      implicit none
      integer:: DegMax(2)
      integer:: TotNumDomain
      integer:: ND1, ND2
      integer:: ND1p, ND2p
      integer:: ierr,step,mg_len

      ND1 = DegMax(1); ND2 = DegMax(2);
      ND1p = ND1 + 1; ND2p = ND2 + 1;
      allocate(   Lx(0:ND1,0:ND1,1:TotNumDomain), &
                  Ly(0:ND2,0:ND2,1:TotNumDomain), &
                  Lxx(0:ND1,0:ND1,1:TotNumDomain,mg_len), &
                  Lyy(0:ND2,0:ND2,1:TotNumDomain,mg_len), &
                  Bx(0:ND1,0:ND1,1:TotNumDomain), &
                  By(0:ND2,0:ND2,1:TotNumDomain), &
                  Bxx(0:ND1,0:ND1,1:TotNumDomain,mg_len), &
                  Byy(0:ND2,0:ND2,1:TotNumDomain,mg_len), &
                  dudx_tmp(0:ND1,0:ND1,1:TotNumDomain), &
                  dudy_tmp(0:ND2,0:ND2,1:TotNumDomain), &
                  lamx(0:ND1,1:TotNumDomain), &
                  lamy(0:ND2,1:TotNumDomain), &
                  lamxx(0:ND1,1:TotNumDomain,mg_len), &
                  lamyy(0:ND2,1:TotNumDomain,mg_len), &
                  bwx(0:4*(ND1+1)*(ND2+1),1:TotNumDomain), &
                  bwy(0:4*(ND1+1)*(ND2+1),1:TotNumDomain), &
                  bwxx(0:4*(ND1+1)*(ND2+1),1:TotNumDomain,mg_len), &
                  bwyy(0:4*(ND1+1)*(ND2+1),1:TotNumDomain,mg_len), &
                  Diagonal(1:(ND1+1)*(ND2+1),1:TotNumDomain), &
                  DDiagonal(1:(ND1+1)*(ND2+1),1:TotNumDomain,mg_len), &
                  scale_c(1:(ND1+1)*(ND2+1),1:TotNumDomain), &
                  Sx_t(0:ND1,0:ND1,1:TotNumDomain), &
                  Sy_t(0:ND2,0:ND2,1:TotNumDomain), &
                  Sx(0:ND1,0:ND1,1:TotNumDomain), &
                  Sy(0:ND2,0:ND2,1:TotNumDomain), &
                  Sxx_t(0:ND1,0:ND1,1:TotNumDomain,mg_len), &
                  Syy_t(0:ND2,0:ND2,1:TotNumDomain,mg_len), &
                  Sxx(0:ND1,0:ND1,1:TotNumDomain,mg_len), &
                  Syy(0:ND2,0:ND2,1:TotNumDomain,mg_len), &
                  Sx_t_norm(0:ND1,0:ND1,1:TotNumDomain), &
                  Sy_t_norm(0:ND2,0:ND2,1:TotNumDomain), &
                  Sx_norm(0:ND1,0:ND1,1:TotNumDomain), &
                  Sy_norm(0:ND2,0:ND2,1:TotNumDomain), &
                  BSx(0:ND1,0:ND1,1:TotNumDomain), &
                  BSy(0:ND2,0:ND2,1:TotNumDomain), &
                  SxBSx(0:ND1,0:ND1,1:TotNumDomain), &
                  SyBSy(0:ND2,0:ND2,1:TotNumDomain), &
                  stat=ierr)

      if (ierr .ne. 0 ) then
         write(*,*)'Cannot allocate memory for Lx and Ly variables'
         write(*,*)'Abort!'
         stop
      endif

      Lx = 0.d0; Ly = 0.d0
      Lxx = 0.d0; Lyy = 0.d0
      Bx = 0.d0; By = 0.d0
      Bxx = 0.d0; Byy = 0.d0
      dudx_tmp = 0.d0; dudy_tmp = 0.d0
      lamx = 0.d0; lamy = 0.d0
      lamxx = 0.d0; lamyy = 0.d0

      bwx = 0.d0; bwy = 0.d0
      bwxx = 0.d0; bwyy = 0.d0
      Diagonal = 0.d0
      DDiagonal = 0.d0

      Sx_t = 0.d0; Sy_t =0.d0
      Sxx_t = 0.d0; Syy_t =0.d0
      Sx =0.d0; Sy = 0.d0
      Sxx =0.d0; Syy = 0.d0
      Sx_norm = 0.d0; Sx_t_norm=0.d0
      Sy_norm = 0.d0; Sy_t_norm=0.d0

      BSx =0.d0; SxBSx=0.d0
      BSy =0.d0; SyBSy=0.d0

      scale_c = 1.d0
      end subroutine
!----------------------------------------------------------------------      
      subroutine Initial_I_U(DegDM,TotNumDomain)

      implicit none
      integer:: DegDM(2,TotNumDomain)
      integer:: TotNumDomain
      integer:: ND1, ND2
      integer:: ND1p, ND2p
      integer:: DDK, n
      integer:: i, j
      
      do DDK = 1, TotNumDomain
         ND1 = DegDM(1,DDK); ND2 = DegDM(2,DDK)
         ND1p = ND1 + 1; ND2p = ND2 + 1
      
         do j = 0, ND2
            do i = 0, ND1
               n = ( j * ND2p + i  + 1)
               I_U(i,j,n,DDK) = 1
            enddo
         enddo
         
      enddo
      
!      open(1123,file='I_U.text')
!
!      do DDK = 1, TotNumDomain
!         ND1 = DegDM(1,DDK); ND2 = DegDM(2,DDK)
!         ND1p = ND1 + 1; ND2p = ND2 + 1
!
!         do n = 1, ND1p*ND2p
!            write(1123,*)n
!            do j = 0, ND2 
!               do i = 0, ND1 
!                  write(1123,*)i,j,I_U(i,j,n,DDK)
!               enddo
!            enddo
!         enddo
!      enddo
!      close(1123)
      end subroutine Initial_I_U
        
      end module Multigrid_Var
