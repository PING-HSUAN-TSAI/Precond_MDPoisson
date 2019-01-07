module GMRES
  implicit none
  integer, parameter :: lgmres=50 !max number of gmres iterations between restarts
  real(kind=8),save,allocatable :: gmres_x(:)
  real(kind=8),save,allocatable :: gmres_r(:)
  real(kind=8),save,allocatable :: gmres_w(:)
  real(kind=8),save,allocatable :: gmres_h(:,:)
  real(kind=8),save,allocatable :: gmres_gamma(:)
  real(kind=8),save,allocatable :: gmres_s(:)
  real(kind=8),save,allocatable :: gmres_c(:)
  real(kind=8),save,allocatable :: gmres_v(:,:)
  real(kind=8),save,allocatable :: gmres_z(:,:)  
  real(kind=8),save,allocatable :: gmres_xx(:,:,:)  
  real(kind=8),save,allocatable :: gmres_rr(:,:,:)  
  real(kind=8),save,allocatable :: gmres_ww(:,:,:)  
  real(kind=8),save,allocatable :: gmres_vv(:,:,:,:)  
  real(kind=8),save,allocatable :: gmres_zz(:,:,:,:)  
contains

  subroutine alloc_gmres_var(n,ND1,ND2,TotNum_DM)
    implicit none
    integer :: n, ierr
    integer:: DegMax(2), TotNum_DM
    integer:: ND1, ND2, N_max
    N_max=maxval(DegMax(1:2))
!    ND1=DegMax(1); ND2=DegMax(2)

    allocate ( gmres_x(1:n), gmres_r(1:n), gmres_w(1:n), &
               gmres_h(1:lgmres,1:lgmres), gmres_gamma(1:lgmres+1), &
               gmres_c(1:lgmres), gmres_s(1:lgmres), &
               gmres_v(1:n,1:lgmres), gmres_z(1:n,1:lgmres), &
               gmres_xx(0:ND1,0:ND2,1:TotNum_DM), &
               gmres_rr(0:ND1,0:ND2,1:TotNum_DM), &
               gmres_ww(0:ND1,0:ND2,1:TotNum_DM), &
               gmres_vv(0:ND1,0:ND2,1:TotNum_DM,1:lgmres), &
               gmres_zz(0:ND1,0:ND2,1:TotNum_DM,1:lgmres), &
               stat=ierr)

    gmres_x = 0.d0; gmres_r=0.d0; gmres_w=0.d0
    gmres_h = 0.d0; gmres_gamma=0.d0; gmres_c=0.d0; gmres_s=0.d0
    gmres_v = 0.d0; gmres_z = 0.d0
    gmres_xx = 0.d0; gmres_rr = 0.d0; gmres_ww =0.d0
    gmres_vv = 0.d0; gmres_zz=0.d0
  endsubroutine alloc_gmres_var

endmodule GMRES

    
