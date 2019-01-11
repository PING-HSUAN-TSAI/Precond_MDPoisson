!---------------------------------------------------------------!
!                                                               !
! This is a module defining variables for 2D Poisson equation   !
!                                                               !
!---------------------------------------------------------------!
module State_Var  ! Variables for Poisson equation
implicit none

! State Variables
real(kind=8), save, allocatable::  v(:,:,:)     !(0:Deg1_Max,0:Deg2_Max,1:TotNum_DM)
real(kind=8), save, allocatable::  dqdx(:,:)     !temporal storage 
real(kind=8), save, allocatable::  dqdy(:,:)
real(kind=8), save, allocatable::  dvdx(:,:,:)
real(kind=8), save, allocatable::  dvdy(:,:,:)

! Material Variables
real(kind=8), save, allocatable::  a(:) !(1:TotNum_DM)
real(kind=8), save, allocatable::  b(:) !(1:TotNum_DM)
real(kind=8), save, allocatable::  a_pts(:,:,:,:)
real(kind=8), save, allocatable::  b_pts(:,:,:,:)

! Boundary variables  ! p= tau * Lag * (Bu - Gu)
integer, save, allocatable:: BC_Type(:,:)
real(kind=8), save::                  strength      !   ( gamma )
real(kind=8), save::                  c_sigma(1:2)       ! constant of dirichlet-type interface penalty
real(kind=8), save::                  c_tau       ! constant of dirichlet penalty
real(kind=8), save, allocatable::     c_bar(:,:,:,:)  ! scaling factor
real(kind=8), save, allocatable::      tauND(:,:,:,:) ! penalty for non-dirichlet
real(kind=8), save, allocatable::      tauND1(:,:,:,:,:) ! penalty for non-dirichlet
real(kind=8), save, allocatable::      tauND2(:,:,:,:,:)
real(kind=8), save, allocatable::      tauD(:,:,:,:,:)! penalty for dirichlet orthogonal term
real(kind=8), save, allocatable::      tauD2(:,:,:,:,:)! penalty for dirichlet nonorthogonal term
real(kind=8), save, allocatable::      tau_tilde(:,:,:,:,:)   !
real(kind=8), save, allocatable::      tau1(:,:,:,:,:)  ! penalty parameters
real(kind=8), save, allocatable::      tau2(:,:,:,:)  !
real(kind=8), save, allocatable::      tau3(:,:,:,:,:)  !
real(kind=8), save, allocatable::     Sigma(:,:,:,:)  !
real(kind=8), save, allocatable::     SigmaHat(:,:,:,:)
real(kind=8), save, allocatable::     Sigma_tild(:,:,:,:)
real(kind=8), save, allocatable::     SigmaSurr(:,:,:,:)  !
real(kind=8), save, allocatable::     SigmaSurr2(:,:,:,:)  !
real(kind=8), save, allocatable:: bnd_alpha(:,:,:,:)  ! alpha
real(kind=8), save, allocatable::  bnd_beta(:,:,:,:)  ! beta
!real(kind=8), save, allocatable::       PBC(:,:,:,:)  !
!real(kind=8), save, allocatable::      PBC1(:,:,:,:)  !
!real(kind=8), save, allocatable::      PBC2(:,:,:,:)  !
real(kind=8), save, allocatable::        G(:,:,:)
!real(kind=8), save, allocatable::        Bu(:,:,:,:)  ! self domain
!real(kind=8), save, allocatable::        Gu(:,:,:,:)  ! connect domain
!real(kind=8), save, allocatable::   JacEdge(:,:,:,:)  ! Edge_Jacobin
real(kind=8), save, allocatable::     bEdge(:,:,:,:)
!real(kind=8), save, allocatable::     vEdge(:,:,:,:)
!real(kind=8), save, allocatable::  dvdnEdge(:,:,:,:)
!real(kind=8), save, allocatable::    BvEdge(:,:,:,:)
!real(kind=8), save, allocatable::     gEdge(:,:,:,:)
real(kind=8), save, allocatable::     e_first(:,:,:)
real(kind=8), save, allocatable::     e_end(:,:,:)
real(kind=8), save, allocatable::   Nor_dot_gradxi(:,:,:,:)
real(kind=8), save, allocatable::   Nor_dot_gradeta(:,:,:,:)
real(kind=8), save, allocatable::   J_xi_eta(:,:,:,:)
real(kind=8), save, allocatable::   J_eta_xi(:,:,:,:)

! allocate variables which are for testing
real(kind=8), save, allocatable::  rhs(:,:,:)
real(kind=8), save, allocatable::  Differx(:,:,:)
real(kind=8), save, allocatable::  laplace_v(:,:,:)
real(kind=8), save, allocatable::  laplace_vx(:,:,:)
real(kind=8), save, allocatable::  laplace_vy(:,:,:)
real(kind=8), save, allocatable::  dFdx(:,:,:)
real(kind=8), save, allocatable::  dFdy(:,:,:)

real(kind=8), save, allocatable::  Appro_u(:,:,:)
real(kind=8), save, allocatable::  potent(:,:,:)

! 1D penalty 
real(kind=8), save, allocatable::  tau(:,:,:,:)
real(kind=8), save, allocatable::  tautil(:,:,:,:)
real(kind=8), save, allocatable::  SigmaSurr2_1d(:,:,:)
real(kind=8), save, allocatable::  SigmaHat1d(:,:,:)
real(kind=8), save, allocatable::  Sigmabar(:,:,:)
real(kind=8), save, allocatable::  Sigmatil_1(:,:,:,:)
real(kind=8), save, allocatable::  Sigmatil_2(:,:,:)



!real(kind=8), save, allocatable::  error_vec1(:,:,:)
!real(kind=8), save :: error, error_gradx, error_grady, error1
!---------------------------------------------------------------------------
! Declare state variable q
! Where q is the unknown of the following PDE
!
!                v_{tt} = a div dot ( b grad v ) 
!
!--------------------------------------------------------------------------------
!
!    [Name]    :: v(:,:,:)
!     ^^^^
!    [Size]    :: (0:Degree1_Max,0:Degree2_Max,1:TotNum_DM)
!     ^^^^
!    [Purpose] :: Value of the PDE State Variable v 
!     ^^^^^^^
!    [Detail]  ::
!     ^^^^^^
!
!
!    [Name]    :: dvdt(:,:,:)
!     ^^^^
!    [Size]    :: (0:Degree1_Max,0:Degree2_Max,1:TotNum_DM)
!     ^^^^
!    [Purpose] :: To Store the value of the PDE State Variable v 
!     ^^^^^^^
!    [Detail]  :: Value of the time derivative of variable v
!     ^^^^^^
!
!--------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------

subroutine Init_State_Variables(DegMax,TotNumDomain,mgl)
  use MD2D_Grid
  implicit none
  ! declare subroutine arguments
  integer:: DegMax      ! Pass in the PolyDegN_Max
  integer:: TotNumDomain
   
  integer:: NDx,NDy
  integer:: ierr, mgl
   
  ! Allocate memory for the PDE state vector (unknown variable to be solve)
  ! by subroutine alloc_mem_state_variables(Degree_Max,TotNumDomain,Storage)
  ! provided by State_Var Module.
  ! 
  ! alloc_mem_state_variables(Degree_Max,TotNumDomain,Storage) !State_Var
  ! storage = 0 => one level 
  ! storage = 1 => two level
  ! storage = 2 => three level
    
  NDx=DegMax; NDy=DegMax;
  allocate(    v(0:NDx,0:NDy,1:TotNumDomain), &
            dvdx(0:NDx,0:NDy,1:TotNumDomain), &
            dvdy(0:NDx,0:NDy,1:TotNumDomain), &
           a_pts(0:NDx,0:NDy,1:TotNumDomain,mgl), &
           b_pts(0:NDx,0:NDy,1:TotNumDomain,mgl), &
           rhs(0:NDx,0:NDy,1:TotNumDomain), &
           laplace_v(0:NDx,0:NDy,1:TotNumDomain), &
           laplace_vx(0:NDx,0:NDy,1:TotNumDomain), &
           laplace_vy(0:NDx,0:NDy,1:TotNumDomain), &
           dFdx(0:NDx,0:NDy,1:TotNumDomain), &
           dFdy(0:NDx,0:NDy,1:TotNumDomain), &
        Appro_u(0:NDx,0:NDy,1:TotNumDomain), &
         potent(0:NDx,0:NDy,1:TotNumDomain), &
        Differx(0:NDx,0:NDy,1:TotNumDomain), &
                stat=ierr)

  if ( ierr .ne. 0 ) then
     write(*,*)'Cannot allocate memory for state variable '
     write(*,*)'Abort!'
     stop
  endif
  v=0.d0;  Differx=0.0d0
  dvdx=0.d0; dvdy=0.d0;
  rhs=0.d0; dFdx=0.d0; dFdy=0.d0;
  laplace_v=0.d0; laplace_vx=0.d0; laplace_vy=0.d0;

  a_pts=0.d0; b_pts=0.d0 
  Appro_u=0.d0; potent=0.d0

  ND1=DegMax
  allocate(dqdx(0:ND1,0:ND1),dqdy(0:ND1,0:ND1), stat=ierr)
  if (ierr .ne. 0) then 
     write(*,*)' Message from State_Var.f90'
     write(*,*)' Cannot allocate memory for variable dqdx'
     write(*,*)' Abort!'
     stop
  endif 
  dqdx=0.d0; dqdy=0.d0;
  allocate( a(1:TotNumDomain), &
            b(1:TotNumDomain), stat=ierr)

  if (ierr .ne. 0) then
     write(*,*)'Message from State_Var.f90'
     write(*,*)'Can not allocate memory for material variables'
     write(*,*)'Abort!'
     stop
  endif
  a=1.d0; b=1.d0;

  c_sigma = 0.d0
  return 
    
end subroutine Init_State_Variables

!-----------------------------------------------------------------
subroutine alloc_mem_BC_Var(N_max,Num_Domain,mgl)
  use MD2D_Grid
  implicit none
  integer:: N_max, Num_Domain, mgl

  integer:: ierr

  ! subroutine begins
  allocate (  BC_Type(4,Num_Domain), &
                c_bar(0:N_max,4,Num_Domain,mgl), &
                 tauD(0:N_max,0:N_max,4,Num_Domain,mgl),&
                tauD2(0:N_max,0:N_max,4,Num_Domain,mgl),&
               tauND1(0:N_max,0:N_max,4,Num_Domain,mgl),&
               tauND2(0:N_max,0:N_max,4,Num_Domain,mgl),&
                tauND(0:N_max,4,Num_Domain,mgl), &
            tau_tilde(0:N_max,0:N_max,4,Num_Domain,mgl), &
                 tau1(0:N_max,0:N_max,4,Num_Domain,mgl), &
                 tau2(0:N_max,4,Num_Domain,mgl), &
                 tau3(0:N_max,0:N_max,4,Num_Domain,mgl), &
                Sigma(0:N_max,4,Num_Domain,mgl), &
             SigmaHat(0:N_max,4,Num_Domain,mgl), &
           Sigma_tild(0:N_max,4,Num_Domain,mgl), &
            SigmaSurr(0:N_max,4,Num_Domain,mgl), &
           SigmaSurr2(0:N_max,4,Num_Domain,mgl), &
            bnd_alpha(0:N_max,4,Num_Domain,mgl), &
             bnd_beta(0:N_max,4,Num_Domain,mgl), &
!                  PBC(0:N_max,4,Num_Domain,mgl), &
!                 PBC1(0:N_max,4,Num_Domain,mgl), &
!                 PBC2(0:N_max,4,Num_Domain,mgl), &
                    G(0:N_max,4,Num_Domain), &
!                   Bu(0:N_max,4,Num_Domain,mgl), &
!                   Gu(0:N_max,4,Num_Domain,mgl), &
!              JacEdge(0:N_max,4,Num_Domain,mgl), &
                bEdge(0:N_max,4,Num_Domain,mgl), &
!                vEdge(0:N_max,4,Num_Domain,mgl), &
!             dvdnEdge(0:N_max,4,Num_Domain,mgl), &
!               BvEdge(0:N_Max,4,Num_Domain,mgl), &
!                gEdge(0:N_Max,4,Num_Domain,mgl), &
              e_first(0:N_max,Num_Domain,mgl), &
                e_end(0:N_max,Num_Domain,mgl), &
       Nor_dot_gradxi(0:N_max,4,Num_Domain,mgl), &
      Nor_dot_gradeta(0:N_max,4,Num_Domain,mgl), &
             J_eta_xi(0:N_max,4,Num_Domain,mgl), &
             J_xi_eta(0:N_max,4,Num_Domain,mgl), &
             tau(0:N_max,4,Num_Domain,mgl), &
             tautil(0:N_max,4,Num_Domain,mgl), &
             SigmaSurr2_1d(4,Num_Domain,mgl), &
             SigmaHat1d(4,Num_Domain,mgl), &
             Sigmabar(4,Num_Domain,mgl), &
             sigmatil_1(0:N_max,4,Num_Domain,mgl), &
             sigmatil_2(4,Num_Domain,mgl), &
               stat=ierr )

  if ( ierr .ne. 0 ) then
     write(*,*)'Error message from State_Var.f90:'
     write(*,*)'Cannot allocate memory for BC Variables'
     write(*,*)'Abort!'
     stop
  endif
    BC_Type=0

    c_bar=1.00;

    tau1=0.0; tau2=0.0; tau3=0.0; tauND=0.0; tauND1=0.0; tauND2=0.0d0;

    tauD=0.0; tau_tilde=0.0; tauD2=0.0;

    bnd_alpha=0.0; bnd_beta=0.0;

!    PBC=0.0; PBC1=0.0; PBC2=0.0; G=0.0; Bu=0.0; Gu=0.0;

!    JacEdge=0.0; bEdge=0.d0;

    e_first=0.0; e_end=0.0;

Sigma=0.0; SigmaSurr=0.0; SigmaSurr2=0.0; SigmaHat=0.0; Sigma_tild=0.0;

    Nor_dot_gradxi = 0.0; Nor_dot_gradeta=0.0;

    J_xi_eta = 0.0; J_eta_xi=0.0;

    tau = 0.0; tautil=0.0
    SigmaSurr2_1d = 0.0; SigmaHat1d = 0.0; Sigmabar = 0.0
    sigmatil_1 = 0.0; sigmatil_2 = 0.0

  return

end subroutine alloc_mem_BC_Var
!-----------------------------------------------------------------
 
end module State_Var
