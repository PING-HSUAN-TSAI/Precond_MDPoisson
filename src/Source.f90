function v_init(x,y,DDK)
use constants
implicit none
! function arguments
real(kind=8):: v_init
real(kind=8):: x, y, t
real(kind=8):: r, theta, r_scale
real(kind=8):: kNum
integer(kind=8):: n,DDK
! local arguments
real(kind=8), parameter :: WaveNum=1.d0
real(kind=8):: dbesjn
COMPLEX :: XJ
integeR:: ModeNum
real(kind=8):: r_in, r_mid, r_out
real(kind=8):: c_in, c_out
real(kind=8):: AngFre, k_in, k_out
real(kind=8):: A_in, B_in, A_out, B_out
real(kind=8):: amp, arg, u

! XJ=(0,1)
!function begin
!------------------------------------------------------------------------------
! Test case for Paul's paper
   v_init = sin(pi*x) * sin(pi*y)
!------------------------------------------------------------------------------

! Test case for hp convergence of continuous material:
! 1domain,4domain,9domain,16domain

!v_init = cos(2*pi*x) * cos(2*pi*y)
! v_init = cos(pi*x) * cos(pi*y)
!    v_init = cos((0.1*pi)*x) * cos((0.1*pi)*y)

! Test case for 1D problem
!      v_init = 1 + cos(pi*x)
!------------------------------------------------------------------------------

!    v_init = log(exp(y)*cos(x)+4)
!--------------------------------------------------------------------
! Test case for Material Discontinuity at x=0
!select case(DDK)
!case(1,6)
!   !if (x .le. 0.0) then
!      v_init = sin(pi*x) * sin(pi*y)
!case(2,3,4,5)
!   !else if (x .ge. 0.0) then
!      v_init = 5.0/4.0 * sin(2.0*pi*x) *sin(pi*y)
!    !end if!
!end select
!select case(DDK)
!case(1,2,7,8,9,10,15,16)
!!if (x .le. 0.0) then
!v_init = sin(pi*x) * sin(pi*y)
!case(3,4,5,6,11,12,13,14)
!!else if (x .ge. 0.0) then
!v_init = 5.0/4.0 * sin(2.0*pi*x) *sin(pi*y)
!!end if!
!end select
!
!------------------------------------------------------------------------------

! Test case for Circular domain solution

! \phi(r,\theta) = (r**5-4*r**3)*cos(3*\theta)
!   v_init = (x**2.0 + y**2.0 -4) * ((x**3.0)-3*x*(y**2.0))

!------------------------------------------------------------------------------
! Test case for RingMode with material discontinuity at r=2.5
!    n=3
!r=sqrt(x**2.0+y**2.0) ; theta=acos(x/r) ; if (y .lt. 0) theta =-theta
!select case(DDK)
!case(1,3,5,7,9,11,13,15)
!    !if (r .le. 2.5d0) then
!        !r quadratic case no theta change
!            !v_init =  (r**2-4)
!        !r linear case with theta change
!            !v_init = 7*(r-2) * cos(n*theta)
!        !r quadratic case with theta change
!!            v_init = 1.76*(r**2.0-4)*cos(n*theta)
!            v_init = ((r**2.0/4.d0)-1) * cos(n*theta)
!
!!        v_init = (r**2.0-4)*cos(n*theta)
!           !v_init = (-r**2.0+5.4*r-7.2)*cos(n*theta)
!          !v_init = (1.584)*r
!         !v_init = r*cos(n*theta)
!        !r_scale = 2.0 * r / 5.0
!        !v_init = (r_scale**3.0-2*r_scale**2.0) * cos(n*theta)
!!(0.4*(x**2.0)+0.4*(y**2.0)-2) * ((x**3.0)-3*x*(y**2.0))
!!v_init = ((r**3-2*r**2) * cos(n*theta) )
!case(2,4,6,8,10,12,14,16)
!    !else
!            !v_init = 45 * (-r**2.0+5.4*r-7.2)
!            !v_init = (3*r-4) * cos(n*theta)
!            !v_init = (-r**2.0+5.4*r-3.29)*cos(n*theta)
!            v_init = 11.25 * (-r**2.0+5.4*r-7.2) * cos(n*theta)
!!        v_init = 45 * (-r**2.0+5.4*r-7.2) * cos(n*theta)
!
!        !v_init = (-r**2.0+5.4*r-3.29)
!        !v_init = r*cos(n*theta)
!        !r_scale = 2.0 * r / 5.0
!        !v_init = -r_scale**3.0 * cos(n*theta)
!!v_init = ((1.75*r**2-0.5*r**3) * cos(n*theta) )
!end select
    !endif
!------------------------------------------------------------------------------



!    if ( abs(x) .le. 1E-15) then
!        if (y .ge. 0.d0 ) theta = pi/2.d0
!        if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else
!        theta = atan(y/x)
!        if (x .le. 0.d0) theta = theta + pi
!    endif

!cos(pi*x)!3-y**2!3-x**2!2
!v_init = (x**2+y**2-5*sqrt(x**2+y**2)+6)*exp(XJ*2*ATANH(x/y))

!sin(pi*x) * cos(pi*y)
  ! function begin
  
!  v_init = sin(WaveNum*pi*x) * sin(WaveNum*pi*y) & 
!         * cos(sqrt(2.0)*WaveNum*pi*t)
         
!  v_init = sin(pi*(x+y-sqrt(2.0)*t))
  
  ! Stationary Wave Solution
!  v_init = cos(WaveNum*pi*x) * cos(WaveNum*pi*y) &
!         * cos(sqrt(2.d0)*WaveNum*pi*t)
         
!  ! Travelling Wave Solution
!   v_init = sin(pi*(x+y-sqrt(2.0)*t))
   

    !v_init = x**5.0 - 2 * x**3.0 * y**2.0 &
    !       - 3 * x * y**4.0 - 9 * x**3.0 &
    !       + 27 * x * y**2.0

!x**3.0 - 3 * x * y**2.0
!r**3.0 * cos(n*theta) 
!(0.5*r**3.0 - r**2) * cos(n*theta)
!(r**4.0-3*r**3.0+2*r**2.0) * sin(n*theta)


!  ! Pulse Circular 
!    r=sqrt(x**2.0+y**2.0)
!    if (r .le. 1.9d0) then
!       v_init = dexp(-(x**2+y**2)/0.1d0)
!    else
!       v_init =0.d0
!    endif
!    !write(*,*)v_init
    
!!-- Circular Mode solution 
!    kNum=3.190080947961992d0 !1.202412778847886d0 !4.88051564990835d0
!    r = sqrt(x**2+y**2)
!    if ( abs(x) .le. 1E-15) then
!        if (y .ge. 0.d0 ) theta = pi/2.d0
!        if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else 
!        theta = atan(y/x)
!        if (x .le. 0.d0) theta = theta + pi
!    endif
!        
!    v_init = dbesjn(3,kNum*r) * cos(3.d0*theta - kNum*1.d0*t) 
!!!    v_init = dbesjn(3,kNum*r) * dcos(3.d0*theta) *cos( kNum*1.d0*t) 
!!-------------------------------------------------------------------
!-- Circular Ring Model
!    
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.8d0; c_out=2.d0;    
!    ModeNum=4; 
!      amp =1.0         
!    AngFre = 3.994461881329555  !4.179564952653371d0
!    A_in  = -0.467501103782679;  !1.135907404274031E-1; 
!    B_in  =  0.728274364663971;  !8.620436140150400E-1;
!    A_out = -0.129812220564521;  !1.580925996430096E-1;
!    B_out =  0.483950364321976;  !4.679579908102411E-1;
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;  
!    c_in = 0.8d0; c_out=2.d0;                      
!    ModeNum= 5;              
!    amp = 1.0     
!    AngFre = 4.179564952653371d0  
!    A_in  =  1.135907404274031E-1;
!    B_in  =  8.620436140150400E-1;
!    A_out =  1.580925996430096E-1;
!    B_out =  4.679579908102411E-1; 
!
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.5d0; c_out=5.d0;   
!    ModeNum= 8
!     amp =7.0E-4
!    AngFre = 3.458269287950603
!    A_in  =  -0.000933627891442
!    B_in  =  -0.003076893577698
!    A_out =  -0.999994830519074
!    B_out =  -0.000000021132401
!
!
!    k_in = AngFre/c_in; k_out = AngFre/c_out;           
!    r =sqrt(x**2+y**2); 
!    theta = acos(x/r) ; if (y .lt. 0) theta =-theta
!    
!    if (r .le. r_mid) then 
!       arg = k_in*r
!       u = A_in*dbesjn(ModeNum,arg) + B_in*dbesyn(ModeNum,arg)
!       v_init = u * cos(ModeNum*theta - k_in*c_in*t) /amp
!       !v_init = u * cos(ModeNum*theta)*cos(-k_in*c_in*t)
!    else 
!       arg = k_out*r
!       u = A_out*dbesjn(ModeNum,arg) + B_out*dbesyn(ModeNum,arg)
!       v_init = u * cos(ModeNum*theta - k_out*c_out*t) /amp
!       !v_init = u * cos(ModeNum*theta)*cos(-k_out*c_out*t)
!    endif
!---Circular Ring Mode    
  return

end function
!----------------------------------------------------------

!------------------------------------------------------------- 
function dFdx_init(x,y,DDK)
use constants
implicit none
! function arguments
real(kind=8):: dFdx_init
real(kind=8):: x, y, t
real(kind=8):: r, theta, r_scale,const1,const2
real(kind=8):: kNum
integer(kind=8):: n,DDK
! local arguments
real(kind=8), parameter :: WaveNum=1.d0
complex :: XJ
  
integer:: ModeNum
real(kind=8):: r_in, r_mid, r_out
real(kind=8):: c_in, c_out
real(kind=8):: AngFre, k_in, k_out
real(kind=8):: A_in, B_in, A_out, B_out
real(kind=8):: u, arg , amp,tmp1,tmp2

!XJ=(0,1)

! function begin
!------------------------------------------------------------------------------
! Test case for Paul's paper
   dFdx_init = pi * cos(pi*x) * sin(pi*y)
!------------------------------------------------------------------------------
! Test case for hp convergence of continuous material:
! 1domain,4domain,9domain,16domain

!dFdx_init = - 2 * pi * sin(2*pi*x) * cos(2*pi*y)
!  dFdx_init = - pi * sin(pi*x) * cos(pi*y)
! dFdx_init = - (0.1*pi) * sin((0.1*pi)*x) * cos((0.1*pi)*y)

! Test case for 1D problem
!      dFdx_init = - pi * sin(pi*x)
!-----------------------------------------------------------------------------
!  dFdx_init = ( -exp(y)*sin(x)) / (exp(y)*cos(x)+4)
!-------------------------------------------------------------------
! Test case for Material Discontinuity at x=0
!select case(DDK)
!case(1,6)
!   !if (x .le. 0.0) then
!      dFdx_init = pi * cos(pi*x) * sin(pi*y)
!case(2,3,4,5)
!   !else if (x .ge. 0.0) then
!      dFdx_init = 5.0/4.0 * 2 * pi * cos(2.0*pi*x) *sin(pi*y)
!   !endif
!end select
!select case(DDK)
!case(1,2,7,8,9,10,15,16)
!!if (x .le. 0.0) then
!dFdx_init = pi * cos(pi*x) * sin(pi*y)
!case(3,4,5,6,11,12,13,14)
!!else if (x .ge. 0.0) then
!dFdx_init = 5.0/4.0 * 2 * pi * cos(2.0*pi*x) *sin(pi*y)
!end select

!------------------------------------------------------------------------------

! Test case for Circular Domain
!    dFdx_init = 5 * (x**4.0) - 6 * (x**2.0) * (y**2.0) &
!              - 3 * (y**4.0) &
!              - 12 * (x**2.0) + 12 * (y**2.0)
!------------------------------------------------------------------------------
! Test case for RingMode with material discontinuity at r=2.5
!    n=3
!r=sqrt(x**2.0+y**2.0) ; theta=acos(x/r) ;!atan2(y,x)
!    if (y .lt. 0) theta =-theta
!select case(DDK)
!case(1,3,5,7,9,11,13,15)
!!    if (r .le. 2.5d0) then
!           !dFdx_init = 2*r*cos(theta)
!           ! dFdx_init = 7 * cos(n*theta) *cos(theta) &
!           !           + n * (7.0-14.0/r) * sin(n*theta) * sin(theta)
!    !dFdx_init = 1.76 * (2*r) * cos(n*theta)*cos(theta) &
!    !          + n * 1.76 * (r-(4/r)) * sin(n*theta) * sin(theta)
!    dFdx_init = (r/2.d0) * cos(n*theta)*cos(theta) &
!        + n * ((r/4.d0)-(1/r)) * sin(n*theta) * sin(theta)
!
!!    dFdx_init = (2*r) * cos(n*theta)*cos(theta) &
!!              + n * (r-(4/r)) * sin(n*theta) * sin(theta)
!
!           ! dFdx_init = (1.584) * cos(theta)
!          !dFdx_init = (-2*r+5.4)*cos(n*theta)*cos(theta) &
!           !       + n * (-r + 5.4 -(7.2/r) ) *sin(n*theta)*sin(theta)
!!          dFdx_init = cos(n*theta) * cos(theta) &
!!                    + n * sin(n*theta) * sin(theta)
!!         dFdx_init = &
!!  (3*const1*r**2.0-4*const2*r)*cos(n*theta)*cos(theta) &
!!+ n*(const1*r**2.0-2*const2*r)*sin(n*theta)*sin(theta)
!!!(2*(x**4.0)-2.4*(x**2.0)*(y**2.0)-1.2*(y**4.0)-6*(x**2.0)+6*(y**2.0))
!!        dFdx_init = ( (3*r**2-4*r) * cos(n*theta) * cos(theta) ) &
!!                  + ( 3*(r**2-2*r) * sin(n*theta) * sin(theta) )
!!    else
!case(2,4,6,8,10,12,14,16)
!            !dFdx_init = 45 * (-2*r+5.4) * cos(theta)
!           ! dFdx_init = 3* cos(n*theta) *cos(theta) &
!           !           + n * (3.0 - 4.0/r) * sin(n*theta) * sin(theta)
!    !dFdx_init = (-2*r+5.4)*cos(n*theta)*cos(theta) &
!    !          + n * (-r+5.4-(3.29/r))*sin(n*theta)*sin(theta)
!
!    dFdx_init = 11.25 * (-2*r+5.4)*cos(n*theta)*cos(theta) &
!              + n * 11.25 * (-r+5.4-(7.2/r))*sin(n*theta)*sin(theta)
!
!
!          ! dFdx_init = (-2*r+5.4) * cos(theta)
!         !dFdx_init = (-2*r+5.4)*cos(n*theta)*cos(theta) &
!         !           + n * (-r + 5.4 -(7.2/r) ) *sin(n*theta)*sin(theta)
!!        dFdx_init = cos(n*theta) * cos(theta) &
!!                  + n * sin(n*theta) * sin(theta)
!!        dFdx_init = (-3*const1*r**2.0)*cos(n*theta)*cos(theta) &
!!            - n*(const1*r**2.0)*sin(n*theta)*sin(theta)
!!(0.5) * (3*(x**2.0)-3*(y**2.0))
!!        dFdx_init = ( (3.5*r-1.5*r**2) * cos(n*theta) * cos(theta) ) &
!!                  + ( 3*(1.75*r-0.5*r**2) * sin(n*theta) * sin(theta) )
!end select
!    endif
!------------------------------------------------------------------------------



!-pi * sin(pi*x)!0!-2*x!0!pi * COS(pi*x) * COS(pi*y)

  !dFdx_init = (2*x - (5 * x / sqrt(x**2+y**2)))* exp(XJ*2*ATANH(x/y)) &
          !  + ((XJ*2*y)/(x**2+y**2)) &
          !  * ( x**2 + y**2 - 5 * sqrt(x**2+y**2) + 6 ) &
          !  * exp(XJ*2*ATANH(x/y))

 ! dvdt_init = -sqrt(2.0) * WaveNum * pi & 
 !           *  sin(WaveNum*pi*x) * sin(WaveNum*pi*y) & 
 !           *  sin(sqrt(2.0)*WaveNum*pi*t)
            
!  dvdt_init = -sqrt(2.0)*pi * cos(pi*(x+y-sqrt(2.0)*t))

  ! Stationary Wave Solution 
!  dvdt_init = -sqrt(2.0) * WaveNum * pi &
!            *  cos(WaveNum*pi*x) * cos(WaveNum*pi*y) &
!            *  sin(sqrt(2.d0)*WaveNum*pi*t)  
                    
  ! Travelling Wave Solution
!  dvdt_init = -sqrt(2.0)*pi * cos(pi*(x+y-sqrt(2.0)*t))
        


!  ! pulse circular
!  dvdt_init = 0.d0
!  
!   Circular domain case
!    n=3
!r=sqrt(x**2.0+y**2.0) ; theta=atan2(y,x)
!
!    if (r .le. 2.5d0) then
!        dFdx_init = (3 * r**2 * cos(n*theta) * cos(theta) &
!                  + 3 * (r**3-8) * sin(n*theta) * sin(theta) /r) &
!                  /(7.625)
!    else
!        dFdx_init = (3 * r**2 * cos(n*theta) * cos(theta) &
!                  + 3 * (r**3-27) * sin(n*theta) * sin(theta) /r) &
!                  / (-11.375)
!    endif


!dFdx_init = 5 * x**4.0 - 6 * x**2.0 * y**2.0 &
!- 3 * y**4.0 - 27 * x**2.0 + 27 * y**2.0
!    if ( abs(x) .le. 1E-15) then
!    if (y .ge. 0.d0 ) theta = pi/2.d0
!    if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else
!    theta = atan(y/x)
!    if (x .le. 0.d0) theta = theta + pi
!    endif
!tmp1 = 24 * r * x**4.0 &
     !- 18 * r**3.0 * x**2.0
!tmp2 = 9 * r * y**4.0 &
    ! - 12 * r**3.0 * y**2.0
!    tmp1 = (1.5*r**2.0-2*r) &
!    * cos(theta) * cos(n*theta)
!    tmp2 = n*(0.5*r**2.0-r) &
!    * sin(n*theta) * sin(theta)
!    tmp1 = (4*r**3.0-9*r**2.0+4*r) &
!         * cos(theta) * sin(n*theta)
!    tmp2 = - n*(r**3.0-3*r**2.0+2*r) &
!         * cos(n*theta) * sin(theta)



!3 * x**2.0 - 3 * y**2.0

  ! Circular mode 
!  kNum=3.190080947961992d0 !1.202412778847886d0 !4.88051564990835d0
!    r = dsqrt(x**2.d0+y**2.d0)
!    if ( abs(x) .le. 1E-15) then
!        if (y .ge. 0.d0 ) theta = pi/2.d0
!        if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else 
!        theta = datan(dble(y/x))
!        if (x .le. 0.d0) theta = theta + pi
!    endif
!        
!    dvdt_init = kNum*1.d0*bessel_jn(3,kNum*r) * dsin(dble(3.d0*theta - kNum*1.d0*t)) 
! !   dvdt_init=0.d0 
 
!!-------------------------------------------------------------------
!-- Circular Ring Model
!    
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.8d0; c_out=2.d0;    
!    ModeNum=4;  
!     amp=1.0;               
!    AngFre = 3.994461881329555;
!    A_in  = -0.467501103782679;
!    B_in  =  0.728274364663971;
!    A_out = -0.129812220564521;
!    B_out =  0.483950364321976;    
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.8d0; c_out=2.d0;
!    ModeNum= 5;
!     amp=1.0;                                   
!    AngFre = 4.179564952653371d0   
!    A_in  =  1.135907404274031E-1; 
!    B_in  =  8.620436140150400E-1; 
!    A_out =  1.580925996430096E-1; 
!    B_out =  4.679579908102411E-1;
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.5d0; c_out=5.d0;   
!    ModeNum= 8
!     amp=7.0E-4
!    AngFre = 3.458269287950603  !3.458269287950603
!    A_in  =  -0.000933627891442
!    B_in  =  -0.003076893577698
!    A_out =  -0.999994830519074
!    B_out =  -0.000000021132401
!    
!    k_in = AngFre/c_in; k_out = AngFre/c_out;       
!    r =sqrt(x**2+y**2);                             
!    theta = acos(x/r) ; if (y .lt. 0) theta =-theta
!    
!    if (r .le. r_mid) then 
!       arg = k_in*r
!       u = A_in*dbesjn(ModeNum,arg) + B_in*dbesyn(ModeNum,arg)
!       dvdt_init = k_in *c_in *  u * sin (ModeNum*theta - k_in*c_in*t)/amp
!    else 
!       arg = k_out*r
!       u = A_out*dbesjn(ModeNum,arg) + B_out*dbesyn(ModeNum,arg)
!       dvdt_init = k_out*c_out * u * sin(ModeNum*theta - k_out*c_out*t)/amp
!    endif
!       !dvdt_init=0.d0
!---Circular Ring Mode     
       
  return

end function 

!-------------------------------------------------------------
function dFdy_init(x,y,DDK)
use constants
implicit none
! function arguments
real(kind=8):: dFdy_init
real(kind=8):: x, y, t
real(kind=8):: r, theta, r_scale,const1,const2
real(kind=8):: kNum
integer(kind=8):: n,DDK
! local arguments
real(kind=8), parameter :: WaveNum=1.d0
complex :: XJ
integer:: ModeNum
real(kind=8):: r_in, r_mid, r_out
real(kind=8):: c_in, c_out
real(kind=8):: AngFre, k_in, k_out
real(kind=8):: A_in, B_in, A_out, B_out
real(kind=8):: u, arg , amp,tmp1,tmp2

!XJ=(0,1)
! function begin
!------------------------------------------------------------------------------
! Test case for Paul's paper
   dFdy_init = (pi) * sin(pi*x) * cos(pi*y)

!------------------------------------------------------------------------------
! Test case for hp convergence of continuous material:
! 1domain,4domain,9domain,16domain

!  dFdy_init = - 2 * pi * cos(2*pi*x) * sin(2*pi*y)
!   dFdy_init = - pi * cos(pi*x) * sin(pi*y)
!dFdy_init = - (0.1*pi) * cos((0.1*pi)*x) * sin((0.1*pi)*y)

! Test case for 1D problem
!      dFdy_init = 0.0
!----------------------------------------------------------------------------
 !  dFdy_init = (exp(y)*cos(x)) / (exp(y)*cos(x)+4)
!-------------------------------------------------------------------
! Test case for Material Discontinuity at x=0
!select case(DDK)
!case(1,6)
!   ! if (x .le. 0.0) then
!        dFdy_init =  pi * sin(pi*x) * cos(pi*y)
!case(2,3,4,5)
!   ! else if (x .ge. 0.0) then
!        dFdy_init =  5.0/4.0 *  pi * sin(2.0*pi*x) *cos(pi*y)
!end select
!select case(DDK)
!case(1,2,7,8,9,10,15,16)
!! if (x .le. 0.0) then
!dFdy_init =  pi * sin(pi*x) * cos(pi*y)
!case(3,4,5,6,11,12,13,14)
!! else if (x .ge. 0.0) then
!dFdy_init =  5.0/4.0 *  pi * sin(2.0*pi*x) *cos(pi*y)
!end select

   ! endif
!------------------------------------------------------------------------------
! Test case for Circular Domain
!    dFdy_init = -4 * (x**3.0) * y - 12 * x * (y**3.0) &
!              + 24 * x *y
!------------------------------------------------------------------------------
! Test case for RingMode with material discontinuity at r=2.5
!    n=3
!    r=sqrt(x**2.0+y**2.0) ; theta=acos(x/r)
!    if (y .lt. 0) theta =-theta!atan2(y,x)
!select case(DDK)
!case(1,3,5,7,9,11,13,15)
!!    const1 = (2.0/5.0)**3.0
!!    const2 = (2.0/5.0)**2.0
!!    if (r .le. 2.5d0) then
!             !dFdy_init = 2 * r* sin(theta)
!            !dFdy_init = 7 * cos(n*theta) *sin(theta) &
!            !          - n * (7.0-14.0/r) *sin(n*theta) * cos(theta)
!!    dFdy_init = 1.76 * (2*r) * cos(n*theta) * sin(theta) &
!!              - n * 1.76 * (r-(4/r)) * sin(n*theta) * cos(theta)
!    dFdy_init = (r/2.d0) * cos(n*theta) * sin(theta) &
!                - n * ((r/4.d0)-(1/r)) * sin(n*theta) * cos(theta)
!
!!     dFdy_init = (2*r) * cos(n*theta) * sin(theta) &
!!               - n * (r-(4/r)) * sin(n*theta) * cos(theta)
!
!          !  dFdy_init = (1.584) * sin(theta)
!          !dFdy_init = (-2*r+5.4)*cos(n*theta)*sin(theta) &
!          !           - n * (-r+5.4-(7.2/r))*sin(n*theta)*cos(theta)
!!          dFdy_init = cos(n*theta) * sin(theta) &
!!                    - n * sin(n*theta) * cos(theta)
!!         dFdy_init = &
!!  (3*const1*r**2.0-4*const2*r)*cos(n*theta)*sin(theta) &
!!- n*(const1*r**2.0-2*const2*r)*sin(n*theta)*cos(theta)
!!(-1.6*(x**3.0)*y - 4.8*x*(y**3.0) + 12*x*y)
!!        dFdy_init = ( (3*r**2-4*r) * cos(n*theta) * sin(theta) ) &
!!                  - ( 3*(r**2-2*r) * sin(n*theta) * cos(theta) )
!!    else
!case(2,4,6,8,10,12,14,16)
!            !dFdy_init = 45 * (-2*r+5.4) * sin(theta)
!            !dFdy_init = 3 * cos(n*theta) * sin(theta) &
!            !          - n * (3.0-4.0/r) * sin(n*theta) * cos(theta)
!!    dFdy_init = (-2*r+5.4)*cos(n*theta)*sin(theta) &
!!              - n * (-r+5.4-(3.29/r)) * sin(n*theta) * cos(theta)
!dFdy_init = 11.25 * (-2*r+5.4)*cos(n*theta)*sin(theta) &
!- n * 11.25 * (-r+5.4-(7.2/r)) * sin(n*theta) * cos(theta)
!
!!     dFdy_init = 45 * (-2*r+5.4)*cos(n*theta)*sin(theta) &
!!               - n * 45 * (-r+5.4-(7.2/r)) * sin(n*theta) * cos(theta)
!
!
!          !dFdy_init = (-2*r+5.4) * sin(theta)
!         ! dFdy_init = (-2*r+5.4)*cos(n*theta)*sin(theta) &
!         !           - n * (-r+5.4-(7.2/r))*sin(n*theta)*cos(theta)
!!        dFdy_init = cos(n*theta) * sin(theta) &
!!                  - n * sin(n*theta) * cos(theta)
!!        dFdy_init = (-3*const1*r**2.0)*cos(n*theta)*sin(theta) &
!!                  + n*(const1*r**2.0)*sin(n*theta)*cos(theta)
!!(0.5) * (-6*x*y)
!
!!        dFdy_init = ( (3.5*r-.15*r**2) * cos(n*theta) * sin(theta) ) &
!!                  - ( 3*(1.75*r-0.5*r**2) * sin(n*theta) * cos(theta) )
!end select
!    endif
!------------------------------------------------------------------------------

!dFdy_init = (2*y - (5 * y / sqrt(x**2+y**2)))* exp(XJ*2*ATANH(x/y)) &
           ! - ((XJ*2*x)/(x**2+y**2)) &
           ! * ( x**2 + y**2 - 5 * sqrt(x**2+y**2) + 6 ) &
           ! * exp(XJ*2*ATANH(x/y))

! dvdt_init = -sqrt(2.0) * WaveNum * pi &
!           *  sin(WaveNum*pi*x) * sin(WaveNum*pi*y) &
!           *  sin(sqrt(2.0)*WaveNum*pi*t)

!  dvdt_init = -sqrt(2.0)*pi * cos(pi*(x+y-sqrt(2.0)*t))

! Stationary Wave Solution
!  dvdt_init = -sqrt(2.0) * WaveNum * pi &
!            *  cos(WaveNum*pi*x) * cos(WaveNum*pi*y) &
!            *  sin(sqrt(2.d0)*WaveNum*pi*t)

! Travelling Wave Solution
!  dvdt_init = -sqrt(2.0)*pi * cos(pi*(x+y-sqrt(2.0)*t))



!  ! pulse circular
!  dvdt_init = 0.d0
!
!   Circular domain case
!    n=3
!    r=sqrt(x**2.0+y**2.0) ; theta=atan2(y,x)
!
!    if (r .le. 2.5d0) then
!        dFdy_init = (3 * r**2 * cos(n*theta) * sin(theta) &
!                  - 3 * (r**3-8) * sin(n*theta) * cos(theta) /r) &
!                  / (7.625)
!    else
!        dFdy_init = (3 * r**2 * cos(n*theta) * sin(theta) &
!                  - 3 * (r**3-27) * sin(n*theta) * cos(theta) /r) &
!                  / (-11.375)
!    endif

!dFdy_init = -4 * x**3.0 * y - 12 * x * y**3.0 &
!+ 54 * x * y

!    if ( abs(x) .le. 1E-15) then
!    if (y .ge. 0.d0 ) theta = pi/2.d0
!    if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else
!    theta = atan(y/x)
!    if (x .le. 0.d0) theta = theta + pi
!    endif
!tmp1 = 24 * r * x**3.0 * y &
!     - 18 * r**3.0 * x * y
!tmp2 = -9 * r * y**3.0 * x &
!     + 12 * r**3.0 * y * x
!    tmp1 = (1.5*r**2.0-2*r) &
!    * cos(n*theta) * sin(theta)
!    tmp2 = - n * (0.5*r**2.0-r) &
!    * sin(n*theta) * cos(theta)

!    tmp1 = (4*r**3.0-9*r**2.0+4*r) &
!         * sin(n*theta) * sin(theta)
!    tmp2 = n * (r**3.0-3*r**2.0+2*r) &
!         * cos(n*theta) * cos(theta)

!-6 * x * y

! Circular mode
!  kNum=3.190080947961992d0 !1.202412778847886d0 !4.88051564990835d0
!    r = dsqrt(x**2.d0+y**2.d0)
!    if ( abs(x) .le. 1E-15) then
!        if (y .ge. 0.d0 ) theta = pi/2.d0
!        if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else
!        theta = datan(dble(y/x))
!        if (x .le. 0.d0) theta = theta + pi
!    endif
!
!    dvdt_init = kNum*1.d0*bessel_jn(3,kNum*r) * dsin(dble(3.d0*theta - kNum*1.d0*t))
! !   dvdt_init=0.d0

!!-------------------------------------------------------------------
!-- Circular Ring Model
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.8d0; c_out=2.d0;
!    ModeNum=4;
!     amp=1.0;
!    AngFre = 3.994461881329555;
!    A_in  = -0.467501103782679;
!    B_in  =  0.728274364663971;
!    A_out = -0.129812220564521;
!    B_out =  0.483950364321976;
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.8d0; c_out=2.d0;
!    ModeNum= 5;
!     amp=1.0;
!    AngFre = 4.179564952653371d0
!    A_in  =  1.135907404274031E-1;
!    B_in  =  8.620436140150400E-1;
!    A_out =  1.580925996430096E-1;
!    B_out =  4.679579908102411E-1;
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.5d0; c_out=5.d0;
!    ModeNum= 8
!     amp=7.0E-4
!    AngFre = 3.458269287950603  !3.458269287950603
!    A_in  =  -0.000933627891442
!    B_in  =  -0.003076893577698
!    A_out =  -0.999994830519074
!    B_out =  -0.000000021132401
!
!    k_in = AngFre/c_in; k_out = AngFre/c_out;
!    r =sqrt(x**2+y**2);
!    theta = acos(x/r) ; if (y .lt. 0) theta =-theta
!
!    if (r .le. r_mid) then
!       arg = k_in*r
!       u = A_in*dbesjn(ModeNum,arg) + B_in*dbesyn(ModeNum,arg)
!       dvdt_init = k_in *c_in *  u * sin (ModeNum*theta - k_in*c_in*t)/amp
!    else
!       arg = k_out*r
!       u = A_out*dbesjn(ModeNum,arg) + B_out*dbesyn(ModeNum,arg)
!       dvdt_init = k_out*c_out * u * sin(ModeNum*theta - k_out*c_out*t)/amp
!    endif
!       !dvdt_init=0.d0
!---Circular Ring Mode

return

end function

!----------------------------------------------------------
function F_init(x,y,DDK)
use constants
implicit none
! function arguments
real(kind=8):: F_init
real(kind=8):: x, y, t
real(kind=8):: r, theta,r_scale,const1,const2
real(kind=8):: kNum
integer(kind=8):: n,DDK
! local arguments
real(kind=8), parameter :: WaveNum=1.d0
complex :: XJ
integer:: ModeNum
real(kind=8):: r_in, r_mid, r_out
real(kind=8):: c_in, c_out
real(kind=8):: AngFre, k_in, k_out
real(kind=8):: A_in, B_in, A_out, B_out
real(kind=8):: u, arg , amp,tmp1,tmp2

!XJ=(0,1)
! function begin
!------------------------------------------------------------------------------
!     Test case for Paul's paper
      F_init = 2 * pi**2 * sin(pi*x) * sin(pi*y) 
!------------------------------------------------------------------------------
! Test case for hp convergence of continuous material:
! 1domain,4domain,9domain,16domain

!   F_init = 8 * pi**2 * cos(2*pi*x)*cos(2*pi*y)
!    F_init = 2 * pi**2 * cos(pi*x) * cos(pi*y)
!F_init = (0.02*pi**2) * cos((0.1*pi)*x)*cos((0.1*pi)*y)

! Test case for 1D problem
!      F_init = pi**2 * cos(pi*x)
!------------------------------------------------------------------------------
!    F_init = (exp(2*y)) &
!           / (exp(y)*cos(x)+4)**2.0
!-------------------------------------------------------------------
! Test case for Material Discontinuity at x=0
!select case(DDK)
!case(1,6)
!    !if (x .le. 0.0) then
!        F_init = 1 * pi**2 * sin(pi*x) * sin(pi*y)
!case(2,3,4,5)
!    !else if (x .ge. 0.0) then
!        F_init = 5.0/4.0 * 1  * pi**2 * sin(2.0*pi*x) *sin(pi*y)
!end select
!select case(DDK)
!case(1,2,7,8,9,10,15,16)
!!if (x .le. 0.0) then
!F_init = 1 * pi**2 * sin(pi*x) * sin(pi*y)
!case(3,4,5,6,11,12,13,14)
!!else if (x .ge. 0.0) then
!F_init = 5.0/4.0 * 1  * pi**2 * sin(2.0*pi*x) *sin(pi*y)
!end select
!endif
!------------------------------------------------------------------------------
! Test case for Circular Domain
!        
!        F_init = -16 * (x**3.0) + 48 * x *(y**2.0)

!------------------------------------------------------------------------------
! Test case for RingMode with material discontinuity at r=2.5
!    n=3
!    r=sqrt(x**2.0+y**2.0) ; theta=acos(x/r)
!    if (y .lt. 0) theta =-theta!atan2(y,x)
!select case(DDK)
!case(1,3,5,7,9,11,13,15)
!!    const1 = (2.0/5.0)**3.0
!!    const2 = (2.0/5.0)**2.0
!!    if (r .le. 2.5d0) then
!            !F_init = -3 * (-56.0/r + 126.0/r**2) * cos(n*theta)
!            !F_init = -4 * 18
!
!    !F_init = -1 * (-8.8+ 1.76*(36/r**2.0))*cos(n*theta)
!    F_init =  -18 * (-1.25 + (9/r**2.0))*cos(n*theta)
!
!!     F_init =  -18 * (-5 + (36/r**2.0))*cos(n*theta)
!            !F_init = -3 * 7.0/r
!           !F_init = -(1.584)/r
!          !F_init = -(5-(43.2/r)+(64.8/r**2.0))*cos(n*theta)
!         !F_init = (8/r)*cos(n*theta)
! !        F_init = (-3) * (5*r*const1*cos(n*theta))
!!(-3) * (6.4*(x**3.0)-19.2*x*(y**2.0))
!!        F_init = (1.0)*(-10)*cos(n*theta)
!!    else
!
!case(2,4,6,8,10,12,14,16)
!          !F_init = -7 * (-24.0/r + 36.0/r**2) *cos(n*theta)
!           !F_init = -45 * 5 * (-4+(5.4/r))
!!    F_init = -22 * (5 - (43.2/r) + (29.61/r**2)) * cos(n*theta)
!    F_init = -5 * 11.25 * (5 - (43.2/r) + (64.8/r**2)) * cos(n*theta)
!
!!     F_init = -5 * (225 - (1944/r) + (2916/r**2)) * cos(n*theta)
!
!          !F_init = -7 * 3.0/r
!         !F_init = -3.96 * (-4 + (5.4/r))
!        !F_init = -(5-(43.2/r)+(64.8/r**2.0))*cos(n*theta)
!        !F_init = (8/r)*cos(n*theta)
! !        F_init = (-1) * (-5*r*const1*cos(n*theta))
!!0!(-44) * (16*(x**3.0)-48*x*(y**2.0))
!!        F_init = (-14.0)*(8.75)*cos(n*theta)
!end select
 !   endif
!------------------------------------------------------------------------------


!F_init = (4- (5/sqrt(x**2+y**2))) * exp(XJ*2*ATANH(x/y)) &
!       - 2**2 * (1-(5/sqrt(x**2+y**2)) + (6/x**2+y**2))  &
!       * exp(XJ*2*ATANH(x/y))

! dvdt_init = -sqrt(2.0) * WaveNum * pi &
!           *  sin(WaveNum*pi*x) * sin(WaveNum*pi*y) &
!           *  sin(sqrt(2.0)*WaveNum*pi*t)

!  dvdt_init = -sqrt(2.0)*pi * cos(pi*(x+y-sqrt(2.0)*t))

! Stationary Wave Solution
!  dvdt_init = -sqrt(2.0) * WaveNum * pi &
!            *  cos(WaveNum*pi*x) * cos(WaveNum*pi*y) &
!            *  sin(sqrt(2.d0)*WaveNum*pi*t)

! Travelling Wave Solution
!  dvdt_init = -sqrt(2.0)*pi * cos(pi*(x+y-sqrt(2.0)*t))

!  ! pulse circular
!  dvdt_init = 0.d0
!
!   Circular domain case
!    n=3
!    r=sqrt(x**2.0+y**2.0) ; theta=atan2(y,x)
!    if (r .le. 2.5d0) then
!        F_init = -(9*8/r**2)*cos(n*theta)
!        !F_init = - ((9*r**2*cos(n*theta))/r &
!        !       - (9*r**3-72)*cos(n*theta)/r**2)
!    else
!        F_init = -(9*27/r**2)*cos(n*theta)
!        !F_init = -((9*r**2*cos(n*theta))/r &
!        !       - (9*r**3-243)*cos(n*theta)/r**2)
!    endif



!    F_init = -20 * x**3.0 + 48 * x * y**2.0 + 4 * x**3.0

!    if ( abs(x) .le. 1E-15) then
!    if (y .ge. 0.d0 ) theta = pi/2.d0
!    if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else
!    theta = atan(y/x)
!    if (x .le. 0.d0) theta = theta + pi
!    endif

!-(7*r**2-10) * sin(n*theta)
!-(16*r**2.0-27*r+8)


! Circular mode
!  kNum=3.190080947961992d0 !1.202412778847886d0 !4.88051564990835d0
!    r = dsqrt(x**2.d0+y**2.d0)
!    if ( abs(x) .le. 1E-15) then
!        if (y .ge. 0.d0 ) theta = pi/2.d0
!        if (y .lt. 0.d0 ) theta =-pi/2.d0
!    else
!        theta = datan(dble(y/x))
!        if (x .le. 0.d0) theta = theta + pi
!    endif
!
!    dvdt_init = kNum*1.d0*bessel_jn(3,kNum*r) * dsin(dble(3.d0*theta - kNum*1.d0*t))
! !   dvdt_init=0.d0

!!-------------------------------------------------------------------
!-- Circular Ring Model
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.8d0; c_out=2.d0;
!    ModeNum=4;
!     amp=1.0;
!    AngFre = 3.994461881329555;
!    A_in  = -0.467501103782679;
!    B_in  =  0.728274364663971;
!    A_out = -0.129812220564521;
!    B_out =  0.483950364321976;
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.8d0; c_out=2.d0;
!    ModeNum= 5;
!     amp=1.0;
!    AngFre = 4.179564952653371d0
!    A_in  =  1.135907404274031E-1;
!    B_in  =  8.620436140150400E-1;
!    A_out =  1.580925996430096E-1;
!    B_out =  4.679579908102411E-1;
!
!    r_in = 2.d0; r_mid = 2.5d0; r_out = 3.d0;
!    c_in = 0.5d0; c_out=5.d0;
!    ModeNum= 8
!     amp=7.0E-4
!    AngFre = 3.458269287950603  !3.458269287950603
!    A_in  =  -0.000933627891442
!    B_in  =  -0.003076893577698
!    A_out =  -0.999994830519074
!    B_out =  -0.000000021132401
!
!    k_in = AngFre/c_in; k_out = AngFre/c_out;
!    r =sqrt(x**2+y**2);
!    theta = acos(x/r) ; if (y .lt. 0) theta =-theta
!
!    if (r .le. r_mid) then
!       arg = k_in*r
!       u = A_in*dbesjn(ModeNum,arg) + B_in*dbesyn(ModeNum,arg)
!       dvdt_init = k_in *c_in *  u * sin (ModeNum*theta - k_in*c_in*t)/amp
!    else
!       arg = k_out*r
!       u = A_out*dbesjn(ModeNum,arg) + B_out*dbesyn(ModeNum,arg)
!       dvdt_init = k_out*c_out * u * sin(ModeNum*theta - k_out*c_out*t)/amp
!    endif
!       !dvdt_init=0.d0
!---Circular Ring Mode

return

end function

!----------------------------------------------------------
function g_ext(x,y,nx,ny,alpha,beta,b,Edge_Num)
use constants
implicit none
integer:: Edge_Num,n
real(kind=8):: g_ext, x, y, nx, ny, alpha, beta, b
real(kind=8):: g0, g1, g2, g3, g4
real(kind=8), parameter:: WaveNum=1.d0
complex :: XJ
real(kind=8):: arg1,arg2, space
real(kind=8):: kx,ky,CFun,SFun, AngFreq
real(kind=8):: u,u_x,u_y
real(kind=8):: Phi,Phi_x,Phi_y,tmp1,tmp2,tmp3,tmp4,r,theta,r_scale,const1,const2

XJ=(0,1)
! subroutine begin
!------------------------------------------------------------------------------
!     Test case for Pual's paper
      arg1= pi*x
      arg2= pi*y
      CFun= cos(arg1)*sin(arg2); SFun = sin(arg1)*cos(arg2)
   
      select case (Edge_Num)
      case(1)
         g_ext = alpha * ( sin(arg1)*sin(arg2) ) &
               + beta * pi * (CFun * nx + SFun * ny)
      case(3)
         g_ext = alpha * ( sin(arg1)*sin(arg2) ) &
               + beta * pi * (CFun * nx + SFun * ny)
      case(2)
         g_ext = alpha * ( sin(arg1)*sin(arg2) ) &
               + beta * pi * (CFun * nx + SFun * ny)
      case(4)
         g_ext = alpha * ( sin(arg1)*sin(arg2) ) &
               + beta * pi * (CFun * nx + SFun * ny)
      end select
!------------------------------------------------------------------------------
! Test case for hp convergence of continuous material:
! 1domain,4domain,9domain,16domain


!   arg1= pi*x
!   arg2= pi*y
!CFun= -sin(arg1); SFun = 0.d0
!
! select case (Edge_Num)
! case(1)
!   g_ext = alpha * (1 + cos(arg1) ) &
!         + beta * pi * (CFun * nx + SFun * ny)
!case(3)
!   g_ext = alpha * (1 + cos(arg1) ) &
!         + beta *  pi * (CFun * nx + SFun * ny)
!case(2)
!   g_ext = alpha * (1 + cos(arg1) ) &
!         + beta *  pi * (CFun * nx + SFun * ny)
!case(4)
!    g_ext = alpha * (1 + cos(arg1) ) &
!          + beta *  pi * (CFun * nx + SFun * ny)
!end select
!arg1= (0.1*pi)*x
!arg2= (0.1*pi)*y
!CFun= -sin(arg1)*cos(arg2); SFun=-cos(arg1)*sin(arg2)
!
! select case (Edge_Num)
! case(1)
!   g_ext = alpha * cos(arg1) * cos(arg2) &
!         + beta * pi * (CFun * nx + SFun * ny)
!case(3)
!   g_ext = alpha * cos(arg1) * cos(arg2) &
!         + beta *  pi * (CFun * nx + SFun * ny)
!case(2)
!   g_ext = alpha * cos(arg1) * cos(arg2) &
!         + beta *  pi * (CFun * nx + SFun * ny)
!case(4)
!    g_ext = alpha * cos(arg1) * cos(arg2) &
!          + beta *  pi * (CFun * nx + SFun * ny)
!end select

!select case (Edge_Num)
!case(1)
!g_ext = alpha * cos(arg1) * cos(arg2) &
!+ beta * (0.1*pi) * (CFun * nx + SFun * ny)
!case(3)
!g_ext = alpha * cos(arg1) * cos(arg2) &
!+ beta * (0.1*pi) * (CFun * nx + SFun * ny)
!case(2)
!g_ext = alpha * cos(arg1) * cos(arg2) &
!+ beta * (0.1*pi) * (CFun * nx + SFun * ny)
!case(4)
!g_ext = alpha * cos(arg1) * cos(arg2) &
!+ beta * (0.1*pi) * (CFun * nx + SFun * ny)
!end select

!   arg1= 2*pi*x
!   arg2= 2*pi*y
! CFun= -sin(arg1)*cos(arg2); SFun=-cos(arg1)*sin(arg2)
!
! select case (Edge_Num)
! case(1)
!   g_ext = alpha * cos(arg1) * cos(arg2) &
!         + beta * 2 * pi * (CFun * nx + SFun * ny)
!case(3)
!   g_ext = alpha * cos(arg1) * cos(arg2) &
!         + beta * 2 * pi * (CFun * nx + SFun * ny)
!case(2)
!   g_ext = alpha * cos(arg1) * cos(arg2) &
!         + beta * 2 * pi * (CFun * nx + SFun * ny)
!case(4)
!    g_ext = alpha * cos(arg1) * cos(arg2) &
!          + beta * 2 * pi * (CFun * nx + SFun * ny)
!end select
!------------------------------------------------------------------------------
! CFun= ( -exp(y)*sin(x))/(exp(y)*cos(x)+4)
! SFun= (  exp(y)*cos(x))/(exp(y)*cos(x)+4)
!
! select case (Edge_Num)
! case(1)
!   g_ext = alpha * log(exp(y)*cos(x)+4) &
!         + beta * (CFun * nx + SFun * ny)
!case(3)
!   g_ext = alpha * log(exp(y)*cos(x)+4) &
!         + beta * (CFun * nx + SFun * ny)
!case(2)
!   g_ext = alpha * log(exp(y)*cos(x)+4) &
!         + beta * (CFun * nx + SFun * ny)
!case(4)
!    g_ext = alpha * log(exp(y)*cos(x)+4) &
!          + beta * (CFun * nx + SFun * ny)
!end select
!-------------------------------------------------------------------
! Test case for Material Discontinuity at x=0
!   Discontinue interface case
!    if ( x .le. 0.0) then
!    arg1= pi*x
!    arg2= pi*y
!    CFun= cos(arg1)*sin(arg2); SFun= sin(arg1)*cos(arg2)
!
!    select case (Edge_Num)
!    case(1)
!        g_ext = alpha * sin(arg1) * sin(arg2) &
!        + beta * pi * (CFun * nx + SFun * ny)
!    case(3)
!        g_ext = alpha * sin(arg1) * sin(arg2) &
!        + beta * pi * (CFun * nx + SFun * ny)
!    case(2)
!        g_ext = alpha * sin(arg1) * sin(arg2) &
!        + beta * pi * (CFun * nx + SFun * ny)
!    case(4)
!        g_ext = alpha * sin(arg1) * sin(arg2) &
!        + beta * pi * (CFun * nx + SFun * ny)
!    end select
!
!    else if (x .ge. 0.0) then
!    arg1= 2.0*pi*x
!    arg2= pi*y
!    CFun= 2.0 * cos(arg1)*sin(arg2);
!    SFun= sin(arg1)*cos(arg2)
!
!    select case (Edge_Num)
!    case(1)
!        g_ext = alpha * 5.0/4.0 * sin(arg1) * sin(arg2) &
!        + beta * 5.0/4.0 * pi * (CFun * nx + SFun * ny)
!    case(3)
!        g_ext = alpha * 5.0/4.0 * sin(arg1) * sin(arg2) &
!        + beta * 5.0/4.0 * pi * (CFun * nx + SFun * ny)
!    case(2)
!        g_ext = alpha * 5.0/4.0 * sin(arg1) * sin(arg2) &
!        + beta * 5.0/4.0 * pi * (CFun * nx + SFun * ny)
!    case(4)
!        g_ext = alpha * 5.0/4.0 * sin(arg1) * sin(arg2) &
!        + beta * 5.0/4.0 * pi * (CFun * nx + SFun * ny)
!    end select
!    endif
!
!------------------------------------------------------------------------------
! Test case for Circular Domain
! Phi = (x**2.0 + y**2.0 -4) * (x**3.0-3*x*(y**2.0))
! Phi_x = 5*(x**4.0) - 6*(x**2.0)*(y**2.0) - 3*(y**4.0) &
!       - 12 * (x**2.0) + 12 * (y**2.0)
! Phi_y = -4 * (x**3.0) * y - 12 * x * (y**3.0) &
!          + 24 * x * y
!    select case (Edge_Num)
!    case(1)
!    g_ext = alpha * Phi &
!         + beta *  (Phi_x * nx + Phi_y * ny)
!    case(3)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(2)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(4)
!    g_ext = alpha * Phi &
!          + beta * (Phi_x * nx + Phi_y * ny)
!    end select
!
!------------------------------------------------------------------------------
! Test case for RingMode with material discontinuity at r=2.5
!    n=3
!    r=sqrt(x**2.0+y**2.0) ; theta=acos(x/r)
!    if (y .lt. 0) theta =-theta!atan2(y,x)
!    const1 = (2.0/5.0)**3.0
!    const2 = (2.0/5.0)**2.0
!    if (r .le. 2.5d0) then
!!            Phi = 7 * (r-2) *cos(n*theta)
!!            Phi_x = 7 * cos(n*theta) *cos(theta) &
!!                  + n * (7.0-14.0/r) * sin(n*theta) * sin(theta)
!!            Phi_y = 7 * cos(n*theta) * sin(theta) &
!!                  - n * (7.0-14.0/r) * sin(n*theta) * cos(theta)
!!           Phi =  (r**2-4)
!!           Phi_x = 2*r*cos(theta)
!!           Phi_y = 2*r*sin(theta)
!!          Phi = 7 * (r-2)!(1.584)*r
!!          Phi_x = 7 * cos(theta) !(1.584)*cos(theta)
!!          Phi_y = 7 * sin(theta) !(1.584)*sin(theta)
!
!!    Phi = 1.76 * (r**2.0-4) * cos(n*theta)
!!    Phi_x = 1.76 * (2*r) * cos(n*theta) * cos(theta) &
!!          + n * 1.76 * (r-(4/r)) * sin(n*theta) *sin(theta)
!!    Phi_y = 1.76 * (2*r) * cos(n*theta) * sin(theta) &
!!          - n * 1.76 * (r-(4/r)) * sin(n*theta) * cos(theta)
!
!Phi = ((r**2.0/4.d0)-1) * cos(n*theta)
!Phi_x = (r/2.d0) * cos(n*theta) * cos(theta) &
!+ n * ((r/4.d0)-(1/r)) * sin(n*theta) *sin(theta)
!Phi_y = (r/2.d0) * cos(n*theta) * sin(theta) &
!- n * ((r/4.d0)-(1/r)) * sin(n*theta) * cos(theta)
!
!!     Phi = (r**2.0-4) * cos(n*theta)
!!     Phi_x = (2*r) * cos(n*theta) * cos(theta) &
!!           + n * (r-(4/r)) * sin(n*theta) *sin(theta)
!!     Phi_y = (2*r) * cos(n*theta) * sin(theta) &
!!           - n * (r-(4/r)) * sin(n*theta) * cos(theta)
!!         Phi = (-r**2.0+5.4*r-7.2)*cos(n*theta)
!!         Phi_x = (-2*r+5.4)*cos(n*theta)*cos(theta) &
!!               + n * (-r + 5.4 -(7.2/r) ) *sin(n*theta)*sin(theta)
!!         Phi_y = (-2*r+5.4)*cos(n*theta)*sin(theta) &
!!               - n * (-r+5.4-(7.2/r))*sin(n*theta)*cos(theta)
!!        Phi = r*cos(n*theta)
!!        Phi_x = cos(n*theta)*cos(theta) &
!!              + n*sin(n*theta)*sin(theta)
!!        Phi_y = cos(n*theta)*sin(theta) &
!!              - n*sin(n*theta)*cos(theta)
!        !r_scale = 2.0 * r / 5.0
!        !Phi   = (r_scale**3.0-2*r_scale**2.0) * cos(n*theta)
!!(0.4*(x**2.0)+0.4*(y**2.0)-2) * ((x**3.0)-3*x*(y**2.0))
!
!        !Phi_x =   (3*const1*r**2.0-4*const2*r)*cos(n*theta)*cos(theta) &
!        !+ n*(const1*r**2.0-2*const2*r)*sin(n*theta)*sin(theta)
!!( (3*r**2.0-4*r) * cos(n*theta) * cos(theta) ) &
!!+ ( n*(r**2.0-2*r) * sin(n*theta) * sin(theta) )
!
!!(2*(x**4.0)-2.4*(x**2.0)*(y**2.0)-1.2*(y**4.0)-6*(x**2.0)+6*(y**2.0))
!        !Phi_y = (3*const1*r**2.0-4*const2*r)*cos(n*theta)*sin(theta) &
!        !- n*(const1*r**2.0-2*const2*r)*sin(n*theta)*cos(theta)
!
!!( (3*r**2.0-4*r) * cos(n*theta) * sin(theta) ) &
!!              - ( n*(r**2.0-2*r) * sin(n*theta) * cos(theta) )
!
!!(-1.6*(x**3.0)*y - 4.8*x*(y**3.0) + 12*x*y)
!
!    select case (Edge_Num)
!    case(1)
!    g_ext = alpha * Phi &
!         + beta *  (Phi_x * nx + Phi_y * ny)
!    case(3)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(2)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(4)
!    g_ext = alpha * Phi &
!          + beta * (Phi_x * nx + Phi_y * ny)
!    end select
!
!    else
!!           Phi = (3*r-4) * cos(n*theta)
!!           Phi_x = 3 * cos(n*theta) * cos(theta) &
!!                 + n * (3.0 - 4.0/r) * sin(n*theta) * sin(theta)
!!           Phi_y = 3 * cos(n*theta) * sin(theta) &
!!                 - n * (3.0 - 4.0/r) * sin(n*theta) * cos(theta)
!!          Phi = (3*r-4)
!!          Phi_x = 3 * cos(theta)
!!          Phi_y = 3 * sin(theta)
!!         Phi = 45 * (-r**2.0+5.4*r-7.2)
!!         Phi_x = 45 * (-2*r+5.4)*cos(theta)
!!         Phi_y = 45 * (-2*r+5.4)*sin(theta)
!
!!    Phi = (-r**2.0+5.4*r-3.29)*cos(n*theta)
!!    Phi_x = (-2*r+5.4)*cos(n*theta)*cos(theta) &
!!          + n * (-r+5.4-(3.29/r)) * sin(n*theta)*sin(theta)
!!    Phi_y = (-2*r+5.4)*cos(n*theta)*sin(theta) &
!!          - n * (-r+5.4-(3.29/r)) * sin(n*theta) * cos(theta)
!Phi = 11.25 * (-r**2.0+5.4*r-7.2)*cos(n*theta)
!Phi_x = 11.25 * (-2*r+5.4)*cos(n*theta)*cos(theta) &
!+ n * 11.25 * (-r+5.4-(7.2/r)) * sin(n*theta)*sin(theta)
!Phi_y = 11.25 * (-2*r+5.4)*cos(n*theta)*sin(theta) &
!- n * 11.25 * (-r+5.4-(7.2/r)) * sin(n*theta) * cos(theta)
!
!!     Phi = 45 * (-r**2.0+5.4*r-7.2)*cos(n*theta)
!!     Phi_x = 45 * (-2*r+5.4)*cos(n*theta)*cos(theta) &
!!           + n * 45 * (-r+5.4-(7.2/r)) * sin(n*theta)*sin(theta)
!!     Phi_y = 45 * (-2*r+5.4)*cos(n*theta)*sin(theta) &
!!           - n * 45 * (-r+5.4-(7.2/r)) * sin(n*theta) * cos(theta)
!!        Phi = (-r**2.0+5.4*r-7.2)*cos(n*theta)
!!        Phi_x = (-2*r+5.4)*cos(n*theta)*cos(theta) &
!!              + n * (-r + 5.4 -(7.2/r) ) *sin(n*theta)*sin(theta)
!!        Phi_y = (-2*r+5.4)*cos(n*theta)*sin(theta) &
!!              - n * (-r+5.4-(7.2/r))*sin(n*theta)*cos(theta)
!
!!        Phi = r*cos(n*theta)
!!        Phi_x = cos(n*theta)*cos(theta) &
!!              + n*sin(n*theta)*sin(theta)
!!        Phi_y = cos(n*theta)*sin(theta) &
!!              - n*sin(n*theta)*cos(theta)
!        !r_scale = 2.0 * r / 5.0
!        !Phi   = (-r_scale**3.0) * cos(n*theta)
!!(0.5) * ((x**3.0)-3*x*(y**2.0))!(1.75*r**2-0.5*r**3) * cos(n*theta)
!        !Phi_x = (-3*const1*r**2.0)*cos(n*theta)*cos(theta) &
!        !    - n*(const1*r**2.0)*sin(n*theta)*sin(theta)
!
!!( (5.4*r**2.0-8*r) * cos(n*theta) * cos(theta) ) &
!!+ ( n*(1.8*r**2.0-4*r) * sin(n*theta) * sin(theta) )
!
!!(0.5) * (3*(x**2.0)-3*(y**2.0))
!              !( (3.5*r-1.5*r**2) * cos(n*theta) * cos(theta) ) &
!              !+ ( 3*(1.75*r-0.5*r**2) * sin(n*theta) * sin(theta) )
!        !Phi_y = (-3*const1*r**2.0)*cos(n*theta)*sin(theta) &
!        !    + n*(const1*r**2.0)*sin(n*theta)*cos(theta)
!
!!( (5.4*r**2.0-8*r) * cos(n*theta) * sin(theta) ) &
!!              - ( n*(1.8*r**2.0-4*r) * sin(n*theta) * cos(theta) )
!
!!(0.5) * (-6*x*y)
!              !( (3.5*r-1.5*r**2) * cos(n*theta) * sin(theta) ) &
!              !- ( 3*(1.75*r-0.5*r**2) * sin(n*theta) * cos(theta) )
!    select case (Edge_Num)
!    case(1)
!    g_ext = alpha * Phi &
!    + beta *  (Phi_x * nx + Phi_y * ny)
!    case(3)
!    g_ext = alpha * Phi &
!    + beta * (Phi_x * nx + Phi_y * ny)
!    case(2)
!    g_ext = alpha * Phi &
!    + beta * (Phi_x * nx + Phi_y * ny)
!    case(4)
!    g_ext = alpha * Phi &
!    + beta * (Phi_x * nx + Phi_y * ny)
!    end select
!endif
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

   ! test cases u= sin(pi*x)*cos(pi*y)
   ! General BCs
!  kx=1.0*pi; ky=1.0*pi;
!  AngFreq=sqrt(kx**2+ky**2)
!  arg1 = kx*x + ky*y - AngFreq*t;
!
!  g0 =                     SFun * alpha + CFun * (nx*kx+ny*ky)*beta*b
!
!CFun= 1!-sin(arg1);
!SFun= 0
!select case (Edge_Num)
!case(1)
!g_ext = alpha * x &
!+ beta * (CFun * nx + SFun * ny)
!case(3)
!g_ext = alpha * x &
!+ beta  * (CFun * nx + SFun * ny)
!case(2)
!g_ext = alpha * x &
!+ beta  * (CFun * nx + SFun * ny)
!case(4)
!g_ext = alpha * x &
!+ beta  * (CFun * nx + SFun * ny)
!end select



!
!arg1= pi*x
!CFun= 0!-2*x!-sin(arg1);
!SFun= -2*y!0
!select case (Edge_Num)
!case(1)
!g_ext = alpha * 3-y**2 &
!+ beta * (CFun * nx + SFun * ny)
!case(3)
!g_ext = alpha * 3-y**2 &
!+ beta  * (CFun * nx + SFun * ny)
!case(2)
!g_ext = alpha * 3-y**2 &
!+ beta  * (CFun * nx + SFun * ny)
!case(4)
!g_ext = alpha * 3-y**2  &
!+ beta  * (CFun * nx + SFun * ny)
!end select
!

!CFun= 0!-sin(arg1); 
!SFun= 0
!select case (Edge_Num)
!case(1)
!g_ext = alpha * 2 &
!+ beta * (CFun * nx + SFun * ny)
!case(3)
!g_ext = alpha * 2 &
!+ beta  * (CFun * nx + SFun * ny)
!case(2)
!g_ext = alpha * 2 &
!+ beta  * (CFun * nx + SFun * ny)
!case(4)
!g_ext = alpha * 2  &
!+ beta  * (CFun * nx + SFun * ny)
!end select
!
!arg1= pi*x
!CFun= -sin(arg1); SFun=0
!select case (Edge_Num)
!case(1)
!g_ext = alpha * cos(arg1) &
!+ beta * pi * (CFun * nx + SFun * ny)
!case(3)
!g_ext = alpha * cos(arg1) &
!+ beta * pi * (CFun * nx + SFun * ny)
!case(2)
!g_ext = alpha * cos(arg1) &
!+ beta * pi * (CFun * nx + SFun * ny)
!case(4)
!g_ext = alpha * cos(arg1)  &
!+ beta * pi * (CFun * nx + SFun * ny)
!end select


!u = (x**2+y**2-5*sqrt(x**2+y**2)+6)*exp(XJ*2*ATANH(x/y))
!
!u_x = (2*x - (5 * x / sqrt(x**2+y**2)))* exp(XJ*2*ATANH(x/y)) &
!    + ((XJ*2*y)/(x**2+y**2)) &
!    * ( x**2 + y**2 - 5 * sqrt(x**2+y**2) + 6 ) &
!    * exp(XJ*2*ATANH(x/y))
!
!u_y = (2*y - (5 * y / sqrt(x**2+y**2)))* exp(XJ*2*ATANH(x/y)) &
!    - ((XJ*2*x)/(x**2+y**2)) &
!    * ( x**2 + y**2 - 5 * sqrt(x**2+y**2) + 6 ) &
!    * exp(XJ*2*ATANH(x/y))
!
!select case (Edge_Num)
!case(1)
!g_ext = alpha * u + beta * (u_x * nx + u_y * ny)
!case(3)
!g_ext = alpha * u + beta * (u_x * nx + u_y * ny)
!case(2)
!g_ext = alpha * u + beta * (u_x * nx + u_y * ny)
!case(4)
!g_ext = alpha * u + beta * (u_x * nx + u_y * ny)
!end select


 !Circular domain case

!    n=3
!    r=sqrt(x**2.0+y**2.0) ; theta=atan2(y,x)
!    if (r .le. 2.5d0) then
!        Phi   = (r**3-8) * cos(n*theta)/(7.625)
!        Phi_x = (3 * r**2 * cos(n*theta) * cos(theta) &
!              + 3 * (r**3-8) * sin(n*theta) * sin(theta) /r) &
!              / (7.625)
!
!        Phi_y = (3 * r**2 * cos(n*theta) * sin(theta) &
!              - 3 * (r**3-8) * sin(n*theta) * cos(theta) /r) &
!              / (7.625)
!    select case (Edge_Num)
!    case(1)
!    g_ext = alpha * Phi &
!         + beta *  (Phi_x * nx + Phi_y * ny)
!    case(3)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(2)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(4)
!    g_ext = alpha * Phi &
!          + beta * (Phi_x * nx + Phi_y * ny)
!    end select
!
!    else
!    Phi   = (r**3-27) * cos(n*theta)/(-11.375)
!    Phi_x = (3 * r**2 * cos(n*theta) * cos(theta) &
!          + 3 * (r**3-27) * sin(n*theta) * sin(theta) /r) &
!          / (-11.375)
!
!    Phi_y = (3 * r**2 * cos(n*theta) * sin(theta) &
!          - 3 * (r**3-27) * sin(n*theta) * cos(theta) /r) &
!          / (-11.375)
!    select case (Edge_Num)
!    case(1)
!    g_ext = alpha * Phi &
!    + beta *  (Phi_x * nx + Phi_y * ny)
!    case(3)
!    g_ext = alpha * Phi &
!    + beta * (Phi_x * nx + Phi_y * ny)
!    case(2)
!    g_ext = alpha * Phi &
!    + beta * (Phi_x * nx + Phi_y * ny)
!    case(4)
!    g_ext = alpha * Phi &
!    + beta * (Phi_x * nx + Phi_y * ny)
!    end select
!endif

!if ( abs(x) .le. 1E-15) then
!if (y .ge. 0.d0 ) theta = pi/2.d0
!if (y .lt. 0.d0 ) theta =-pi/2.d0
!else
!theta = atan(y/x)
!if (x .le. 0.d0) theta = theta + pi
!endif

!tmp1 =  24 * r * x**4.0 &
!     -  18 * r**3.0 * x**2.0
!tmp2 =  9 * r * y**4.0 &
!     -  12 * r**3.0 * y**2.0

!tmp3 =  24 * r * x**3.0 * y &
!     -  18 * r**3.0 * x * y
!tmp4 =  -9 * r * y**3.0 * x &
!     +  12 * r**3.0 * y * x

!    tmp1 = (4*r**3.0-9*r**2.0+4*r) &
!         * cos(theta) * sin(n*theta)
!    tmp2 = - n*(r**3.0-3*r**2.0+2*r) &
!         * cos(n*theta) * sin(theta)
!    tmp3 = (4*r**3.0-9*r**2.0+4*r)  &
!         * sin(n*theta) * sin(theta)
!    tmp4 = n * (r**3.0-3*r**2.0+2*r) &
!         * cos(n*theta) * cos(theta)


!!    if(x == 0.0 .and. y==0.0) then
!!    write(*,*)"r",r
!!    write(*,*)"theta",theta
!!    write(*,*)"tmp1",tmp1
!!    write(*,*)"tmp3",tmp3
!!    write(*,*)"phi",Phi
!!    write(*,*)"phi_x",Phi_x
!!    write(*,*)"phi_y",Phi_y
!!    endif
!
!    select case (Edge_Num)
!    case(1)
!    g_ext = alpha * Phi &
!         + beta *  (Phi_x * nx + Phi_y * ny)
!    case(3)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(2)
!    g_ext = alpha * Phi &
!         + beta * (Phi_x * nx + Phi_y * ny)
!    case(4)
!    g_ext = alpha * Phi &
!          + beta * (Phi_x * nx + Phi_y * ny)
!    end select

!
!


!  ! stationary solution
!  arg1=WaveNum*pi; arg2=WaveNum*pi*sqrt(2.d0)
!  g0 = sin(arg1*x) * sin(arg1*y) * cos(arg2*t)
!  g1 = - arg2* sin(arg1*x) * sin(arg1*y) * sin(arg2*t)       
!  g2 = - arg2**2 * sin(arg1*x) * sin(arg1*y) * cos(arg2*t) 
!  g3 =   arg2**3 * sin(arg1*x) * sin(arg1*y) * sin(arg2*t)
!  g4 =   arg2**4 * sin(arg1*x) * sin(arg1*y) * cos(arg2*t)
  
!  ! travelling wave with Dirichlet BCs
!  arg1 = pi*(x+y-sqrt(2.0)*t)
!  g0 =                       sin(arg1)
!  g1 =  (-sqrt(2.0)*pi)    * cos(arg1)
!  g2 = -(-sqrt(2.0)*pi)**2 * sin(arg1)
!  g3 = -(-sqrt(2.0)*pi)**3 * cos(arg1)
!  g4 =  (-sqrt(2.0)*pi)**4 * sin(arg1)

!  ! stationary wave with general BCs
!  arg1 = WaveNum*pi; arg2= WaveNum*pi*sqrt(2.0)
!  space = alpha * cos(arg1*x) * cos(arg1*y) &
!        - beta * b * arg1* (  nx*sin(arg1*x)*cos(arg1*y) &
!                            + ny*cos(arg1*x)*sin(arg1*y) )
!  
!  g0 =              space * cos(arg2*t)
!  g1 = -(arg2)    * space * sin(arg2*t)
!  g2 = -(arg2)**2 * space * cos(arg2*t)
!  g3 =  (arg2)**3 * space * sin(arg2*t)
!  g4 =  (arg2)**4 * space * cos(arg2*t)
!  
  ! travelling wave with general Robin BCs
!  kx=1.0*pi; ky=1.0*pi; 
!  AngFreq=sqrt(kx**2+ky**2)
!  arg1 = kx*x + ky*y - AngFreq*t; 
!  CFun= cos(arg1); SFun=sin(arg1)
!
!  g0 =                     SFun * alpha + CFun * (nx*kx+ny*ky)*beta*b
!  g1 = (-AngFreq)**1 * (   CFun * alpha - SFun * (nx*kx+ny*ky)*beta*b )
!  g2 = (-AngFreq)**2 * ( - SFun * alpha - CFun * (nx*kx+ny*ky)*beta*b )
!  g3 = (-AngFreq)**3 * ( - CFun * alpha + SFun * (nx*kx+ny*ky)*beta*b )
!  g4 = (-AngFreq)**4 * (   SFun * alpha + CFun * (nx*kx+ny*ky)*beta*b )
!  
!!  ! standing wave & circular pulse with neumann BC
!  g0 = 0.0
!  g1 = 0.0
!  g2 = 0.0 
!  g3 = 0.0
!  g4 = 0.0
!  
!  select case (ks)
!  case (1)
!    g_ext = g0
!  case (2,3)
!    g_ext = g0 + dt/2.0*g1 + (dt**2)/8.0*g2
!  case (4)  
!    g_ext = g0 + dt*g1 + (dt**2)/2.0*g2 &
!          + (dt**3)/4.d0*g3 + (dt**4)/16.d0*g4
!  end select
  return
end function


!subroutine compute_PBC2
!use MD2D_Grid
!use State_Var
!use RKN_Var
!implicit none
!integer:: i,j
!real(kind=8):: Bq_loc, x_coord, y_coord
!real(kind=8):: nx, ny, alpha, beta, b_loc, g_ext
!
!  do DDK=1,TotNum_DM
!     ND1 = PolyDegN_DM(1,DDK)
!     ND2 = PolyDegN_DM(2,DDK)
!     
!     Edge_Num =1 
!     select case (BC_Type(Edge_Num,DDK))
!
!     case(1,2,3)     
!        do j=0,ND2
!           x_coord= x1(0,j,DDK); 
!           y_coord= x2(0,j,DDK);
!           nx     = NorVec_x1(j,Edge_Num,DDK)
!           ny     = NorVec_x2(j,Edge_Num,DDK)
!           alpha  = bnd_alpha(j,Edge_Num,DDK)
!           beta   =  bnd_beta(j,Edge_Num,DDK) 
!           b_loc  =     bEdge(j,Edge_Num,DDK)
!           
!           Bq_loc = alpha * q_tmp(0,j,DDK) &
!                  + beta * b_loc * ( nx * dvdx(0,j,DDK) &
!                                   + ny * dvdy(0,j,DDK) )
!           PBC(j,Edge_Num,DDK) = Bq_loc &
!                  - g_ext(x_coord,y_coord,time,dt,&
!                          nx,ny,alpha,beta,b_loc,ks) 
!        enddo 
!        
!!     case (0)
!!       DDK_Connect  = DM_Connect(1,DDK)
!!       Edge_Connect = DM_Connect(2,DDK)
!!       ND1_Connect  = PolyDegN_DM(1,DDK_Connect)
!!       ND2_Connect  = PolyDegN_DM(2,DDK_Connect)
!!       
!!       do j=0,ND2
!!           x_coord= x1(0,j,DDK); 
!!           y_coord= x2(0,j,DDK);
!!           nx     = NorVec_x1(j,Edge_Num,DDK)
!!           ny     = NorVec_x2(j,Edge_Num,DDK)
!!           alpha  = bnd_alpha(j,Edge_Num,DDK)
!!           beta   =  bnd_beta(j,Edge_Num,DDK) 
!!           b_loc  =     bEdge(j,Edge_Num,DDK)
!!           
!!           Bq_loc = alpha * q_tmp(0,j,DDK) &
!!                  + beta * b_loc * ( nx * dvdx(0,j,DDK) &
!!                                   + ny * dvdy(0,j,DDK) )
!!                                  
!!!           alpha  = bnd_alpha(j,Edge_Connect,DDK_Connect)
!!!           beta   =  bnd_beta(j,Edge_Connect,DDK_Connect)                     
!!!            g_loc = alpha * q_tmp(0,ND2_Connect,DDK_Connect) 
!!                     
!!           PBC(j,Edge_Num,DDK) = Bq_loc &
!!                  - g_ext(x_coord,y_coord,time,dt,&
!!                          nx,ny,alpha,beta,b_loc,ks)  
!!                          
!!        enddo      
!     end select  
!     
!     Edge_Num=2
!     select case (BC_Type(Edge_Num,DDK))
!
!     case(1,2,3)     
!        do j=0,ND2
!           x_coord= x1(ND1,j,DDK); 
!           y_coord= x2(ND1,j,DDK);
!           nx     = NorVec_x1(j,Edge_Num,DDK)
!           ny     = NorVec_x2(j,Edge_Num,DDK)
!           alpha  = bnd_alpha(j,Edge_Num,DDK)
!           beta   =  bnd_beta(j,Edge_Num,DDK) 
!           b_loc  =     bEdge(j,Edge_Num,DDK) 
!           
!           Bq_loc = alpha * q_tmp(ND1,j,DDK) &
!                  + beta * b_loc * ( nx * dvdx(ND1,j,DDK) &
!                                   + ny * dvdy(ND1,j,DDK) )
!           PBC(j,Edge_Num,DDK) = Bq_loc &
!                  - g_ext(x_coord,y_coord,time,dt,&
!                          nx,ny,alpha,beta,b_loc,ks)            
!        enddo        
!     end select 
!     
!     Edge_Num=3
!     select case (BC_Type(Edge_Num,DDK))
!
!     case(1,2,3)     
!        do i=0,ND1
!           x_coord= x1(i,0,DDK); 
!           y_coord= x2(i,0,DDK);
!           nx     = NorVec_x1(i,Edge_Num,DDK)
!           ny     = NorVec_x2(i,Edge_Num,DDK)
!           alpha  = bnd_alpha(i,Edge_Num,DDK)
!            beta  =  bnd_beta(i,Edge_Num,DDK)
!            b_loc =     bEdge(i,Edge_Num,DDK) 
!                       
!           Bq_loc = alpha *q_tmp(i,0,DDK) &
!                  + beta * b_loc * ( nx * dvdx(i,0,DDK) &
!                                   + ny * dvdy(i,0,DDK) )
!           
!           PBC(i,Edge_Num,DDK) = Bq_loc &
!                  - g_ext(x_coord,y_coord,time,dt,&
!                          nx,ny,alpha,beta,b_loc,ks) 
!      
!        enddo 
!
!     end select      
!        
!     Edge_Num=4
!     select case (BC_Type(Edge_Num,DDK))
!
!     case(1,2,3)
!        do i=0,ND1
!           x_coord= x1(i,ND2,DDK); 
!           y_coord= x2(i,ND2,DDK);
!           nx     = NorVec_x1(i,Edge_Num,DDK)
!           ny     = NorVec_x2(i,Edge_Num,DDK)
!           alpha  = bnd_alpha(i,Edge_Num,DDK)
!           beta   =  bnd_beta(i,Edge_Num,DDK) 
!           b_loc  =     bEdge(i,Edge_Num,DDK)          
!           Bq_loc = alpha *q_tmp(i,ND2,DDK) &
!                  + beta * b_loc * ( nx * dvdx(i,ND2,DDK) &
!                                   + ny * dvdy(i,ND2,DDK) )
!           PBC(i,Edge_Num,DDK) = Bq_loc &
!                  - g_ext(x_coord,y_coord,time,dt,&
!                          nx,ny,alpha,beta,b_loc,ks)
!        enddo   
!     end select      
!  enddo 
!
!end subroutine 
