module constants
implicit none
real(kind=8), save :: pi, pi2

contains

subroutine init_constants
  implicit none
  
  pi=4.d0*datan(1.d0) 
  pi2=2.d0*pi
  return 
end subroutine init_constants

end module constants
