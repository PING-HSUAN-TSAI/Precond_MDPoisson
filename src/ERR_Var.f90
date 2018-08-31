      module ERR_Var
      implicit none

      real(kind=8), save, allocatable::  error_vec(:,:,:)
      real(kind=8), save, allocatable::  error_vec1(:,:,:)
      real(kind=8), save :: error, error_gradx, error_grady, error_max,error1,error_max1
      contains
!--------------------------------------------------------------------
      subroutine Init_ERR(DegMax,TotNum_DM)
      implicit none

      integer:: DegMax(2), TotNum_DM
      integer:: ND1, ND2
      integer:: ierr 
      
      ND1=DegMax(1); ND2=DegMax(2)
      allocate( error_vec(0:ND1,0:ND2,1:TotNum_DM), &
               error_vec1(0:ND1,0:ND2,1:TotNum_DM), &
                                           stat=ierr)
      
      if ( ierr .ne. 0 ) then
         write(*,*)'Cannot allocate memory for error variable '
         write(*,*)'Abort!'
         stop
      endif

      error_vec=0.0d0;  error_vec1=0.0d0;

      return
      end subroutine Init_ERR
      
      end module ERR_VAR
