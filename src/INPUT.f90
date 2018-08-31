      module INPUT
      implicit none 
      real(kind=8), save, allocatable::  PARAM(:)
      integer, save :: democase, NPARAM 
      contains
!-----------------------------------------------------
      subroutine setupinvar
      implicit none
      integer :: ierr

      allocate( PARAM(1:200),&
                stat=ierr)
      if (ierr .ne. 0) then
         write(*,*)'Cannot allocate memory for Input variables'
         write(*,*)'at INPUT module, setupinvar'
         write(*,*)'Abort!'
         stop
      endif

      PARAM=0.d0
      return
      end subroutine 
      end module
