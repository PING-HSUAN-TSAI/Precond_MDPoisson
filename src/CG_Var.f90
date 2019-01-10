      module CG_Var
      implicit none

      ! Variable for iterative solver
      real(kind=8), save, allocatable::  x_init(:,:,:)

      contains
!---------------------------------------------------------
      subroutine Init_CG(DegMax,TotNum_DM,iterNum)
      implicit none

      integer:: DegMax(2), TotNum_DM,iterNum
      integer:: ND1, ND2, N_max
      integer:: ierr

      N_max=maxval(DegMax(1:2))
      ND1=DegMax(1); ND2=DegMax(2)

      allocate( x_init(0:ND1,0:ND2,1:TotNum_DM),stat=ierr)

      if (ierr .ne. 0) then
         write(*,*)'Cannot allocate memory for CG variables'
         write(*,*)'at CG_Var module, Init_CG'
         write(*,*)'Abort!'
         stop
      endif

      x_init=0.0d0

      
      return
      end subroutine Init_CG
      end module CG_Var
