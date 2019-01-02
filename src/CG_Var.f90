      module CG_Var
      implicit none

      ! parameter for CG
!      real(kind=8), save  :: beta_cg, alpha_cg
!      real(kind=8), save  :: rsold, rsnew, pAp, error_cg,rsnew1

      real(kind=8), save, allocatable::  x_init(:,:,:)
!      real(kind=8), save, allocatable::  Ap(:,:,:)
!      real(kind=8), save, allocatable::  Ax(:,:,:)
!      real(kind=8), save, allocatable::  r(:,:,:)
!      real(kind=8), save, allocatable::  p(:,:,:),z(:,:,:)
      contains
      !---------------------------------------------------------
      subroutine Init_CG(DegMax,TotNum_DM,iterNum)
      implicit none

      integer:: DegMax(2), TotNum_DM,iterNum
      integer:: ND1, ND2, N_max
      integer:: ierr

      N_max=maxval(DegMax(1:2))
      ND1=DegMax(1); ND2=DegMax(2)

      allocate( x_init(0:ND1,0:ND2,1:TotNum_DM), &
!                       Ax(0:ND1,0:ND2,1:TotNum_DM), &
!                       Ap(0:ND1,0:ND2,1:TotNum_DM), &
!                        r(0:ND1,0:ND2,1:TotNum_DM), &
!                        p(0:ND1,0:ND2,1:TotNum_DM), &
!                        z(0:ND1,0:ND2,1:TotNum_DM), &
                                    stat=ierr)
      if (ierr .ne. 0) then
         write(*,*)'Cannot allocate memory for CG variables'
         write(*,*)'at CG_Var module, Init_CG'
         write(*,*)'Abort!'
         stop
      endif

      x_init=0.0d0

!      Ax=0.0d0;
!      rsnew=0.d0; rsold=0.d0; r=0.0d0; 
!      p=0.0d0; Ap=0.0d0; z=0.0d0

      
!      ND1=maxval(DegMax(1:2))
!      allocate(cg_dqdx(0:ND1,0:ND1),cg_dqdy(0:ND1,0:ND1), stat=ierr)
!      if (ierr .ne. 0) then
!         write(*,*)' Message from State_Var.f90'
!         write(*,*)' Cannot allocate memory for variable dqdx'
!         write(*,*)' Abort!'
!         stop
!      endif
!      
!      cg_dqdx=0.d0; cg_dqdy=0.d0
      
      return
      end subroutine Init_CG
      end module CG_Var
