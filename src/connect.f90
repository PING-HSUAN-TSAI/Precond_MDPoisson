      subroutine rdparam(UnitNum)

      use INPUT

      implicit none

      integer :: I, ERR, UnitNum 

!      OPEN (UNIT=9,FILE='pspm2d.rea',STATUS='OLD',iostat=ierr)
      READ(UnitNum,*)
      READ(UnitNum,*) 
      READ(UnitNum,*) 
      READ(UnitNum,*) NPARAM

      DO 20 I=1,NPARAM
         READ(UnitNum,*) PARAM(I)
!         write(*,*)PARAM(I)
   20   CONTINUE

      return 
      end subroutine
