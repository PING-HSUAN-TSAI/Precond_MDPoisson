      function glsc2(LD1,LD2,x,y,level)
      use MD2D_Grid    ! MD2D_Grid.f90
      implicit none

      integer:: i,j, LD1, LD2, level
      real(kind=8) :: x(0:LD1,0:LD2,1), y(0:LD1,0:LD2,1)
      real(kind=8) :: ds
      real(kind=8) :: glsc2

      ds = 0.0

      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level)
         ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               ds = ds + (x(i,j,DDK) * y(i,j,DDK))
            enddo
         enddo
      enddo

      glsc2 = ds
      return
      end function
!=============================================================================================
      subroutine ADD3S2(A,B,C,C1,C2,LD1,LD2,level)
      use MD2D_Grid
      implicit none

      integer:: i,j, LD1, LD2, level
      real(kind=8) :: C1, C2
      real(kind=8) :: A(0:LD1,0:LD2,1), B(0:LD1,0:LD2,1), C(0:LD1,0:LD2,1)
      
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level)
         ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               A(i,j,DDK) = C1 * B(i,j,DDK) + C2 * C(i,j,DDK)
            enddo
         enddo
      enddo
      return 
      end subroutine
!=============================================================================================
      subroutine ADD3S3(A,B,C,D,C1,C2,C3,LD1,LD2,level)
      use MD2D_Grid
      implicit none

      integer:: i,j, LD1, LD2, level
      real(kind=8) :: C1, C2, C3
      real(kind=8) :: A(0:LD1,0:LD2,1), B(0:LD1,0:LD2,1)
      real(kind=8) :: C(0:LD1,0:LD2,1), D(0:LD1,0:LD2,1)

      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level)
         ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               A(i,j,DDK) = C1 * B(i,j,DDK) + C2 * C(i,j,DDK) & 
                          + C3 * D(i,j,DDK)
            enddo
         enddo
      enddo
      return
      end subroutine
!=============================================================================================
      subroutine COPY(A,B,LD1,LD2,level)
      use MD2D_Grid
      implicit none
      
      integer:: i,j, LD1, LD2, level
      real(kind=8) :: A(0:LD1,0:LD2,1), B(0:LD1,0:LD2,1)

      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level)
         ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               A(i,j,DDK) = B(i,j,DDK)
            enddo
         enddo
      enddo

      return
      end subroutine
!=============================================================================================
      subroutine ADD2S2(A,B,C1,LD1,LD2,level)
      use MD2D_Grid
      implicit none

      integer:: i,j, LD1, LD2, level
      real(kind=8) :: C1
      real(kind=8) :: A(0:LD1,0:LD2,1), B(0:LD1,0:LD2,1)

      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level)
         ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               A(i,j,DDK) = A(i,j,DDK) + C1 * B(i,j,DDK)
            enddo
         enddo
      enddo

      return
      end subroutine
!=============================================================================================
      function GLMAX(A,LD1,LD2,level)
      use MD2D_Grid
      implicit none
      
      integer:: i,j, LD1, LD2, level
      real(kind=8) :: A(0:LD1,0:LD2,1)
      real(kind=8) :: TMAX,TMP
      real(kind=8) :: glmax
      
      TMAX = -99.0e20
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level)
         ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               TMAX = max(TMAX,A(i,j,DDK))
            enddo
         enddo
      enddo

      TMP = TMAX

      GLMAX = TMP
      
      return
      end function
!=============================================================================================
      function GLAMAX(A,LD1,LD2,level)
      use MD2D_Grid
      implicit none

      integer:: i,j, LD1, LD2, level
      real(kind=8) :: A(0:LD1,0:LD2,1)
      real(kind=8) :: TMAX,TMP
      real(kind=8) :: glamax

      TMAX = 0.0
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level)
         ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               TMAX = max(TMAX,abs(A(i,j,DDK)))
            enddo
         enddo
      enddo

      TMP = TMAX

      GLAMAX = abs(TMP)

      return
      end function
!=============================================================================================
      subroutine chk_amax(s3,A,LD1,LD2,level)
      use MD2D_Grid
      implicit none 

      integer:: i,j, LD1, LD2, level
      character*3 s3 
      real(kind=8) :: amx, glamax
      real(kind=8) :: A(0:LD1,0:LD2,1)

      amx = glamax(A,LD1,LD2,level)

      write(10,*)"check2:",s3,amx

      return

      end subroutine
