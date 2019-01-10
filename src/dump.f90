      subroutine dumpopt(opt,level,name3)

      use MD2D_Grid

      implicit none
      integer :: i,j,level,ii
      real(kind=8) :: opt(1)
      character*3 name3
      character*7 filename

      write(filename,fmt="(A3)") name3
      open(999,file=filename)

      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level);ND2=PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               ii = (DDK-1)*(ND2+1)*(ND1+1) + (i+1) + j*(ND1+1)
               write(999,1005) opt(ii)
            enddo ! i
         enddo ! j
      enddo ! DDK
      close(999)

1005  format(9e24.15)
      end subroutine
