      subroutine dumpopt(opt,level,name3)

      use MD2D_Grid

      implicit none
      integer :: i,j,level,ii,E
      real(kind=8) :: opt(1)
      character*3 name3
      character*7 filename

      write(filename,fmt="(A3)") name3
      open(999,file=filename)

      do E=1,TotNum_DM
         ND1=PolyDegN_DM(1,E,level);ND2=PolyDegN_DM(2,E,level)
         do j=0,ND2
            do i=0,ND1
               ii = (E-1)*(ND2+1)*(ND1+1) + (i+1) + j*(ND1+1)
               write(999,1005) opt(ii)
            enddo ! i
         enddo ! j
      enddo ! E
      close(999)

1005  format(9e24.15)
      end subroutine
