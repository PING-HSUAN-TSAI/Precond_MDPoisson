      subroutine Field_Final(level)
      use MD2D_Grid
      use CG_Var
      use State_Var
      implicit none
      
      integer::level
      
      if (level .eq. 1) then

         call VTK_output(level)
         call DAT_output(level)

      else if ( level .ge. 2) then
      
      endif
      end subroutine Field_Final
!!======================================================================
      subroutine VTK_output(level)
      use MD2D_Grid
      use CG_Var
      use State_Var
      implicit none 

      integer:: i,j,points,cells,cells_size
      integer:: indice1,indice2,indice3,indice4
      integer:: level

      points=0.0
      cells=0.0


      !-output field as .vtk
      do DDK=1,TotNum_DM
        points = points + (PolyDegN_DM(1,DDK,level)+1) &
                        * (PolyDegN_DM(2,DDK,level)+1)
        cells = cells + (PolyDegN_DM(1,DDK,level) &
                      *  PolyDegN_DM(2,DDK,level))
      enddo

      cells_size = cells + cells*4

      open (1112,file='Field.vtk')
      write(1112, fmt="(A26,1x)" )'# vtk DataFile Version 2.0'
      write(1112, fmt="(A10,1x)" )'Field Data'
      write(1112, fmt="(A5,1x)"  )'ASCII'
      write(1112, fmt="(A25,1x)" )'DATASET UNSTRUCTURED_GRID'
      write(1112, fmt="(A6,1x,I4,1x,A5)")'POINTS',points,'float'

      do DDK=1,TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level)
         ND2 = PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               write(1112,fmt="(f13.5,f13.5,f13.5)")x1(i,j,DDK,level),x2(i,j,DDK,level),0.d0
            enddo
         enddo
      enddo

      write(1112,fmt="(A5,1x,I4,1x,I5)")'CELLS',cells,cells_size
      do DDK=1,TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level)
         ND2 = PolyDegN_DM(2,DDK,level)
         do j=0,ND2-1
            do i=0,ND1-1
               indice1 = (i+(j*(ND1+1))) &
                       + ((DDK-1)*(ND1+1)*(ND2+1))
               indice2 = (i+1+(j*(ND1+1))) &
                       + ((DDK-1)*(ND1+1)*(ND2+1))
               indice3 = (i+ND1+1+(j*(ND1+1))) &
                       + ((DDK-1)*(ND1+1)*(ND2+1))
               indice4 = (i+ND1+2+(j*(ND1+1))) &
                       + ((DDK-1)*(ND1+1)*(ND2+1))
               write(1112,fmt="(I2,1x,I4,1x,I4,1x,I4,1x,I4)")4,indice1,indice2,indice4,indice3
            enddo
         enddo
      enddo

      write(1112,fmt="(A10,1x,I4)")'CELL_TYPES',cells
      do DDK=1,TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level)
         ND2 = PolyDegN_DM(2,DDK,level)
         do j=0,ND2-1
            do i=1,ND1
               write(1112,*)9
            enddo
         enddo
      enddo

      write(1112,fmt="(A10,1x,I4)")'POINT_DATA',points
      write(1112,fmt="(A19)")'SCALARS Field float'
      write(1112,fmt="(A20)")'LOOKUP_TABLE default'
      do DDK=1,TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level)
         ND2 = PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND1
               write(1112,fmt="(f13.5)")potent(i,j,DDK)
            enddo
         enddo
      enddo

      close(1112)
     !write(1112,fmt="(A5,A9,I4)")'FIELD','fieldData',1
      !write(1112,fmt="(A5,I2,I4,1x,A5)")'field',3,points,'float'
      !do DDK=1,TotNum_DM
      !do j=0,ND2
      !do i=0,ND1
      !write(1112,fmt="(3f13.5)")0.0,0.0,U_CG(i,j,DDK)
      !enddo
      !enddo
      !enddo

      !write(1112,fmt="(A10,1x,I4)")'POINT_DATA',points
      !write(1112,fmt="(A7,1x,A9,1x,A5)")'VECTORS','fieldData','float'
      !do DDK=1,TotNum_DM
      !do j=0,ND2
      !do i=0,ND1
      !write(1112,fmt="(3f13.5)")0.0,0.0,U_CG(i,j,DDK)
      !enddo
      !enddo
      !enddo
      !end of output Field.vtk
      !----------------------------------------------------------
      end subroutine VTK_output
!!======================================================================
      subroutine DAT_output(level)
      use MD2D_Grid
      use CG_Var
      use State_Var
      implicit none
      integer:: i,j,k,ND1Tot,ND2Tot
      integer ::level

      open(999,file='field.dat')
      write(999,*)'VARIABLES = "X", "Y", "Z", "U"'

      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level);ND2=PolyDegN_DM(2,DDK,level)
         write(999,fmt="(A7,I5,A4,I5,A4,I5,A9,1x)") 'ZONE I=', ND1+1, &
         ', J=', ND2+1, ', K=', 1,', F=POINT'
         do j=0,ND2
            do i=0,ND1
               write(999,1001) x1(i,j,DDK,level),x2(i,j,DDK,level), &
                               potent(i,j,DDK),potent(i,j,DDK)
            enddo ! i
         enddo ! j
      enddo ! DDK
      close(999)


!      open(555,file='error.dat')
!      write(555,*)'VARIABLES = "X", "Y", "Z", "U"'
!
!      do DDK=1,TotNum_DM
!         ND1=PolyDegN_DM(1,DDK,level);ND2=PolyDegN_DM(2,DDK,level)
!         write(555,fmt="(A7,I5,A4,I5,A4,I5,A9,1x)") 'ZONE I=', ND1+1, ', J=', ND2+1, ', K=', 1,', F=POINT'
!         do j=0,ND2
!            do i=0,ND1
!               write(555,1001) x1(i,j,DDK,level),x2(i,j,DDK,level),(v(i,j,DDK) - potent(i,j,DDK)),(v(i,j,DDK) - potent(i,j,DDK))
!            enddo ! i
!         enddo ! j
!      enddo ! DDK
!      close(555)
1000 format(9e24.15)
1001 format(9e24.15,1x)

      end subroutine 
