      subroutine Error_Computing(mgl)
      use Legendre
      use State_Var
      use MD2D_Grid
      use ERR_Var
      use CG_Var
      implicit none
      
      integer :: i,j,k,error_loc,gridmax,grid
      integer :: mgl
      real(kind=8) :: U_CG(0:PolyDegN_DM(1,1,mgl),0:PolyDegN_DM(2,1,mgl),1:TotNum_DM) 
      !  real(kind=8) :: error_iter(1:count)
      
      gridmax=0.0
      error_max = 0.0
      error_max1 =0.0
      
      ! open(1107,file="cg_error.csv")
      !write(1107,*)'"iteration number", "cg_error"'
      open(93,file="error.text",position="append",action="write")
      
      !   write(*,*)'open error.text sucessfully'
      !open(92,file="error.csv",position="append",action="write")
      !write(92,*)'"grid number", "error"'
      
      !open(74,file='u_cg.text')
      !open(90,file="errorvec.text")
      !open(1125,file='dvdx_dFdx.text')
      !open(1126,file='dvdy_dFdy.text')
      !
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,mgl);ND2=PolyDegN_DM(2,DDK,mgl)
      !
      !    error_gradx = maxval(abs(dvdx(0:ND1,0:ND2,DDK)-dFdx(0:ND1,0:ND2,DDK)))
      !    error_grady = maxval(abs(dvdy(0:ND1,0:ND2,DDK)-dFdy(0:ND1,0:ND2,DDK)))
      !    do j=0,ND2
      !      do i=0,ND1
      !      write(1125,*)i,j,dvdx(i,j,DDK),dFdx(i,j,DDK),abs(dvdx(i,j,DDK)-dFdx(i,j,DDK))
      !      write(1126,*)i,j,dvdy(i,j,DDK),dFdy(i,j,DDK),abs(dvdy(i,j,DDK)-dFdy(i,j,DDK))
      !      enddo
      !    enddo
      !
      !    write(*,*)'error_gradx',error_gradx
      !    write(*,*)'error_grady',error_grady
      
      !      error_vec(0:ND1,0:ND2,DDK) = ( F(0:ND1,0:ND2,DDK)) * Jacobin(0:ND1,0:ND2,DDK) &
      !                                - laplace_v(0:ND1,0:ND2,DDK)
      
      
         error_vec(0:ND1,0:ND2,DDK) = ( v(0:ND1,0:ND2,DDK)) &
                                    - potent(0:ND1,0:ND2,DDK)
      
           do j=0,ND2
              do i=0,ND1
                  write(*,*) i,j,potent(i,j,DDK)
              enddo
          enddo
      !fmt="(2I3,1x,3E11.5e3)"
          
      
      !do j=0,ND2
      !do i=0,ND1
      !write(74,1000)x1(i,j,DDK),x2(i,j,DDK),U_CG(i,j,DDK),v(i,j,DDK)
      !!write(*,*)U_CG(i,j,DDK)
      !enddo
      !enddo
      
         grid  = max(ND1,ND2)
         gridmax = max(grid,gridmax)
      
         error = maxval(abs(error_vec(0:ND1,0:ND2,DDK)))
         error_max = max(error,error_max)
         write(10,1003)'Domain:',DDK,'Error:',error
      
      enddo
      
      !close(74)
      !close(90)
      
      write(10,fmt="(A10,es24.15)")' Max error:',error_max
      write(93,1001) gridmax,error_max
      close(93)
      
      open(91,file="out.text")
      write(91,fmt="(I3,1x,E11.5e3)") gridmax,error_max
      close(91)
      
      
      !    do k=1,count
      !    do DDK=1,TotNum_DM
      !
      !    ND1=PolyDegN_DM(1,DDK,1)
      !    ND2=PolyDegN_DM(2,DDK,1)
      !
      !    error_vec1(0:ND1,0:ND2,DDK) = ( v(0:ND1,0:ND2,DDK)) &
      !    - u_cg_iter(0:ND1,0:ND2,DDK,k)
      !
      !    error1 = maxval(abs(error_vec1(0:ND1,0:ND2,DDK)))
      !    error_max1 = MAX(error1,error_max1)
      !    enddo
      !    error_iter(k)=error_max1
      !    error_max1=0.0
      !    !write(*,*) error_iter(k)
      !    write(1107,1002) k,',',error_iter(k)
      !    enddo
      !
      !close(1007)
1000 format(9e23.15)      
1001 format(I3,es24.15)
!1002 format(I5,A2,e24.15)
1003 format(' ',A6,1x,I0.3,1x,A6,es24.15)
     return
     end subroutine Error_Computing
