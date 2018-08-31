subroutine Init_Physical_Grid_2D(mgl)
use constants
use MD2D_Grid
use Legendre
  implicit none

  ! Declare local arguments
  ! Boundary grid for setting physical space grid by transfinite blending
  ! Function : Need to deallocate before leaving this subroutine
  
  real(kind=8), allocatable :: Edge_x(:,:,:), Edge_y(:,:,:)
  integer ierr
  integer N, i, j, k, l
  integer level, mgl
  
  !Linear mapping parameters
  real(kind=8) :: x_start, x_end
  real(kind=8) :: y_start, y_end
  !Circular arc mapping parameters
  real(kind=8) :: x_cen, y_cen, radius
  real(kind=8) :: v1x, v1y, Lv1, v2x, v2y, Lv2
  real(kind=8) :: angle_ref, cross_ref, cross, angle, t  
  
  !-----------------------------------------------------
  integer:: ND_Max
  real(kind=8) :: xi1, xi2 
  
  write(*,*)'mgl',mgl
  
  !  Start subroutine
  !  Find the max of degree in all directions
!  ND1=PolyDegN_Max(1); ND2=PolyDegN_Max(2);
  ND1=PolyDegN_DM(1,1,mgl); ND2=PolyDegN_DM(2,1,mgl);
    write(*,*)ND1,ND2
!  ND_Max=max(ND1,ND2,ND3)
  
  allocate( Edge_x(0:ND1,1:12,mgl), &
            Edge_y(0:ND2,1:12,mgl), stat= ierr)
  
  if ( ierr .ne. 0 ) then
     write(*,*) 'Grid2D_Pack.f90:'
     write(*,*) 'Can not allocate Edge_x Edge_y andin Grid2D_Pack'
     write(*,*) 'Abort!'
     stop
  endif
  
! Construct Computational Grids by Transfinite Blending Mapping
! Step 1: Construct Points on Egdes (Totally 4 Edges)
! Step 2: Construct Points on Surfaces (Total 6 Surfaces) based on the Edges
!
! Start looping over all domain
open(1127,file='fine_grid.text')
open(1128,file='coarse_grid.text')
  do l = mgl,1,-1
    do DDK=1,TotNum_DM
     
    ND1=PolyDegN_DM(1,DDK,l);
    ND2=PolyDegN_DM(2,DDK,l);
    write(*,*)ND1,ND2

     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     !
     ! Set up grid points on Edges 
     !
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

     ! Edge 1: v1 (-1,-1) to v2 ( 1,-1)
     ! ^^^^^^ 
     !...................... Side 1 ..............................!
     select case(DMShapeType(1,DDK)) ! Side 1
     case (0) ! straight side linear transform

     x_start = DM_Vertex(1,1,DDK) 
     y_start = DM_vertex(1,2,DDK) 
     
     x_end  = DM_Vertex(2,1,DDK)
     y_end  = DM_vertex(2,2,DDK) 
     
     ND = PolyDegN_DM(1,DDK,l)

     Edge_x(0:ND,1,l) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,1,l) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
          
     case (1) ! circular arc
     side_num = 1
	 x_cen	= DMShapePara(1,side_num,DDK)
	 y_cen	= DMShapePara(2,side_num,DDK)
	 radius = DMShapePara(3,side_num,DDK)

	 v1x = DM_Vertex(1,1,DDK) - x_cen
	 v1y = DM_Vertex(1,2,DDK) - y_cen
	 Lv1 = sqrt(v1x**2+v1y**2)
	 v2x = DM_Vertex(2,1,DDK) - x_cen
	 v2y = DM_Vertex(2,2,DDK) - y_cen
	 Lv2 = sqrt(v2x**2+v2y**2)
	 angle_ref = acos(v1x/Lv1)
	 cross_ref = v1y
	 if (cross_ref .lt. 0.0d0) angle_ref = 2.0d0*pi - angle_ref
	 cross = v1x*v2y - v1y*v2x
	 angle = acos((v1x*v2x+v1y*v2y)/(Lv1*Lv2))
	 if (cross .lt. 0.0d0) angle=-angle
	 if (abs(angle) .lt. 1.0d-6) angle=2.0d0*pi
  
     ND=PolyDegN_DM(1,DDK,l)
	 do i=0,ND
	    t	= (1.0d0+LGLCoord(i,ND))/2.0d0*angle + angle_ref
	    Edge_x(i,1,l)= x_cen + radius * cos(t)
	    Edge_y(i,1,l)= y_cen + radius * sin(t)
	 end do    
     
     end select ! Side 1
     
     !-------------------------------------------------------------!     
     ! Edge 2: v2 ( 1,-1) to v3 ( 1, 1)
     ! ^^^^^^ 
     select case(DMShapeType(2,DDK)) ! Side 2
     case (0) ! straight side linear transform
     x_start = DM_Vertex(2,1,DDK)
     y_start = DM_vertex(2,2,DDK) 
     
     x_end  = DM_Vertex(3,1,DDK)
     y_end  = DM_vertex(3,2,DDK) 
     
     ND=PolyDegN_DM(2,DDK,l)
     
     Edge_x(0:ND,2,l) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,2,l) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
     
     case (1) ! Circular arr
	 side_num = 2
	 x_cen	= DMShapePara(1,side_num,DDK)
	 y_cen	= DMShapePara(2,side_num,DDK)
	 radius = DMShapePara(3,side_num,DDK)

	 v1x = DM_Vertex(2,1,DDK) - x_cen
	 v1y = DM_Vertex(2,2,DDK) - y_cen
	 Lv1 = sqrt(v1x**2+v1y**2)
	 v2x = DM_Vertex(3,1,DDK) - x_cen
	 v2y = DM_Vertex(3,2,DDK) - y_cen
	 Lv2 = sqrt(v2x**2+v2y**2)
	 angle_ref = acos(v1x/Lv1)
	 cross_ref = v1y
	 if (cross_ref.lt.0.0d0) angle_ref = 2.0d0*pi - angle_ref
	 cross = v1x*v2y - v1y*v2x
	 angle = acos((v1x*v2x+v1y*v2y)/(Lv1*Lv2))
	 if (cross.lt.0.0d0) angle=-angle
	 if (abs(angle).lt.1.0d-6) angle=2.0d0*pi

     ND=PolyDegN_DM(2,DDK,l)
	 do i=0,ND
	    t	= (1.0d0+LGLCoord(i,ND))/2.0d0*angle + angle_ref
	    Edge_x(i,2,l)= x_cen + radius * cos(t)
	    Edge_y(i,2,l)= y_cen + radius * sin(t)
	 end do     
          
     end select ! Side 2
     !-------------------------------------------------------------!     
     ! Edge 3: v4 (-1, 1) to v3 ( 1, 1) 
     ! ^^^^^^ 
     select case(DMShapeType(3,DDK)) ! Side 3
     case (0) ! straight side linear transform
     x_start = DM_Vertex(4,1,DDK)
     y_start = DM_vertex(4,2,DDK) 
     
     x_end  = DM_Vertex(3,1,DDK)
     y_end  = DM_vertex(3,2,DDK) 
     
     
     ND=PolyDegN_DM(1,DDK,l)
     
     Edge_x(0:ND,3,l) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,3,l) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
     
     case (1) ! circular arc
	 side_num = 3
	 x_cen	= DMShapePara(1,side_num,DDK)
	 y_cen	= DMShapePara(2,side_num,DDK)
	 radius = DMShapePara(3,side_num,DDK)

	 v1x = DM_Vertex(4,1,DDK) - x_cen
	 v1y = DM_Vertex(4,2,DDK) - y_cen
	 Lv1 = sqrt(v1x**2+v1y**2)
	 v2x = DM_Vertex(3,1,DDK) - x_cen
	 v2y = DM_Vertex(3,2,DDK) - y_cen
	 Lv2 = sqrt(v2x**2+v2y**2)
	 angle_ref = acos(v1x/Lv1)
	 cross_ref = v1y
	 if (cross_ref .lt. 0.0d0) angle_ref = 2.0d0*pi - angle_ref
	 cross = v1x*v2y - v1y*v2x
	 angle = acos((v1x*v2x+v1y*v2y)/(Lv1*Lv2))
	 if (cross .lt. 0.0d0) angle=-angle
	 if (abs(angle) .lt. 1.0d-6) angle=2.0d0*pi

     ND=PolyDegN_DM(1,DDK,l)
	 do i=0,ND
	    t	= (1.0d0+LGLCoord(i,ND))/2.0d0*angle + angle_ref
	    Edge_x(i,3,l)= x_cen + radius * cos(t)
	    Edge_y(i,3,l)= y_cen + radius * sin(t)
	 end do

     !end select ! Finish side 3 grid a     
     
     end select ! Side 3
     !-------------------------------------------------------------!     
     ! Edge 4: v1 (-1,-1) to v4 (-1, 1)
     ! ^^^^^^ 
     select case(DMShapeType(4,DDK)) ! Side 4
     case (0) ! straight side linear transform
     x_start = DM_Vertex(1,1,DDK)
     y_start = DM_vertex(1,2,DDK) 
     
     x_end  = DM_Vertex(4,1,DDK)
     y_end  = DM_vertex(4,2,DDK) 
     
     ND=PolyDegN_DM(2,DDK,l)
     
     Edge_x(0:ND,4,l) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,4,l) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
 
     case (1) ! Circular arc
     side_num = 4
	 x_cen	= DMShapePara(1,side_num,DDK)
	 y_cen	= DMShapePara(2,side_num,DDK)
	 radius = DMShapePara(3,side_num,DDK)

	 v1x = DM_Vertex(1,1,DDK) - x_cen
	 v1y = DM_Vertex(1,2,DDK) - y_cen
	 Lv1 = sqrt(v1x**2+v1y**2)
	 v2x = DM_Vertex(4,1,DDK) - x_cen
	 v2y = DM_Vertex(4,2,DDK) - y_cen
	 Lv2 = sqrt(v2x**2+v2y**2)
	 angle_ref = acos(v1x/Lv1)
	 cross_ref = v1y
!	 write(*,*)angle_ref, cross_ref
	 if (cross_ref .lt. 0.0d0) angle_ref = 2.0d0*pi - angle_ref
	 cross = v1x*v2y - v1y*v2x
	 angle = acos((v1x*v2x+v1y*v2y)/(Lv1*Lv2))
!	 write(*,*)cross,angle 
	 if (cross .lt. 0.0d0) angle=-angle
	 if (abs(angle) .lt. 1.0d-6) angle=2.0d0*pi


     ND=PolyDegN_DM(2,DDK,l)
	 do i=0,ND
	    t	= (1.0d0+LGLCoord(i,ND))/2.0d0*angle + angle_ref
	    Edge_x(i,4,l)= x_cen + radius * cos(t)
	    Edge_y(i,4,l)= y_cen + radius * sin(t)
!	    write(*,*)Edge_x(i,4),Edge_y(i,4)
	 end do
!	    pause

     end select ! Side 4
     
     !-------------------------------------------------------------!
     
     !-------------------------------------------------------------!
     ! This completes the points on Edges
     !-------------------------------------------------------------!
     
     !=============================================================!
     
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     !
     ! Set up grid points on Surfaces 
     !
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     
     do j=0,ND2
        
        xi2 = LGLCoord(j,ND2)
        
        do i=0,ND1
           
           xi1 = LGLCoord(i,ND1)
           
           ! Original Two D version
           x1(i,j,DDK,l) = &
                 (1.d0-xi2)/2.d0 * Edge_x(i,1,l) + &
                 (1.d0+xi2)/2.d0 * Edge_x(i,3,l) + &
                 (1.d0-xi1)/2.d0 * Edge_x(j,4,l) + &
                 (1.d0+xi1)/2.d0 * Edge_x(j,2,l) - &
                 (1.d0-xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(1,1,DDK)- &
                 (1.d0+xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(2,1,DDK)- &
                 (1.d0+xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(3,1,DDK)- &
                 (1.d0-xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(4,1,DDK)

           x2(i,j,DDK,l) = &
                 (1.d0-xi2)/2.d0 * Edge_y(i,1,l) + &
                 (1.d0+xi2)/2.d0 * Edge_y(i,3,l) + &
                 (1.d0-xi1)/2.d0 * Edge_y(j,4,l) + &
                 (1.d0+xi1)/2.d0 * Edge_y(j,2,l) - &
                 (1.d0-xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(1,2,DDK)- &
                 (1.d0+xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(2,2,DDK)- &
                 (1.d0+xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(3,2,DDK)- &
                 (1.d0-xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(4,2,DDK)

        enddo ! enddo i
        
     enddo ! enddo j
     
     
     !-------------------------------------------------------------!
     ! This completes the points on Quadalateron
     !-------------------------------------------------------------!
     
     !=============================================================!
     
  enddo ! enddo DDK
  
  enddo ! l

  do DDK = 1, TotNum_DM
  ND1 = PolyDegN_DM(1,DDK,mgl)
  ND2 = PolyDegN_DM(2,DDK,mgl)
    do j = 0 ,ND2
      do i = 0, ND1
        write(1127,*)i,j,x1(i,j,DDK,mgl),x2(i,j,DDK,mgl)
      enddo
    enddo
  enddo
close(1127)

do DDK = 1, TotNum_DM
ND1 = PolyDegN_DM(1,DDK,mgl-1)
ND2 = PolyDegN_DM(2,DDK,mgl-1)
do j = 0 ,ND2
do i = 0, ND1
write(1128,*)i,j,x1(i,j,DDK,mgl-1),x2(i,j,DDK,mgl-1)
enddo
enddo
enddo
close(1128)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!  Deallocate variable Edge_x and Edge_y
  if ( allocated( Edge_x ) ) Deallocate( Edge_x, stat=ierr )
  if ( ierr .ne. 0 ) then
     write(*,*) 'Can not deallocate Bnd_x and in Grid'
     write(*,*) 'Abort'
     stop
  endif
  
  if ( allocated( Edge_y ) ) Deallocate( Edge_y, stat=ierr )
  if ( ierr .ne. 0 ) then
     write(*,*) 'Can not deallocate Bnd_y and in Grid'
     write(*,*) 'Abort'
     stop
  endif
  
  return
  
end subroutine Init_Physical_Grid_2D

subroutine Init_Metric(mgl)
  use Legendre
  use MD2D_Grid
  implicit none

  integer:: i, j, k
  integer:: l,mgl

  call alloc_mem_metric_variables(PolyDegN_DM(1,1,mgl),&
                                  PolyDegN_DM(2,1,mgl),&
                                  TotNum_DM,mgl)

  ! compute dx/dxi
  do l =mgl,1,-1
  do DDK=1,TotNum_DM
     ND1=PolyDegN_DM(1,DDK,l); 
     ND2=PolyDegN_DM(2,DDK,l); 
     
     dx1_dxi1(0:ND1,0:ND2,DDK,l)=&
          Matmul( Diff_xi1(0:ND1,0:ND1,ND1), x1(0:ND1,0:ND2,DDK,l))

     dx2_dxi1(0:ND1,0:ND2,DDK,l)=&
          Matmul( Diff_xi1(0:ND1,0:ND1,ND1), x2(0:ND1,0:ND2,DDK,l))

     dx1_dxi2(0:ND1,0:ND2,DDK,l)=&
          Matmul( x1(0:ND1,0:ND2,DDK,l) , Diff_xi2(0:ND2,0:ND2,ND2) )

     dx2_dxi2(0:ND1,0:ND2,DDK,l)=&
          Matmul( x2(0:ND1,0:ND2,DDK,l) , Diff_xi2(0:ND2,0:ND2,ND2) )

     ! compute Jacobin and dxi/dx 
     Jacobin(0:ND1,0:ND2,DDK,l) = &
        dx1_dxi1(0:ND1,0:ND2,DDK,l) * dx2_dxi2(0:ND1,0:ND2,DDK,l)   &
                              - &
        dx2_dxi1(0:ND1,0:ND2,DDK,l) * dx1_dxi2(0:ND1,0:ND2,DDK,l) 
     
     dxi1_dx1(0:ND1,0:ND2,DDK,l) = &
          dx2_dxi2(0:ND1,0:ND2,DDK,l)/Jacobin(0:ND1,0:ND2,DDK,l)

     dxi2_dx1(0:ND1,0:ND2,DDK,l) = &
         -dx2_dxi1(0:ND1,0:ND2,DDK,l)/Jacobin(0:ND1,0:ND2,DDK,l) 

     dxi1_dx2(0:ND1,0:ND2,DDK,l) = &
         -dx1_dxi2(0:ND1,0:ND2,DDK,l)/Jacobin(0:ND1,0:ND2,DDK,l) 

     dxi2_dx2(0:ND1,0:ND2,DDK,l) = &
          dx1_dxi1(0:ND1,0:ND2,DDK,l)/Jacobin(0:ND1,0:ND2,DDK,l)

  enddo ! DDK 
  enddo ! l
!  write(10,*)'Jaociban'
!do DDK=1,TotNum_DM
!     ND1=PolyDegN_DM(1,DDK,mgl); 
!     ND2=PolyDegN_DM(2,DDK,mgl); 
!do j=0,ND2
!do i=0,ND1
!write(10,*)i,j,DDK,Jacobin(i,j,DDK,mgl),dx1_dxi1(i,j,DDK,mgl),dx2_dxi2(i,j,DDK,mgl),&
!         dxi1_dx1(i,j,DDK,mgl),dxi2_dx2(i,j,DDK,mgl)
!enddo
!enddo
!enddo

  ! we need to check wheather the computed jacobin is positive
  do l = mgl,1,-1
  do DDK=1,TotNum_DM 
     ND1=PolyDegN_DM(1,DDK,l); 
     ND2=PolyDegN_DM(2,DDK,l); 

     do j=0,ND2
        do i=0,ND1 
           if (Jacobin(i,j,DDK,l) .le. 0.d0) then
              write(*,*)'Message from Metric_Pack.f90'
              write(*,*)'Jacobin is not positive:',i,j,Jacobin(i,j,DDK,l)
              write(*,*)'Abort!'
           endif
        enddo
     enddo
  enddo ! DDK
  enddo !l

  ! compute grid distortion for determining dt
  do l = mgl,1,-1
  do DDK=1,TotNum_DM

     ND1=PolyDegN_DM(1,DDK,l); 
     ND2=PolyDegN_DM(2,DDK,l);

     ! compute local grid distance dxi1 and dxi2
     dxi1(0,0:ND2,DDK,l)=LGLCoord(1,ND1)-LGLCoord(0,ND1)
     do i=1,ND1-1
        dxi1(i,0:ND2,DDK,l)=(LGLCoord(i+1,ND1)-LGLCoord(i-1,ND1))/2.d0
     enddo
     dxi1(ND1,0:ND2,DDK,l)=LGLCoord(ND1,ND1)-LGLCoord(ND1-1,ND1)
     
     
     dxi2(0:ND1,0,DDK,l)=LGLCoord(1,ND2)-LGLCoord(0,ND2)
     do j=1,ND2-1
        dxi2(0:ND1,j,DDK,l)=(LGLCoord(j+1,ND2)-LGLCoord(j-1,ND2))/2.d0
     enddo
     dxi2(0:ND1,ND2,DDK,l)=LGLCoord(ND2,ND2)-LGLCoord(ND2-1,ND2)


     ! compute local grid distortion vector Cal X dot Cal X = dtrans

     ! Let 
     !           | grad xi1 | = ( |dxi1/dx1|, |dxi1/dx2| )
     !           | grad xi2 | = ( |dxi2/dx1|, |dxi2/dx2| )
     !           
     !           X = | grad xi1 |/dxi1 + | grad xi2 |/dxi2
     !           dtrans = sqrt( X dot X ) 

     dtrans(0:ND1,0:ND2,DDK,l) = sqrt(&
       ( abs(dxi1_dx1(0:ND1,0:ND2,DDK,l)) / dxi1(0:ND1,0:ND2,DDK,l) + &
         abs(dxi2_dx1(0:ND1,0:ND2,DDK,l)) / dxi2(0:ND1,0:ND2,DDK,l) ) ** 2 &
      +( abs(dxi1_dx2(0:ND1,0:ND2,DDK,l)) / dxi1(0:ND1,0:ND2,DDK,l) + &
         abs(dxi2_dx2(0:ND1,0:ND2,DDK,l)) / dxi2(0:ND1,0:ND2,DDK,l) ) ** 2 )

  enddo ! DDK
  enddo ! l

  write(10,*)'Complete Initializing Metric Variables'
  return 
  
end subroutine Init_Metric

subroutine Init_Normal_Vector(mgl)
!use legendre 
  use MD2D_Grid
  implicit none

  integer i, j, l,mgl
  ! subroutine begin
  ! allocate memory for Normal Vector Veriables

!  call alloc_mem_norvec(maxval(PolyDegN_Max(1:2)),TotNum_DM)
  call alloc_mem_norvec(PolyDegN_DM(1,1,mgl),TotNum_DM,mgl)


  ! Compute Normal Vectors of all the surface of a cube
  do l = mgl,1,-1
  do DDK=1, TotNum_DM

     ND1=PolyDegN_DM(1,DDK,l)
     ND2=PolyDegN_DM(2,DDK,l)
     
     do Edge_Num=1,4
        
        select case(Edge_Num)
        case(1) ! n = - grad xi_2
              
           NorVec_mg(0:ND1,Edge_Num,DDK,l) = sqrt( &
                dxi2_dx1(0:ND1,0,DDK,l)**2 + dxi2_dx2(0:ND1,0,DDK,l)**2 )
           
           NorVec_x1(0:ND1,Edge_Num,DDK,l) =  - &
                ( dxi2_dx1(0:ND1,0,DDK,l) / NorVec_mg(0:ND1,Edge_Num,DDK,l) )
           
           NorVec_x2(0:ND1,Edge_Num,DDK,l) =  - &
                ( dxi2_dx2(0:ND1,0,DDK,l) / NorVec_mg(0:ND1,Edge_Num,DDK,l) )
                
           JacNorVec(0:ND1,Edge_Num,DDK,l) = Jacobin(0:ND1, 0 ,DDK,l) &
                                         * NorVec_mg(0:ND1,Edge_Num,DDK,l)
                   
        case(2) ! n = + grad xi_1

           NorVec_mg(0:ND2,Edge_Num,DDK,l) = sqrt( &
                dxi1_dx1(ND1,0:ND2,DDK,l)**2 + dxi1_dx2(ND1,0:ND2,DDK,l)**2 )
           
           NorVec_x1(0:ND2,Edge_Num,DDK,l) =   &
                ( dxi1_dx1(ND1,0:ND2,DDK,l) / NorVec_mg(0:ND2,Edge_Num,DDK,l) )
           
           NorVec_x2(0:ND2,Edge_Num,DDK,l) =   &
                ( dxi1_dx2(ND1,0:ND2,DDK,l) / NorVec_mg(0:ND2,Edge_Num,DDK,l) )
                
           JacNorVec(0:ND2,Edge_Num,DDK,l) = Jacobin(ND1,0:ND2,DDK,l) &
                                         * NorVec_mg(0:ND2,Edge_Num,DDK,l)
           
        case(3) ! n = + grad xi_2

           NorVec_mg(0:ND1,Edge_Num,DDK,l) = sqrt( &
                dxi2_dx1(0:ND1,ND2,DDK,l)**2 + dxi2_dx2(0:ND1,ND2,DDK,l)**2 )
           
           NorVec_x1(0:ND1,Edge_Num,DDK,l) =    &
                ( dxi2_dx1(0:ND1,ND2,DDK,l) / NorVec_mg(0:ND1,Edge_Num,DDK,l) )
           
           NorVec_x2(0:ND1,Edge_Num,DDK,l) =    &
                ( dxi2_dx2(0:ND1,ND2,DDK,l) / NorVec_mg(0:ND1,Edge_Num,DDK,l) )
                
           JacNorVec(0:ND1,Edge_Num,DDK,l) = Jacobin(0:ND1,ND2,DDK,l) &
                                         * NorVec_mg(0:ND1,Edge_Num,DDK,l)     

        case(4) ! n = -grad xi_1

           NorVec_mg(0:ND2,Edge_Num,DDK,l) = sqrt( &
                dxi1_dx1(0,0:ND2,DDK,l)**2 + dxi1_dx2(0,0:ND2,DDK,l)**2 )
           
           NorVec_x1(0:ND2,Edge_Num,DDK,l) = - &
                ( dxi1_dx1(0,0:ND2,DDK,l) / NorVec_mg(0:ND2,Edge_Num,DDK,l) )
           
           NorVec_x2(0:ND2,Edge_Num,DDK,l) = - &
                ( dxi1_dx2(0,0:ND2,DDK,l) / NorVec_mg(0:ND2,Edge_Num,DDK,l) )

           JacNorVec(0:ND2,Edge_Num,DDK,l) = Jacobin( 0 ,0:ND2,DDK,l) &
                                         * NorVec_mg(0:ND2,Edge_Num,DDK,l)
!        do i=0,ND1
!           write(*,*)i,NorVec_x1(i,Edge_Num,DDK,l),NorVec_x2(i,Edge_Num,DDK,l)
!        enddo
        end select
       
!        do i=0,ND1
!           write(*,*)'JacNorvec',NorVec_mg(i,Edge_Num,DDK)
!        enddo
!        pause 
     enddo ! Edge_Num

  enddo ! DDK
  enddo ! l
  write(10,*)'Complete Initializing Normal Vector'

10000 format(2i4,3f15.7)
             
end subroutine Init_Normal_Vector


