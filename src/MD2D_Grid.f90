module MD2D_Grid
!--------------------------------------------------------------------------------
!
! This module contains functions and subroutines to initialize required 
! geometrical parameters for 3 dimensional multidomain pseudospectral 
! scheme
!
!--------------------------------------------------------------------------------
implicit none
!Declare Domain Variable
integer, save :: TotNum_DM
integer, save :: PolyDegN_max(2)
integer, save, allocatable      :: PolyDegN_DM(:,:,:) !(2,TotNum_DM)
real(kind=8), save, allocatable :: DM_Vertex(:,:,:) !(4,2,TotNum_DM)
real(kind=8), save, allocatable :: x1(:,:,:,:)      !(0:PolyDegN_Max(1:2),TotNum_DM)
real(kind=8), save, allocatable :: x2(:,:,:,:)      !(0:PolyDegN_Max(1:2),TotNum_DM)
integer,      save, allocatable :: DM_Connect(:,:,:)!(2,4,TotNum_DM)

integer :: DDK, ND, ND1, ND2, ND3, Surf_Num, Edge_Num  ! [ND: Degree Index, 
                                                       !DDK: Domain Index]

integer :: DDK_Connect, Surf_Connect, Patch_Type, Edge_Connect
integer :: ND1_Connect, ND2_Connect

integer :: side_num

! Domain Boundary Shape Variables
integer, save, allocatable:: DMShapeType(:,:)      ! (1:4,TotNum_DM)
real(kind=8), save, allocatable:: DMShapePara(:,:,:) ! (1:7,1:4,TotNum_DM)

! Metric Variables
real(kind=8), save, allocatable :: DM_Metric(:,:,:,:)
real(kind=8), save, allocatable :: Jacobin(:,:,:,:)   !((0:N_max)^3,TotNum_DM)
real(kind=8), save, allocatable :: dxi1(:,:,:,:)      !((0:N_max)^3,TotNum_DM)
real(kind=8), save, allocatable :: dxi2(:,:,:,:)
real(kind=8), save, allocatable :: dxi3(:,:,:,:)

! Declare variable dxi/dx
real(kind=8), save, allocatable :: dxi1_dx1(:,:,:,:)
real(kind=8), save, allocatable :: dxi1_dx2(:,:,:,:)
real(kind=8), save, allocatable :: dxi2_dx1(:,:,:,:)
real(kind=8), save, allocatable :: dxi2_dx2(:,:,:,:)

! Declare variable dx/dxi
real(kind=8), save, allocatable :: dx1_dxi1(:,:,:,:)
real(kind=8), save, allocatable :: dx1_dxi2(:,:,:,:)
real(kind=8), save, allocatable :: dx2_dxi1(:,:,:,:)
real(kind=8), save, allocatable :: dx2_dxi2(:,:,:,:)

! Declare Grid distortion variables
real(kind=8), save, allocatable :: dtrans(:,:,:,:)        !(0:N_max,TotNumDomain)

! Declare Normal vector Variables
real(kind=8), save, allocatable:: NorVec_x1(:,:,:,:)
real(kind=8), save, allocatable:: NorVec_x2(:,:,:,:)
real(kind=8), save, allocatable:: NorVec_mg(:,:,:,:)
real(kind=8), save, allocatable:: JacNorVec(:,:,:,:)

!Declare variable for multigrid level
integer , parameter:: MG_level = 2

! Description of variables
!================================================================================
! [Name] :: PolyDegN_DM(Index1,Index2)
!  ^^^^
!    [Size]    :: (1:2,1:TotNum_DM)
!     ^^^^
!    [Purpose] :: Store Polynomial Degree used in each domain
!     ^^^^^^^
!    [Detail]  :: Index2 is used to represent domain number
!     ^^^^^^
!                 Index1 is used to distinguish the number of Polynomial
!                 Degree used in x1 and x2 directions
!
!     Example :: PolyDeg_DM(1,DDK) stores the following:
!                    The degree of the polynomial used in the x1 direction 
!                    of Domain DDK.
!                  
!                PolyDeg_DM(2,DDK) stores the following:
!                    The degree of the polynomial used in the x2 direction
!                    of Domain DDK. 
! 
! [Name] :: DM_Vertex(Index1,Index2,Index3)
! ^^^^
!     Size    :: (1:5,1:2,1:TotNum_DM) 
!     ^^^^
!     Purpose :: To store the 4 physical coordinates of each domain  
!     ^^^^^^^
!     Detail  :: The first index indicate the four vertex 
!                the fifth one
!     ^^^^^^     is a copy of the first one for coding convience
!                The second index stands for x and y (1 for x, and 2 for y)
!                The third index is the number of the Domain.
!
! [Name] :: DM_Connect(index1,index2,index3) for 2D problem
! ^^^^
!     Size    ::  (3,4,1:TotNum_DM)
!     ^^^^
!     Purpose :: This variable is used to store the required info for  
!     ^^^^^^^    patching bc between two connecting domains
!              
!     Detail  :: Index3 is used to denoted the domain number 
!     ^^^^^^     
!                Index2 is used to denote the 4 edges: 
!                    1: xi_2= -1
!                    2: xi_1=  1
!                    3: xi_2=  1
!                    4: xi_1= -1
!                    
!
!                Index1 is used to denote the store info:
!                    1 for the number of the connected domain
!                    2 for the connecting edge
!                    3 for patching direction
!
!     Example :: DM_Connect(1,1,DDK) stores the following:
!                    The number of the connected domain on the 
!                    1st edge of domain DDK
!                  
!                DM_Connect(2,1,DDK) stores the following:
!                    The edge number of the connected domain on the 
!                    1st edge of domain DDK
!
!                DM_Connect(3,1,DDK) stores the following:
!                    The patching direction of the connected domain on the 
!                    1st edge of domain DDK. 1 for the same and -1 for opposite
!================================================================================

contains
  

  !==============================================================================
  subroutine Input_Domain_Grid_Parameters(mgl)
    implicit none
    !  Declare local variable
    !  integer :: DDK
    integer:: mgl
    
    !read in domain limit parameter
    open(20,file='domain.in',form='formatted',status='old')
    read(20, * ) !'====================== Input Parameters ==================='
    read(20, * )  TotNum_DM
    read(20, * )  PolyDegN_max(1)
    read(20, * )  PolyDegN_max(2)
    ! 
    !---allocate memory for Geometric Parameter
    call alloc_mem_domain_grid_variables(mgl)
    !                                        
    !---Read in Domain Geometric Parameter
    do DDK=1,TotNum_DM
!       read(20, * ) !'-----------------------------------------------------------'
!       read(20, * ) !'Parameters of Domain ### ----------------------------------'
!       read(20, * ) !'-----------------------------------------------------------'
!       read(20,200)  PolyDegN_DM(1,DDK), PolyDegN_DM(2,DDK)
!       read(20, * ) !'-----------------------------------------------------------' 
!       read(20,202)  DM_Vertex(1,1,DDK), DM_Vertex(1,2,DDK)
!       read(20,202)  DM_Vertex(2,1,DDK), DM_Vertex(2,2,DDK)
!       read(20,202)  DM_Vertex(3,1,DDK), DM_Vertex(3,2,DDK)
!       read(20,202)  DM_Vertex(4,1,DDK), DM_Vertex(4,2,DDK)
!       read(20, * ) !'-----------------------------------------------------------' 
!       read(20,204)  DM_Connect(1,1,DDK),DM_Connect(2,1,DDK) 
!       read(20,204)  DM_Connect(1,2,DDK),DM_Connect(2,2,DDK)
!       read(20,204)  DM_Connect(1,3,DDK),DM_Connect(2,3,DDK)
!       read(20,204)  DM_Connect(1,4,DDK),DM_Connect(2,4,DDK)
!       read(20, * ) !'-----------------------------------------------------------'
!       do side_num=1,4
!       read(20, * )  DMShapeType(side_num,DDK),   DMShapePara(1,side_num,DDK), &
!                     DMShapePara(2,side_num,DDK), DMShapePara(3,side_num,DDK), &
!                     DMShapePara(4,side_num,DDK), DMShapePara(5,side_num,DDK), &
!                     DMShapePara(6,side_num,DDK), DMShapePara(7,side_num,DDK)

       read(20, * ) !'-----------------------------------------------------------'
       read(20, * ) !'Parameters of Domain ### ----------------------------------'
       read(20, * ) !'-----------------------------------------------------------'
       read(20, * )  PolyDegN_DM(1,DDK,mgl), PolyDegN_DM(2,DDK,mgl)
       read(20, * ) !'-----------------------------------------------------------'
       read(20, * )  DM_Vertex(1,1,DDK), DM_Vertex(1,2,DDK)
       read(20, * )  DM_Vertex(2,1,DDK), DM_Vertex(2,2,DDK)
       read(20, * )  DM_Vertex(3,1,DDK), DM_Vertex(3,2,DDK)
       read(20, * )  DM_Vertex(4,1,DDK), DM_Vertex(4,2,DDK)
       read(20, * ) !'-----------------------------------------------------------'
       read(20, * )  DM_Connect(1,1,DDK),DM_Connect(2,1,DDK)
       read(20, * )  DM_Connect(1,2,DDK),DM_Connect(2,2,DDK)
       read(20, * )  DM_Connect(1,3,DDK),DM_Connect(2,3,DDK)
       read(20, * )  DM_Connect(1,4,DDK),DM_Connect(2,4,DDK)
       read(20, * ) !'-----------------------------------------------------------'
       do side_num=1,4
       read(20, * )  DMShapeType(side_num,DDK),   DMShapePara(1,side_num,DDK), &
                     DMShapePara(2,side_num,DDK), DMShapePara(3,side_num,DDK), &
                     DMShapePara(4,side_num,DDK), DMShapePara(5,side_num,DDK), &
                     DMShapePara(6,side_num,DDK), DMShapePara(7,side_num,DDK)
       enddo
    enddo
       read(20, * ) !'-----------------------------------------------------------'
    close(20)
    ! finish domain parameters computation
!    do DDK=1,TotNum_DM 
!       do Edge_Num=1,4 
!          Write(*,*)DM_Connect(1,Edge_Num,DDK),DM_Connect(2,Edge_Num,DDK)
!       enddo
!    enddo
!    pause
    return
200 format(1x,i5,1x,i5)
202 format(3x,f9.4,4x,f9.4)
204 format(1x,i5,1x,i5)
206 format(1x,i3,2x,f15.7,2x,f15.7,2x,f15.7,2x,f15.7,&
                 2x,f15.7,2x,f15.7,2x,f15.7)
  end subroutine Input_Domain_Grid_Parameters

  !==============================================================================
  Subroutine Geometric_Parameters_On_Screen(mgl)
    implicit none
    integer:: mgl,l
    
    ! Print Computation control parameters on screen
    write(10,* ) '====================== Input Parameters ==================='
    write(10,204)'ToTal Number of Domain : ', TotNum_DM
    write(10,205)'Maxium Degree of Polynomial :' , PolyDegN_max(1),&
                                                  PolyDegN_max(2)
    
    ! Set all approximation polynomial to same degree
    !   do DDK=2,TotNumDomain
    !      PolyDegreeN_Domain(DDK)=PolyDegreeN_Domain(1)
    !   enddo
    
    do DDK=1,TotNum_DM
       write(10, *)  '-----------------------------------------------------------'
       write(10, 203)'---------- Parameters of Domain', DDK,'--------------------'
       write(10, *)  '-----------------------------------------------------------'
       write(10, 207)  PolyDegN_DM(1,DDK,mgl), PolyDegN_DM(2,DDK,mgl)
       write(10, *)  '-----------------------------------------------------------'
       write(10, 201)  DM_Vertex(1,1,DDK), DM_Vertex(1,2,DDK)
       write(10, 201)  DM_Vertex(2,1,DDK), DM_Vertex(2,2,DDK)
       write(10, 201)  DM_Vertex(3,1,DDK), DM_Vertex(3,2,DDK)
       write(10, 201)  DM_Vertex(4,1,DDK), DM_Vertex(4,2,DDK)
       write(10, *)  '-----------------------------------------------------------'  
       write(10, 207)  DM_Connect(1,1,DDK),DM_Connect(2,1,DDK)
       write(10, 207)  DM_Connect(1,2,DDK),DM_Connect(2,2,DDK)
       write(10, 207)  DM_Connect(1,3,DDK),DM_Connect(2,3,DDK)
       write(10, 207)  DM_Connect(1,4,DDK),DM_Connect(2,4,DDK)
    enddo
       write(10, *)  '-----------------------------------------------------------'

    PolyDegN_DM(1,2:TotNum_DM,mgl)=PolyDegN_DM(1,1,mgl)
    PolyDegN_DM(2,2:TotNum_DM,mgl)=PolyDegN_DM(2,1,mgl)

      do l=mgl-1,1,-1
         do DDK = 1, TotNum_DM

            PolyDegN_DM(1,DDK,l) = PolyDegN_DM(1,DDK,mgl)/2
            PolyDegN_DM(2,DDK,l) = PolyDegN_DM(2,DDK,mgl)/2
            write(10,208)l,DDK,PolyDegN_DM(1,DDK,l),PolyDegN_DM(2,DDK,l)
         enddo
      enddo

    ! Check consistency of PolyDegN_DM and  PolyDegN_max
    do DDK=1,TotNum_DM
       if ( (PolyDegN_DM(1,DDK,mgl) .gt. PolyDegN_Max(1)) .or. &
            (PolyDegN_DM(2,DDK,mgl) .gt. PolyDegN_Max(2)) ) then
   
            write(*,*)'MD3D_Grid.f90:'
            write(*,*)'PolyDegN_DM > PolyDegN_Max for DDK=',DDK
            write(*,*)'Abort!'
            stop

       endif
    enddo

    return

!200 format(1x,i5)
201 format(3x,f8.4,4x,f8.4)
203 format(A32,1x,i3,2x,A22)
204 format(A32,1x,i4)
205 format(A32,1x,i4,1x,i4,1x,i4)
207 format(1x,i5,1x,i5)
208 format(' ',' ','l ',I0.3,' DDK: ',I0.3,1x,'Nc: ',I0.3,1x,I0.3)

  end subroutine Geometric_parameters_On_Screen

  !==============================================================================

  subroutine Init_patch_direction()
    implicit none
     
    DM_Connect(3,1:4,1:TotNum_DM)=1
    
    do DDK=1,TotNum_DM
     
       do Edge_Num=1,2
       
          select case (DM_Connect(2,Edge_Num,DDK)) ! side_number 1,2 of DDK
          case(1,2) ! connect to side_num=1,2: reverse
             DM_Connect(3,Edge_Num,DDK)=-1
             
          end select
          
       enddo
       
       do Edge_Num=3,4
          
          select case (DM_Connect(2,Edge_Num,DDK)) ! side_number 3,4 of DDK
          case(3,4) ! connect to side_num=3,4: reverse
             DM_Connect(3,Edge_Num,DDK)=-1
             
          end select
          
       enddo
       
    enddo
    
    return 
    
  end subroutine Init_patch_direction
  
  !==============================================================================

  subroutine Init_patch_direction2()
    implicit none
     
    DM_Connect(3,1:4,1:TotNum_DM)=1
    
    do DDK=1,TotNum_DM
     
       do Edge_Num=1,4
       
          if (DM_Connect(2,Edge_Num,DDK) .eq. Edge_Num) then 
             DM_Connect(3,Edge_Num,DDK)=-1
          endif
          
          if (Edge_Num .eq. 1) then 
             if (DM_Connect(2,Edge_Num,DDK) .eq. 4) then 
                DM_Connect(3,Edge_Num,DDK)=-1
             endif 
          endif 
        
          if (Edge_Num .eq. 2) then 
             if (DM_Connect(2,Edge_Num,DDK) .eq. 3) then 
                DM_Connect(3,Edge_Num,DDK)=-1
             endif 
          endif         
          
          if (Edge_Num .eq. 3) then 
             if (DM_Connect(2,Edge_Num,DDK) .eq. 2) then 
                DM_Connect(3,Edge_Num,DDK)=-1
             endif 
          endif
          
          if (Edge_Num .eq. 4) then 
             if (DM_Connect(2,Edge_Num,DDK) .eq. 1) then 
                DM_Connect(3,Edge_Num,DDK)=-1
             endif 
          endif     
          
!          write(*,*)'Patch',DM_Connect(3,Edge_Num,DDK)        
       enddo ! Edge_Num
          
    enddo ! DDK
    
    
    return 
    
  end subroutine Init_patch_direction2
  
 
  !==============================================================================

  subroutine alloc_mem_domain_grid_variables(mgl)
    !--Declare subroutine argument
    implicit none
    !   integer PolyDegreeN_max,TotNumDomain
    !
    !--Declare local argument
    integer ierr
    integer:: mgl
    !
    !--Start subroutine
    !allocate required memory for basic domain parameters
    allocate(PolyDegN_DM(2,TotNum_DM,mgl), stat=ierr)
    if ( ierr .ne. 0) then
       write(*,*)'Can not allocatable PolyDegreeN_Domain variables'
    endif
    PolyDegN_DM=0
    
    allocate(DM_Vertex(4,2,TotNum_DM), stat=ierr)
    if ( ierr .ne. 0) then
       write(*,*)'Can not allocatable DM_Vertex variables'
    endif
    DM_Vertex=0.d0
     
    allocate( x1(0:PolyDegN_Max(1),0:PolyDegN_Max(2),1:TotNum_DM,mgl), &
              x2(0:PolyDegN_Max(1),0:PolyDegN_Max(2),1:TotNum_DM,mgl), stat=ierr)

    if (ierr .ne. 0) then
       write(*,*)'Can not allocate variable x and y '
    endif
    x1=0.d0; x2=0.d0

    allocate(DM_Connect(1:3,1:4,1:TotNum_DM), stat=ierr)
    if (ierr .ne. 0) then 
       write(*,*)'Can not allocate veriable DM_Connect'
    endif
    DM_Connect=0
    
    allocate(DMShapeType(1:4,1:TotNum_DM), &
             DMShapePara(1:7,1:4,1:TotNum_DM), stat=ierr)
    DMShapeType=0
    DMShapePara=0.d0
    
    return
    
  end subroutine alloc_mem_domain_grid_variables


!==============================================================================
  
  subroutine alloc_mem_metric_variables(Deg1_Max,Deg2_Max,TotNumDomain,mgl)
    implicit none
    integer:: Deg1_Max, Deg2_Max, Deg3_Max, TotNumDomain, mgl

    integer:: ierr
    allocate( &
         DM_Metric(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
           Jacobin(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
              dxi1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
              dxi2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dxi1_dx1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dxi1_dx2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dxi2_dx1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dxi2_dx2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dx1_dxi1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dx1_dxi2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dx2_dxi1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
          dx2_dxi2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), &
            dtrans(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain,mgl), stat=ierr)

    if (ierr .ne. 0) then 
       write(*,*)'Cannot allocate memory for Metric Variables'
       write(*,*)'Abort!'
       stop
    endif

    DM_Metric=0.d0
    Jacobin=0.d0

    dxi1=0.d0; dxi2=0.d0;

    dxi1_dx1=0.d0; dxi1_dx2=0.d0;
    dxi2_dx1=0.d0; dxi2_dx2=0.d0;

    dx1_dxi1=0.d0; dx1_dxi2=0.d0;
    dx2_dxi1=0.d0; dx2_dxi2=0.d0;

    dtrans=0.d0

    return 
  end subroutine alloc_mem_metric_variables

!=================================================================
  subroutine alloc_mem_norvec(N_max,Num_Domain,mgl)
    implicit none
    integer:: N_max
    integer:: Num_Domain, mgl

    integer:: ierr

    ! subroutine begins
    allocate (NorVec_x1(0:N_max,4,Num_Domain,mgl), &
              NorVec_x2(0:N_max,4,Num_Domain,mgl), &
              NorVec_mg(0:N_max,4,Num_Domain,mgl), &
              JacNorVec(0:N_max,4,Num_Domain,mgl), &
              stat=ierr )

    if (ierr .ne. 0) then
       write(*,*)'Error message from MD2D_Grid.f90:'
       write(*,*)'Cannot allocate memory for Normal Vector Variables'
       write(*,*)'abort!'
    endif

    NorVec_x1=0.d0
    NorVec_x2=0.d0
    NorVec_mg=0.d0
    JacNorVec=0.d0
    
    return
  end subroutine alloc_mem_norvec
!=================================================================

end module MD2D_Grid
