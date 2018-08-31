module Legendre
implicit none
real(kind=8), save, allocatable:: LG_grids(:),LG_weights(:)
real(kind=8), save, allocatable:: ValuesOfPolyNatGrids(:)
real(kind=8), save, allocatable:: LG_DiffMat(:,:)

real(kind=8), save, allocatable:: D1_DiffMat(:,:,:) !(0:N_max,0:N_max,Degree_max)
real(kind=8), save, allocatable::  D1_Vertex(:,:,:) !(0:N_max,2,Degree_max)
real(kind=8), save, allocatable:: D2_DiffMat(:,:,:)
real(kind=8), save, allocatable::   LGLCoord(:,:)   !(0:N_max,Degree_Max)
real(kind=8), save, allocatable:: LGLWeights(:,:)   !(0:N_max,Degree_max)


! The following variables are for three dimensional problems
real(kind=8), save, allocatable:: Diff_xi1(:,:,:)
real(kind=8), save, allocatable:: Diff_xi2(:,:,:)
! real(kind=8), save, allocatable:: Diff_xi3(:,:,:)
! To save memory we use Diff_xi3 as same as Diff_xi1
! Diff_xi2 is the transpose of Diff_xi

contains

!==============================================================================
 subroutine Init_LGL_Diff_Matrix_23D(Degree1_Max,Degree2_Max,Degree3_Max)
 implicit none 
 integer:: Degree1_Max, Degree2_Max, Degree3_Max
 
 integer:: N_max, N13_max
 integer:: N
 integer:: ierr,i,j

 N_max=max(Degree1_Max,Degree2_Max,Degree3_Max)
 allocate(  LGLCoord(0:N_max,N_max), &
          LGLWeights(0:N_max,N_max), stat=ierr)
 if (ierr .ne. 0) then
    write(*,*)'Can not allocate memory for LGLCoord and LGLWeights'
    stop
 endif
 LGLCoord=0.d0
 LGLWeights=0.d0

 ! allocate memory for differential matrix in xi1 direction
 N_max=max(Degree1_Max,Degree2_Max,Degree3_Max)
 allocate(Diff_xi1(0:N_Max,0:N_Max,N_Max), &
          Diff_xi2(0:N_Max,0:N_Max,N_Max), stat=ierr)

 if (ierr .ne. 0) then
    write(*,*)'Legendre.f90:' 
    write(*,*)'Fail to allocate memory for Diff_xi1 and Diff_xi2'
    write(*,*)'Abort!'
    stop 
 endif

 Diff_xi1=0.d0
 Diff_xi2=0.d0

 !allocate variables for LGL_grid and initialize it

 call Alloc_Mem_LG_grid(N_max) 

 call Alloc_Mem_Leg_Diff_Mat(N_max)
 
 ! set up LGL grid and differential matrix computational space
    
 do N=1,N_max
       
    call init_LGL_grids(N)
       
      LGLCoord(0:N,N) =   LG_grids(0:N)
      LGLWeights(0:N,N) = LG_weights(0:N)

!      do j=0,N
!         write(10,*)N,j,LGLWeights(j,N)
!      enddo

    call init_diff_matrix(N,N_max)  ! out put LG_DiffMat
 
      Diff_xi1(0:N,0:N,N)=LG_DiffMat(0:N,0:N)

      Diff_xi2(0:N,0:N,N)=transpose(LG_DiffMat(0:N,0:N))
      
 enddo           

  
    return

  end subroutine Init_LGL_Diff_Matrix_23D
  !==============================================================================
  subroutine Init_LGL_Diff_Matrix_1D(Degree_Max)
    implicit none
    integer:: Degree_Max
    
    integer:: N       ! Degree of the Polynomial
    integer:: ierr
    
    ! Differential Matrix of all Degree
    allocate(D1_DiffMat(0:Degree_Max,0:Degree_Max,Degree_Max),stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for D1_DiffMat'
       stop
    endif
    D1_DiffMat=0.d0
    
    allocate(D2_DiffMat(0:Degree_Max,0:Degree_Max,Degree_Max),stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for D2_DiffMat'
       stop
    endif
    D2_DiffMat=0.d0
    
    allocate(D1_Vertex(0:Degree_Max,1:2,Degree_Max),stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for D1_DiffMat_Vertex'
       stop
    endif
    D1_Vertex=0.d0

    allocate(  LGLCoord(0:Degree_Max,Degree_Max), &
             LGLWeights(0:Degree_Max,Degree_Max), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for LGLCoord and LGLWeights'
       stop
    endif
    LGLCoord=0.d0; write(*,*)'here'
    LGLWeights=0.d0
    
    !allocate variables for LGL_grid and initialize it

    call Alloc_Mem_LG_grid(Degree_Max)

    call Alloc_Mem_Leg_Diff_Mat(Degree_Max)
    
    ! set up LGL grid and differential matrix computational space
    
    do N=1,Degree_Max
       
       call init_LGL_grids(N)
       
         LGLCoord(0:N,N) =   LG_grids(0:N);

       LGLWeights(0:N,N) = LG_weights(0:N)
       
       call init_diff_matrix(N,Degree_Max)  ! out put LG_DiffMat
       
       ! set up D1 matrix and  boundary operator
      
       D1_DiffMat( 0:N , 0:N , N )=LG_DiffMat(0:N,0:N)
        D1_Vertex( 0:N ,  1  , N )=LG_DiffMat(0,0:N)
        D1_Vertex( 0:N ,  2  , N )=LG_DiffMat(N,0:N)
      
            
      ! set up D2 matrix

       D2_DiffMat( 0:N , 0:N , N ) = Matmul(LG_DiffMat(0:N,0:N),&
                                            LG_DiffMat(0:N,0:N))

      
   
    enddo

    return
    
  end subroutine Init_LGL_Diff_Matrix_1D

  !==============================================================================
  subroutine init_diff_matrix(N,NM)
    implicit none
    integer N,NM !NM is the Maximum Degree
    
    !     subroutine begins
    call DMLEGL(N,NM,LG_grids(0),ValuesOfPolyNatGrids(0),LG_DiffMat(0,0))
    
    return
    
  end subroutine init_diff_matrix

  !==============================================================================
  subroutine init_LGL_grids(N)
    implicit none
    integer N
    
    !     subroutine begins
    call ZELEGL(N, LG_grids, ValuesOfPolyNatGrids)
    call WELEGL(N, LG_grids, ValuesOfPolyNatGrids, LG_weights)
    
    return
    
  end subroutine init_LGL_grids

  !==============================================================================
  subroutine Alloc_Mem_LG_grid(N)
    implicit none
    integer N ! degree of polynomial N and total number of grids is N+1
    
    integer ierr
    
    !     subroutine begin
    allocate(LG_grids(0:N), LG_weights(0:N), ValuesOfPolyNatGrids(0:N), &
         stat=ierr)
    if ( ierr .ne. 0 ) then
       write(*,*)'Can not allocate Legendre grids'
       write(*,*)'Abort! '
       stop
    endif

    return
    
  end subroutine Alloc_Mem_LG_grid
  !==============================================================================
  subroutine Alloc_Mem_Leg_Diff_Mat(N)
    implicit none
    integer N
    
    integer ierr
    
    !     subroutine begins
    allocate(LG_DiffMat(0:N,0:N), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)' Can not allocatable Memory of variables LG_DiffMat '
       write(*,*)' Abort! '
       stop
    endif
    
  end subroutine Alloc_Mem_Leg_Diff_Mat

end module Legendre
