      Program Poisson2D
      use constants    ! constants.f90
      use MD2D_Grid    ! MD3D_Grid.f90
      use Legendre     ! Legendre.f90
      use State_Var    ! State_Var.f90
      use ERR_Var      ! ERR_Var.f90
      use CG_Var       ! CG_Var.f90
      use Multigrid_Var! Multigrid_Var.f90
      use GMRES        ! GMRES.f90
      use INPUT
      implicit none
      
      integer :: i, j, sovflag, ierr
      integer :: iterNum, tot_size, LD1, LD2, n

      OPEN (UNIT=10,FILE='logfile')
      OPEN (UNIT=9,FILE='pspm2d.rea',STATUS='OLD',iostat=ierr)
      
      call setupinvar
      call rdparam(9)
      
      democase = PARAM(1)
      sovflag = PARAM(2)
      mg_lmax = PARAM(7)
      write(10,*)'mg_lmax',mg_lmax
      
!     Initialize stage
!-----Initialize universal constants
      call init_constants

!-----Initialize Domain Parameter 
      call Input_Domain_Grid_Parameters(mg_lmax)        !-- MD2D_Grid.f90
      write(10,*)'Complete Initializing Grid Parameters'

      call Geometric_Parameters_On_Screen(mg_lmax)     !-- MD2D_Grid.f90
      
!-----Initialize Legendre Pseudospectral Module
!     1. Gauss-Lobatto-Legendre Grid
!     2. Gauss-Lobatto-Legendre Quadurature Weights
!     3. Legenedre Differentiation Matrix
!
      call Init_LGL_Diff_Matrix_23D(PolyDegN_DM(1,1,mg_lmax), &
                                    PolyDegN_DM(2,1,mg_lmax), &
                                    0  )  !-- Legendre.f90
      write(10,*)'Complete Initializing Legendre &  
                 Pseudospectral Operators'
      
!-----Initialize computational nodes in physical space
!     by transfinite blending function
      
!     Initialize edge patch directions between adjacent domains
      call Init_Patch_Direction  !-- MD2D_Grid.f90
      write(10,*)'Complete Initializing Patching Direction'

!     Initialize physical grid points 
      call Init_Physical_Grid_2D(mg_lmax) !-- Grid2D_Pack.f90
      write(10,*)'Complete Initializing Grid points'
      
!     Initialize Metric variables
      call Init_Metric(mg_lmax) !-- Grid2D_Pack.f90
           
!     Initialize Normal Vector
      call Init_Normal_Vector(mg_lmax) !-- Grid2D_Pack.f90
           
!-----Initialize State Variables
      call Init_State_Variables(PolyDegN_DM(1,1,mg_lmax),TotNum_DM,mg_lmax)  !-- State_Var.f90
      write(10,*)'Complete Initializing State Variables'
      
!-----Initialize Material Parameters
      call Init_Material_Parameters(mg_lmax) !-- State_Pack.f90
      write(10,*)'Complete Initializing Material Variables'
      
!-----Initialize BC Parameters
      call Init_BC_Variables(mg_lmax) !-- State_Pack.f90
      write(10,*)'Complete Initializing Variables for &
                   Boundary Condition'

      call Initial_Field(mg_lmax) !--State_Pack.f90
      write(10,*)'Complete Initializing Initial Field'

      call FluxComp(mg_lmax) !--State_Pack.f90
      write(10,*)'Complete Flux computing'

!     Compute the iteration number
      iterNum = 0
      do DDK=1,TotNum_DM
         iterNum = iterNum + (PolyDegN_DM(1,DDK,mg_lmax)+1) &
                           * (PolyDegN_DM(2,DDK,mg_lmax)+1)
      enddo
      tot_size = (PolyDegN_DM(1,1,mg_lmax)+1) &
               * (PolyDegN_DM(2,1,mg_lmax)+1) &
               * TotNum_DM

      LD1 = PolyDegN_DM(1,1,mg_lmax); LD2 = PolyDegN_DM(2,1,mg_lmax)
      n = (PolyDegN_DM(1,1,mg_lmax) + 1)**2 * TotNum_DM



      if (democase .eq. 1) then
         
         call Init_CG(PolyDegN_DM(1,1,mg_lmax),TotNum_DM,iterNum)
!         write(10,*)'Complete Initializing Iterative Solver variables'

         if (sovflag .eq. 1) then

      
            call CG (potent,x_init,rhs,mg_lmax,iterNum,1e-16) !--CG_Pack.f90
            write(10,*)'Complete Operating Conjugate Gradient'

         else if (sovflag .eq. 2) then

            call alloc_gmres_var(tot_size,PolyDegN_DM(1,1,mg_lmax)&
            ,PolyDegN_DM(2,1,mg_lmax),TotNum_DM)
            write(10,*)'Complete Initializing GMRES variables'

            call HMH_GMRES(potent(0:ND1,0:ND2,1:TotNum_DM),x_init(0:ND1,0:ND2,1:TotNum_DM)&
            ,rhs(0:ND1,0:ND2,1:TotNum_DM),tot_size,1e-8,mg_lmax) !--HMH_GMRES.f90

         endif
      

!     Iterative solver with zero initial guess requires less iteration to converge then using Jacobi smoothing - Fixed
      else if (democase .eq. 2) then

         call Construct_ML_operator(LD1,LD2,n,mg_lmax) !--Precondition.f90
         write(10,*)'Complete Constructing ML Operator'

         call alloc_mem_Interpolatematrix_var(PolyDegN_DM(1,1,mg_lmax),TotNum_DM)
         call Interp_mat(mg_lmax)
         write(10,*)'Complete calling Interp_mat which set up &
         Jhx and Jhy'

         call Smoothing_Pack(LD1,LD2,mg_lmax)
         write(10,*)'Complete Smoothing Residue'

         call Init_CG(PolyDegN_DM(1,1,mg_lmax),TotNum_DM,iterNum)

         if (sovflag .eq. 1) then

            call CG (potent,x_vc,rhs,mg_lmax,iterNum,1e-16) !--CG_Pack.f90
            write(10,*)'Complete Operating Conjugate Gradient'

         else if (sovflag .eq. 2) then

            call alloc_gmres_var(tot_size,PolyDegN_DM(1,1,mg_lmax),PolyDegN_DM(2,1,mg_lmax),TotNum_DM)
            write(10,*)'Complete Initializing GMRES Variables'

            call HMH_GMRES(potent(0:ND1,0:ND2,1:TotNum_DM),x_vc(0:ND1,0:ND2,1:TotNum_DM),&
               rhs(0:ND1,0:ND2,1:TotNum_DM),tot_size,1e-8,mg_lmax) !--HMH_GMRES.f90

         endif      

!------------------------------------------------------------------
      else if (democase .eq. 3 ) then

         call alloc_mem_Lx_Ly_var(PolyDegN_DM(1,1,mg_lmax),TotNum_DM,mg_lmax)
         call Construct_Lx_Ly_operator(mg_lmax) !-Precondition.f90
         write(10,*)'Complete Constructing Lx Ly operator'

         call alloc_mem_Interpolatematrix_var(PolyDegN_DM(1,1,mg_lmax),TotNum_DM)
         call Interp_mat(mg_lmax)
         write(10,*)'Complete calling Interp_mat which set up & 
         Jhx and Jhy'

         call alloc_mem_ovlapsmooth_var(PolyDegN_DM(1,1,mg_lmax),TotNum_DM)
         call Smoothing_Pack_Overlapping(LD1,LD2,mg_lmax)
         write(10,*)'Complete Smoothing Residue'

         call CG (potent,x_vc,rhs,mg_lmax,iterNum,1e-16) !--CG_Pack.f90
         write(10,*)'Complete Operating Conjugate Gradient'
!----------------------------------------------------------------
      else if (democase .eq. 4) then

         Nk = param(5)
         write(10,*)'Total number of projection iteration:',Nk
      
         call alloc_mem_Lx_Ly_var(PolyDegN_DM(1,1,1),TotNum_DM,mg_lmax)
         call Construct_Lx_Ly_operator(1) !-Precondition.f90
         write(10,*)'Complete Constructing Lx Ly operator'

         call alloc_mem_Interpolatematrix_var(PolyDegN_DM(1,1,1),TotNum_DM)      
         call Interp_mat
         write(10,*)'Complete calling Interp_mat & 
         which set up Jhx and Jhy'

!         call hsmg_setup

         call alloc_mem_wrapper_var(PolyDegN_DM(1,1,1),PolyDegN_DM(1,1,2),TotNum_DM,Nk)
         call Projection_WRAPPER(1,Nk)
         write(10,*)'Complete calling Projection & 
         method with preconditioning'
      
         call copy(potent,x_precond,PolyDegN_DM(1,1,1),PolyDegN_DM(2,1,1),1)
!------------------------------------------------------------------
      else if (democase .eq. 5) then

         call alloc_mem_Lx_Ly_var(PolyDegN_DM(1,1,mg_lmax),TotNum_DM,mg_lmax)
         call Construct_Lx_Ly_operator(mg_lmax) !-Precondition.f90
!         write(10,*)'Complete Constructing Lx Ly operator'

         call alloc_mem_Interpolatematrix_var(PolyDegN_DM(1,1,mg_lmax),TotNum_DM)
         call Interp_mat(mg_lmax)
         write(10,*)'Complete calling Interp_mat &
            which set up Jhx and Jhy'

         call alloc_gmres_var(tot_size,PolyDegN_DM(1,1,mg_lmax),PolyDegN_DM(2,1,mg_lmax),TotNum_DM)
         write(10,*)'Complet initializing GMRES variables'

!         call hsmg_setup
         call alloc_mem_wrapper_var(PolyDegN_DM(1,1,mg_lmax),&
         PolyDegN_DM(1,1,mg_lmax-1),TotNum_DM,1)

         write(10,*)'Gmreswrapper',ND1,ND2
         call GMRES_WRAPPER(potent(0:ND1,0:ND2,1:TotNum_DM)&
                           ,rhs(0:ND1,0:ND2,1:TotNum_DM),tot_size,1e-8,mg_lmax)
!         call GMRES_CG_WRAPPER(potent(0:ND1,0:ND2,1:TotNum_DM)&
!                           ,rhs(0:ND1,0:ND2,1:TotNum_DM),tot_size,1e-8)
         write(10,*)'Complete calling GMRES with preconditioning'

      end if ! democase

      close(9)
      call Field_Final(1)

      call Init_ERR(PolyDegN_DM(1,1,mg_lmax),TotNum_DM) !--ERR_Var.f90
      write(10,*)'Complete Initializing Error Variables'

      call Error_Computing(mg_lmax)!--Error_pack.f90
      write(10,*)'Complete Error computing'      
      close(10)
      
      end Program Poisson2D
