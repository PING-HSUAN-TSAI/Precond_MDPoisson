      subroutine Projection_WRAPPER(LD1,LD2,l,step)

      use Legendre
      use MD2D_Grid
      use State_Var
      use CG_Var
      use Multigrid_Var

      implicit none
      
      integer :: step !-Gram-Schmidt step
      integer :: index_pro, index_k, Nx, Ny
      integer :: k, i, j, ND1p, ND2p,l, ND1c, ND2c
      integer :: LD1,LD2
      real(kind=8):: tmp_pro, glsc2, rnorm, err, glamax, tol
      real(kind=8) :: error_precond(0:LD1,0:LD2,1:TotNum_DM)

      Nx = PolyDegN_DM(1,1,l); Ny = PolyDegN_DM(2,1,l)
      
      tol = 1e-8
      
      call copy(r_wrapper,rhs,Nx,Ny,l)
      call copy(b_wrapper,r_wrapper,Nx,Ny,l)
      
      x_precond = 0.d0
      
      call axhelm2(Nx,Ny,x_precond,Ax_precond,l)
      call add3s2(r_wrapper,b_wrapper,Ax_precond,1.0,-1.0,Nx,Ny,l)
      
!     Projection start   
      do k = 1,step
      
      !-This have to be uncommnet so that the projection method 
      !-itself can work
      !   do DDK = 1, TotNum_DM
      !     ND1 = PolyDegN_DM(1,DDK,l)
      !     ND2 = PolyDegN_DM(2,DDK,l)
      !        do j=0,ND2
      !           do i=0,ND1
      !              p_cond(i,j,DDK) = r_wrapper(i,j,DDK) 
      !           enddo
      !        enddo
      !   enddo     
         call vcycle_projection(Nx,Ny,&
                                r_wrapper,p_cond,l)

!        Start Projection method ( fixme)
         if (k > 1) then
      
            do index_pro = 1, k-1
      
               call inner_product(Nx,Ny,p_cond,&
                                  (Proj_AP(:,:,:,index_pro)),pc_proAp(index_pro),l)
      
            enddo !index_pro
      
            do DDK = 1, TotNum_DM
               ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
               do j = 0,ND2
                  do i = 0, ND1
                     tmp_pro = 0
                     do index_k = 1, k-1
                        tmp_pro = tmp_pro &
                                + Proj_P(i,j,DDK,index_k) &
                                * pc_proAp(index_k)
                     enddo
                     Pap_pcond(i,j,DDK) = tmp_pro
                  enddo
               enddo
            enddo
      
            call add2s2(p_cond,Pap_pcond,-1.0,Nx,Ny,l)
         
         end if
      
      
         call axhelm2(Nx,Ny,p_cond,w,l)
      
         pAp_wrapper = glsc2(Nx,Ny,p_cond,w,l)
      
      
!        Normalize
         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
      
            Proj_P(0:ND1,0:ND2,DDK,k) = p_cond(0:ND1,0:ND2,DDK) & 
                                       / sqrt(pAp_wrapper)
         enddo
      
!        Normalize
         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
      
            Proj_AP(0:ND1,0:ND2,DDK,k) = w(0:ND1,0:ND2,DDK) & 
                                       / sqrt(pAp_wrapper)
         enddo

         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
            r_tmp(0:ND1,0:ND2,DDK) = r_wrapper(0:ND1,0:ND2,DDK) & 
                                   / (pAp_wrapper)
         enddo
      
         alpha_wrapper = glsc2(Nx,Ny,p_cond,r_tmp,l)
      
      
         call add2s2(x_precond,p_cond,alpha_wrapper,Nx,Ny,l)
         call add3s2(error_precond,v,x_precond,1.0,-1.0,Nx,Ny,l)
         err = glamax(error_precond,Nx,Ny,l)
      
         call add2s2(r_wrapper,w,-alpha_wrapper,Nx,Ny,l)
      
         rnorm = glsc2(Nx,Ny,r_wrapper,r_wrapper,l)
         rnorm = sqrt(rnorm)

         write(10,9998)k,pAp_wrapper,alpha_wrapper,rnorm,err

         if (abs(rnorm) .lt. tol) then
            exit
         endif
      

      enddo ! end of the projection method
      write(10,9997)k,err,tol
      
9998 format(' ',' ', 'Projection step',i5,4es24.15)
9997 format(' ',' ', 'End:',i5,'/',es10.4,1x,es10.4)
      end subroutine Projection_WRAPPER
!----------------------------------------------------------------
      subroutine vcycle_projection(LD1,LD2,ri,po,l)

      use Legendre
      use MD2D_Grid
      use State_Var
      use CG_Var
      use Multigrid_Var
      use Input

      implicit none
      
      integer :: i, j, LD1,LD2, vcycle
      integer :: m_vcycle, m_smooth, TotN, TotNc
      integer :: k, n, ND1p, ND2p,l, ND1c, ND2c, ind_JS
      integer :: Nx,Ny
      real(kind=8) :: smoothpar, err_sm, r_sm, err_vc, max_ef,glamax

!     vcycle variables
      real(kind=8) :: ri(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: r(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: po(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: x_vcp(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: errvec_vc(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: errvec_sm(0:LD1,0:LD2,1:TotNum_DM)

!     working arrays
      real(kind=8) :: tmp(0:LD1,0:LD2,1:TotNum_DM),tmp1(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp2(0:LD1,0:LD2,1:TotNum_DM),tmp3(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp4(0:LD1,0:LD2,1:TotNum_DM),tmp5(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp6(0:LD1,0:LD2,1:TotNum_DM)      
      real(kind=8) :: tmp7(0:LD1,0:LD2,1:TotNum_DM)      
      real(kind=8) :: tmp8(0:LD1,0:LD2,1:TotNum_DM)      


      Nx = PolyDegN_DM(1,1,l); Ny = PolyDegN_DM(2,1,l)

      call chk_amax('ri',ri,Nx,Ny,l) 
      call chk_amax('po',po,Nx,Ny,l)
      call copy(r,ri,Nx,Ny,l)

      m_vcycle = 1
      m_smooth = 1
      smoothpar = param(6)!2.d0/3.d0
      
      x_vcp = 0.d0
      
      do vcycle = 1, m_vcycle
      
!        Additive Jacobi Schwartz
         do k = 1,m_smooth

            z_ov = 0*x_vcp
      
!           Multiplying the inverse of the BlockL operator.
            do ind_JS = 1,1
               do DDK = 1, TotNum_DM
                  ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
      
                  tmp2(0:ND1,0:ND2,DDK) = &
                     Matmul(r(0:ND1,0:ND2,DDK), Sy_norm(0:ND2,0:ND2,DDK))
      
                  tmp3(0:ND1,0:ND2,DDK) = &
                     Matmul(Sx_t_norm(0:ND1,0:ND1,DDK), tmp2(0:ND1,0:ND2,DDK))
      
                  do j=0,ND2
                     do i=0,ND1
!                        n = (DDK-1)*(ND1+1)*(ND2+1) + i + j*(ND2+1) +1
                        n = i + j*(ND2+1) +1
                        tmp4(i,j,DDK) = Diagonal(n,DDK) * tmp3(i,j,DDK)
                     enddo
                  enddo
      
                  tmp5(0:ND1,0:ND2,DDK) = &
                     Matmul(tmp4(0:ND1,0:ND2,DDK), Sy_t_norm(0:ND2,0:ND2,DDK))
      
                  z_ov_sol(0:ND1,0:ND2,DDK) = &
                     Matmul(Sx_norm(0:ND1,0:ND1,DDK), tmp5(0:ND1,0:ND2,DDK))
      
                  z_ov(0:ND1,0:ND2,DDK) = z_ov(0:ND1,0:ND2,DDK) &
                                        + z_ov_sol(0:ND1,0:ND2,DDK)
      
               enddo! DDK
            enddo
               
!            do DDK=1,TotNum_DM
!               ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
!               do j=0,ND2
!                  do i=0,ND1
!                     n = i+1+j*(ND2+1)
!                     z_ov(i,j,DDK) = z_ov(i,j,DDK)/scale_c(n,DDK)
!                  enddo
!               enddo
!            enddo
!            call chk_amax('zas',z_ov,PolyDegN_Max(1),PolyDegN_Max(2),l)
      
!           update solution
            call add2s2(x_vcp,z_ov,smoothpar,Nx,Ny,l)

!           compute smoothing error
            call add3s2(errvec_sm,v,x_vcp,1.0,-1.0,Nx,Ny,l)
            err_sm = glamax(errvec_sm,Nx,Ny,l)
            call chk_amax('err',errvec_sm,Nx,Ny,l)

            call axhelm2(Nx,Ny,x_vcp,Ax_vc,l)
      
            call add2s2(r,Ax_vc,-1.0,Nx,Ny,l)
            r_sm = glamax(r,Nx,Ny,l)

            write(10,777)vcycle,m_smooth,err_sm,r_sm      

         enddo ! smooth

777  format(' ',' ', 'vcycle:' ,i2,1x,'smoothing:',1x,i2,2e10.3)

!     Coarse-grid restriction x <--- x + e, where e is approximated on coarse grid
      
         do DDK = 1 ,TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
            ND1c = PolyDegN_DM(1,DDK,l-1); ND2c = PolyDegN_DM(2,DDK,l-1)
      
            tmp(0:ND1,0:ND2c,DDK) = &
               Matmul(r(0:ND1,0:ND2,DDK), Ihy(0:ND2,0:ND2c,DDK))
      
            rc_smooth(0:ND1c,0:ND2c,DDK) = &
               Matmul(Ihx_transpose(0:ND1c,0:ND1,DDK),tmp(0:ND1,0:ND2c,DDK))
         enddo

!         call hsmg_tnsr(rc_smooth(0:mg_nx(2),0:mg_ny(2),1:TotNum_DM),mg_nh(2)&
!         ,r(0:mg_nx(3),0:mg_ny(3),1:TotNum_DM),mg_nh(3),mg_jht(:,:,2),mg_jh(:,:,2)&
!         ,PolyDegN_DM(1,1,1),PolyDegN_DM(2,1,1))

         TotNc = 0; TotN = 0
      
         do DDK=1,TotNum_DM
            TotNc = TotNc + (PolyDegN_DM(1,DDK,l-1)+1) &
                          * (PolyDegN_DM(2,DDK,l-1)+1)
            TotN = TotN + (PolyDegN_DM(1,DDK,l)+1) &
                        * (PolyDegN_DM(2,DDK,l)+1)
         enddo
      
         call chk_amax('xc1',xc_in,PolyDegN_DM(1,1,l-1),PolyDegN_DM(2,1,l-1),l-1)
         write(10,*)'check point: before CG'
         xc_in = 0
         call chk_amax('xci',xc_in,PolyDegN_DM(1,1,l-1),PolyDegN_DM(2,1,l-1),l-1)

!        Call CG to solve for ec
!         call CG(ec,xc_in,rc_smooth,l-1,6000,1e-20)
         call CG(ec(0:ND1c,0:ND2c,1:TotNum_DM),xc_in(0:ND1c,0:ND2c,1:TotNum_DM),&
         rc_smooth(0:ND1c,0:ND2c,1:TotNum_DM),l-1,TotNc,1e-16)
      
         do DDK = 1 ,TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
            ND1c = PolyDegN_DM(1,DDK,l-1); ND2c = PolyDegN_DM(2,DDK,l-1)
     
!     Interpolate ec to ef which is from coarse to fine
!--------------------------------------------------------------------------
            tmp1(0:ND1c,0:ND2,DDK) = &
               Matmul(ec(0:ND1c,0:ND2c,DDK),Ihy_transpose(0:ND2c,0:ND2,DDK))
      
            ef(0:ND1,0:ND2,DDK) = &
               Matmul(Ihx(0:ND1,0:ND1c,DDK),tmp1(0:ND1c,0:ND2,DDK))
!---------------------------------------------------------------------------
         enddo


!         call hsmg_intp(ef(0:mg_nx(3),0:mg_ny(3),1:TotNum_DM),ec(0:mg_nx(2),0:mg_ny(2),1:TotNum_DM),&
!         2,PolyDegN_DM(1,1,1),PolyDegN_DM(2,1,1))

         max_ef = glamax(ef,Nx,Ny,l)
         call add2s2(x_vcp,ef,1.0,Nx,Ny,l)
      
!        Compare the x_vc with exact solution after doing smoothing and coarse correction
         call add3s2(errvec_vc,v,x_vcp,1.0,-1.0,Nx,Ny,l)
         err_vc = glamax(errvec_vc,Nx,Ny,l)

         write(10,1003)vcycle,err_vc,max_ef

1003 format(' ',' ', 'vcycle:',i2,1x,'err_vc:',es10.3,1x,'ef:',es10.3)
      
      enddo ! vcycle
      
      call copy(po,x_vcp,Nx,Ny,l)
      end subroutine
