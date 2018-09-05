      subroutine GMRES_WRAPPER(phi,res,n,tol,l)
      use Legendre
      use MD2D_Grid
      use GMRES
      implicit none

      integer:: n, outer, Nx, Ny
      integer:: i, j, k, s, iconv, iter, l, m
      real(kind=8):: phi(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8):: res(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8):: tmp(0:PolyDegN_DM(1,1,l),0:PolyDegN_DM(2,1,l),1:TotNum_DM)
      real(kind=8):: tol,alpha,temp, rnorm, tolps, tolpss, ll
      real(kind=8):: glsc2

      tmp = 0.d0

      iter  = 0
      m     = lgmres ! FIXEME: You can define the size of m (Krylov subspce dimension)

      tolps = tol
      tolpss= tolps

      iconv = 0

      Nx=PolyDegN_DM(1,1,l); Ny=PolyDegN_DM(2,1,l)

      call chk_amax('res',res,Nx,Ny,l)

!     Initial guess and upper-Hessenberg form of H
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l);ND2=PolyDegN_DM(2,DDK,l)
         gmres_xx(0:ND1,0:ND2,DDK) = 0
      enddo 

      do j=1,m
         do i =1,m
            gmres_h(i,j) = 0
         enddo
      enddo

!     GMRES outer iterations
      outer = 0

      do while (iconv.eq.0.and.iter.lt.100)
         outer = outer+1
         if(iter.eq.0) then
            call copy  (gmres_rr,res,Nx,Ny,l) !-r = res
            call chk_amax('it0',gmres_rr,Nx,Ny,l)
         else
!     update residual
            call copy  (gmres_rr,res,Nx,Ny,l) !-r = res
            call chk_amax('it1',gmres_rr,Nx,Ny,l)

            call axhelm2(Nx,Ny,gmres_xx,gmres_ww,l) !-w = Ax
            call chk_amax('it2',gmres_ww,Nx,Ny,l)

!-Compute initial residual r = r - w
            call add2s2(gmres_rr,gmres_ww,-1.0,Nx,Ny,l)

         endif

!     Compute initial residual norm
         gamma(1) = sqrt(glsc2(Nx,Ny,gmres_rr,gmres_rr,l)) ! gamma  = sqrt{ (r,r) }
         write(10,*)'gamma',gamma(1)

!     Check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)

         call chk_amax('it3',gmres_rr,Nx,Ny,l)

!     Define first Krylov vector
!     v  = r / gamma
         do DDK = 1, TotNum_DM
            ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
            gmres_vv(0:ND1,0:ND2,DDK,1) = gmres_rr(0:ND1,0:ND2,DDK) &
                                   / (gamma(1))
         enddo

         call chk_amax('it4',gmres_rr,Nx,Ny,l)
         call chk_amax('it5',gmres_vv(:,:,:,1),Nx,Ny,l)

!     GMRES inner iteration
         do j=1,m
            iter = iter+1
            call copy  (gmres_ww,gmres_vv(:,:,:,j),Nx,Ny,l)
            write(*,*)'j',j

            do DDK=1,TotNum_DM
               ND1=PolyDegN_DM(1,DDK,l);ND2=PolyDegN_DM(2,DDK,l)
               do k=0,ND2
                  do s=0,ND1
                     tmp(s,k,DDK) = 0.d0 !gmres_zz(s,k,DDK,j)
                  enddo
               enddo
            enddo

            call chk_amax('it6',gmres_zz(:,:,:,j),Nx,Ny,l)

            call vcycle_projection(Nx,Ny,gmres_ww,tmp,l)
            call copy  (gmres_zz(:,:,:,j),tmp,Nx,Ny,l)

            call chk_amax('it7',gmres_zz(:,:,:,j),Nx,Ny,l)
            call chk_amax('it8',gmres_ww,Nx,Ny,l)
!     Matrix-vector product
            call axhelm2(Nx,Ny,gmres_zz(:,:,:,j),gmres_ww,l)
            call chk_amax('it9',gmres_ww,Nx,Ny,l)
!------------------------------------------------------------------------
!     Modified Gram-Schmidt orthogonalization
            do i=1,j
!     h    = (w,v)
               gmres_h(i,j) = glsc2(Nx,Ny,gmres_ww,gmres_vv(:,:,:,i),l)
!     w = w - h    v
               call add2s2(gmres_ww,gmres_vv(:,:,:,i),-gmres_h(i,j),Nx,Ny,l)
            enddo
!-------------------------------------------------------------------------
!     Apply Givens rotations to new column of Hessenberg H
            do i=1,j-1
               temp = gmres_h(i,j)
               gmres_h(i  ,j) =  gmres_c(i)*temp &
                              +  gmres_s(i)*gmres_h(i+1,j)
               gmres_h(i+1,j) = -gmres_s(i)*temp &
                              +  gmres_c(i)*gmres_h(i+1,j)
            enddo
            
            alpha = glsc2(Nx,Ny,gmres_ww,gmres_ww,l)
            alpha = sqrt(alpha)

            if(alpha.eq.0.) goto 900  !converged

            ll = sqrt(gmres_h(j,j)*gmres_h(j,j)+alpha*alpha)
            temp = 1./ll
            gmres_c(j) = gmres_h(j,j) * temp
            gmres_s(j) = alpha  * temp
            gmres_h(j,j) = ll
            gamma(j+1) = -gmres_s(j) * gamma(j)
            gamma(j)   =  gmres_c(j) * gamma(j)


            rnorm = abs(gamma(j+1))!*norm_fac
            write (10,66) iter,tolpss,rnorm
66    format(i5,1p2e12.5,' gmres_mg')

            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart

            temp = 1./alpha

!     Define next Krylov vector
!     v = w / alpha
            do DDK = 1, TotNum_DM
               ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
               gmres_vv(0:ND1,0:ND2,DDK,j+1) = gmres_ww(0:ND1,0:ND2,DDK) * temp
            enddo

         enddo !inner iteration
900      iconv = 1
1000  continue

!        back substitution
!             -1
!        c = H   gamma
!        write(6,*) 'start solving least squre problem'

         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - gmres_h(k,i)*gmres_c(i)
            enddo
               gmres_c(k) = temp/gmres_h(k,k)
         enddo

!     Sum up Arnoldi vectors form approximate solution
         do i=1,j
            call add2s2(gmres_xx,gmres_zz(:,:,:,i),gmres_c(i),Nx,Ny,l)
            ! x = x + c  z
            !          i  i
         enddo
      enddo !do while

9000 continue

      call copy  (phi,gmres_xx,Nx,Ny,l)

      write(10,999) iter,tolpss
 999 format(' ','gmres_mg: iteration#',i5,1es12.4)

      end
