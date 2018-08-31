      subroutine hmh_gmres(phi,x,res,n,tol,mgl)
!-----------------------------------------------------------------------
      use GMRES
      use MD2D_GRID
!-----------------------------------------------------------------------
!(potent,rhs,h1,h2,mult,dmask(1,1),isd,imsh,npts,tol)
!..... phi: return value
!..... res: redisual    
!..... h1 & h2 from :  Helmholtz operator AU = h1*[A]u + h2*[B]u
!..... n: total size of res
!..... tol: tolerance         
      implicit none
      integer:: n,outer,isd,imsh
      real(kind=8):: phi(n),res(n),h1(n),h2(n),x(n)
      real(kind=8):: tol,alpha,l,temp, rnorm, tolps, tolpss
      integer:: i, j, k, iconv, iter, m,test, mgl
      
      iter  = 0
      m     = lgmres ! FIXEME: You can define the size of m (Krylov subspce dimension)
      
      tolps = tol
      tolpss= tolps
      
      iconv = 0
      ND1=PolyDegN_DM(1,1,mgl); ND2=PolyDegN_DM(2,1,mgl)

      write(*,*)'here'
!     Initial guess and upper-Hessenberg form of H
      do i=1,n
         gmres_x(i) = x(i)
      enddo
      
      do j=1,m
         do i =1,m
            gmres_h(i,j) = 0
         enddo
      enddo
      
!     GMRES outer iterations
      outer = 0
      
      do while (iconv.eq.0.and.iter.lt.10000)
         outer = outer+1
         if(iter.eq.0) then
            do i=1,n
               gmres_r(i) = res(i) ! r = res
            enddo                
          else
!     update residual
            do i =1,n
               gmres_r(i) = res(i) ! r = res
            enddo
            call axhelm2(ND1,ND2,gmres_x(1:n),&
                  gmres_w(1:n),mgl) ! w = Ax
!     Compute initial residual
!     r = r - w
            do i=1,n
               gmres_r(i) = gmres_r(i) - gmres_w(i)
            enddo
         endif

!     gamma  = (r,r)
!     Compute initial residual norm

         call inner_product(ND1,ND2,gmres_r(1:n),gmres_r(1:n),gamma(1),mgl)
         gamma(1) = sqrt(gamma(1))  ! gamma  = sqrt{ (r,r) }
                                    !      1
      
!     check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
      
!     Define first Krylov vector
!     v  = r / gamma
         do i =1,n
            gmres_v(i,1) = gmres_r(i) * temp
         enddo
      
!     GMRES inner iteration
         do j=1,m
            iter = iter+1
!     Matrix-vector product
            call axhelm2(ND1,ND2,gmres_v(1:n,j),gmres_w(1:n),mgl)
!------------------------------------------------------------------------
!     Modified Gram-Schmidt orthogonalization
            do i=1,j
!     h    = (w,v)
               call inner_product(ND1,ND2,&
                     gmres_w,gmres_v(1:n,i),gmres_h(i,j),mgl)
!     w = w - h    v
               do k=1,n
                  gmres_w(k) = gmres_w(k) &
                             - gmres_h(i,j) * gmres_v(k,i)
               enddo
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
      
            call inner_product(ND1,ND2,&
                   gmres_w(1:n),gmres_w(1:n),alpha,mgl)
            alpha = sqrt(alpha) ! alpha =  \/ (w,w)
      
            if(alpha.eq.0.) goto 900  !converged
      
            l = sqrt(gmres_h(j,j)*gmres_h(j,j)+alpha*alpha)
            temp = 1./l
            gmres_c(j) = gmres_h(j,j) * temp
            gmres_s(j) = alpha  * temp
            gmres_h(j,j) = l
            gamma(j+1) = -gmres_s(j) * gamma(j)
            gamma(j)   =  gmres_c(j) * gamma(j)
      
      
            rnorm = abs(gamma(j+1))!*norm_fac
      
            if (rnorm .lt. tolps) goto 900 !converged
            if (j.eq.m) goto 1000 !not converged, restart
      
            temp = 1./alpha
      
!     Define next Krylov vector
!     v = w / alpha
            do i =1,n
               gmres_v(i,j+1) = gmres_w(i) * temp
            enddo
      
         enddo !inner iteration
900   iconv = 1
1000  continue

!     back substitution
!             -1
!        c = H   gamma
      
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - gmres_h(k,i)*gmres_c(i)
            enddo
            gmres_c(k) = temp/gmres_h(k,k)
         enddo

!     sum up Arnoldi vectors form approximate solution
         do i=1,j
!     x = x + c  z
!             i  i
            do k=1,n
               gmres_x(k) = gmres_x(k) + gmres_c(i) * gmres_v(k,i)
            enddo
         enddo
      enddo !do while
      
9000 continue
      
      do i = 1,n
         phi(i) = gmres_x(i)
      enddo
      
      write(10,9999) iter,tolpss
      
9999 format(' ',' ',' gmres   : iteration#',i5,1pes12.4)
      return
      end subroutine
      
