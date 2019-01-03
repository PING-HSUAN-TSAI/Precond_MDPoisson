      subroutine Construct_ML_operator(LD1,LD2,lpts,l)
      use constants
      use Legendre
      use MD2D_Grid
      use State_Var
      use Multigrid_Var

      implicit none
      integer:: ND1p, ND2p, n, lpts, LD1, LD2
      integer:: i, j, l
      
      real(kind=8) :: dudx(0:LD1,0:LD2,1:lpts), dudy(0:LD1,0:LD2,1:lpts)
      real(kind=8) :: temp1(0:LD1,0:LD2,1:lpts,1:TotNum_DM), temp2(0:LD1,0:LD2,1:lpts,1:TotNum_DM)
      real(kind=8) :: BvEdge(0:LD1,1:4,1:lpts,1:TotNum_DM)
      real(kind=8) :: BITF1(0:LD1,1:4,1:lpts,1:TotNum_DM), BITF2(0:LD1,1:4,1:lpts,1:TotNum_DM)
      
      
      call alloc_mem_MLoperator_var(PolyDegN_DM(1,1,l),TotNum_DM)
!      call alloc_mem_Multigrid_var(PolyDegN_DM(1,1,1),TotNum_DM,1)
      
      call Initial_I_U(PolyDegN_DM(1:2,1:TotNum_DM,l),TotNum_DM)
      
        
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do n = 1, (ND1p * ND2p)
!     Compute
            dudx(0:ND1,0:ND2,n) = Matmul(Diff_xi1(0:ND1,0:ND1,ND1), &
                                              I_U(0:ND1,0:ND2,n,DDK))
      
            dudy(0:ND1,0:ND2,n) = Matmul( I_U(0:ND1,0:ND2,n,DDK),&
                                     Diff_xi2(0:ND2,0:ND2,ND2))
      
      
            temp1(0:ND1,0:ND2,n,DDK) = &
               dxi1_dx1(0:ND1,0:ND2,DDK,l) * dudx(0:ND1,0:ND2,n) &
             + dxi2_dx1(0:ND1,0:ND2,DDK,l) * dudy(0:ND1,0:ND2,n)
                                  
            temp2(0:ND1,0:ND2,n,DDK) = &
               dxi1_dx2(0:ND1,0:ND2,DDK,l) * dudx(0:ND1,0:ND2,n) &
             + dxi2_dx2(0:ND1,0:ND2,DDK,l) * dudy(0:ND1,0:ND2,n)
         enddo
      enddo
      
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do n = 1, (ND1p*ND2p)
!     Compute divergence b*F  (div dot bF)
            dudx(0:ND1,0:ND2,n) = ( dxi1_dx1(0:ND1,0:ND2,DDK,l) & 
                                *      temp1(0:ND1,0:ND2,n,DDK)  &
                                +   dxi1_dx2(0:ND1,0:ND2,DDK,l) & 
                                *      temp2(0:ND1,0:ND2,n,DDK) ) &
                                *    Jacobin(0:ND1,0:ND2,DDK,l) &
                                *      b_pts(0:ND1,0:ND2,DDK,l)
      
            dudy(0:ND1,0:ND2,n) = ( dxi2_dx1(0:ND1,0:ND2,DDK,l) &
                                *      temp1(0:ND1,0:ND2,n,DDK)  &
                                +   dxi2_dx2(0:ND1,0:ND2,DDK,l) & 
                                *      temp2(0:ND1,0:ND2,n,DDK) ) &
                                *    Jacobin(0:ND1,0:ND2,DDK,l) &
                                *      b_pts(0:ND1,0:ND2,DDK,l)
      
      
            ML(0:ND1,0:ND2,n,DDK) = - Matmul(Diff_xi1(0:ND1,0:ND1,ND1), &
                                                 dudx(0:ND1,0:ND2,n)) &
                                    - Matmul(    dudy(0:ND1,0:ND2,n), &
                                             Diff_xi2(0:ND2,0:ND2,ND2))
         enddo
      enddo
       
!     Check global A matrix
      !  open(1237,file='ML.text')
      !     do DDK = 1, TotNum_DM
      !        ND1=PolyDegN_DM(1,DDK,l)
      !        ND2=PolyDegN_DM(2,DDK,l)
      !        ND1p = ND1 + 1; ND2p = ND2 + 1
      !          do n = 1, ND1p*ND2p
      !            write(1237,*)n
      !              do j = 0, ND2 
      !                do i = 0, ND1 
      !                  write(1237,*)i,j,ML(i,j,n,DDK)
      !                enddo
      !              enddo
      !           enddo
      !       enddo
      !
      !       close(1237)
      call compute_PBC_MG(LD1,LD2,lpts,temp1,temp2,BvEdge,BITF1,BITF2,l)
      
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do n = 1, (ND1p*ND2p)

!     Start adding penalty
!     Side 4
            Edge_Num=4
            select case (BC_Type(Edge_Num,DDK))
               case(1)
               do j = 0,ND2
                  ML(0:ND1,j,n,DDK) = ML(0:ND1,j,n,DDK) &
                                    + (    tauD(0:ND1,j,Edge_Num,DDK,l) &
                                    - tau_tilde(0:ND1,j,Edge_Num,DDK,l)) &
                                    *  BvEdge(j,Edge_Num,n,DDK)
       
               enddo
      !do j=0,ND2
      !do i=0,ND1
      !write(185,fmt="(2I3,1x,E11.3e3)") i,j,b_cg(i,j,DDK)
      !enddo
      !enddo
               case(2,3)
      
               do j = 0,ND2
                  ML(0:ND1,j,n,DDK) = ML(0:ND1,j,n,DDK) &
                                    +    tauND1(0:ND1,j,Edge_Num,DDK,l) &
                                    *  BvEdge(j,Edge_Num,n,DDK)
       
               enddo
      
      ! write(*,*)'PBC4',PBC(0:ND2,Edge_Num,DDK)
      !write(*,*)'BvEdge',BvEdge(0:ND2,Edge_Num,DDK)
      !write(*,*)'tau4',tauND(0:ND2,Edge_Num,DDK)
      
      
               case(0)
!     Add interface penalty : first set of dirichlet
               do j=0,ND2
                  ML(0:ND1,j,n,DDK) = ML(0:ND1,j,n,DDK) &
                                    + tau1(0:ND1,j,Edge_Num,DDK,l)  &
                                    * BITF1(j,Edge_Num,n,DDK)
                  ML(0,0:ND2,n,DDK) = ML(0,0:ND2,n,DDK) &
                                    + tau3(j,0:ND2,Edge_Num,DDK,l) &
                                    * BITF1(j,Edge_Num,n,DDK)     
               enddo
               ML(0,0:ND2,n,DDK) = ML(0,0:ND2,n,DDK) &
                                 + tau2(0:ND2,Edge_Num,DDK,l) &
                                 * BITF1(0:ND2,Edge_Num,n,DDK)
!     Add interface penalty : second set of dirichlet
               ML(0,0:ND2,n,DDK) = ML(0,0:ND2,n,DDK) &
                                 + Sigma_tild(0:ND2,Edge_Num,DDK,l) &
                                 *    BITF1(0:ND2,Edge_Num,n,DDK)
!     Add interface penalty : neumann
               ML(0,0:ND2,n,DDK) = ML(0,0:ND2,n,DDK) &
                                 + SigmaHat(0:ND2,Edge_Num,DDK,l) &
                                 *  BITF2(0:ND2,Edge_Num,n,DDK)
      
            end select
!     Side 2
            Edge_Num=2
            select case (BC_Type(Edge_Num,DDK))
               case(1)
               do j=0,ND2
                  ML(0:ND1,j,n,DDK)  = ML(0:ND1,j,n,DDK) &
                                     + (    tauD(0:ND1,j,Edge_Num,DDK,l) &
                                     - tau_tilde(0:ND1,j,Edge_Num,DDK,l)) &
                                     *  BvEdge(j,Edge_Num,n,DDK)
               enddo
               case(2,3)
               do j=0,ND2
                  ML(0:ND1,j,n,DDK) = ML(0:ND1,j,n,DDK) &
                                    +    tauND1(0:ND1,j,Edge_Num,DDK,l) &
                                    *  BvEdge(j,Edge_Num,n,DDK)
               enddo
      !write(*,*)'PBC2',PBC(0:ND2,Edge_Num,DDK)
      
      !write(*,*)'tau2',tauND(0:ND2,Edge_Num,DDK)
      
      
               case(0)
      
!     Add interface penalty : first set of dirichlet
               do j=0,ND2
                  ML(0:ND1,j,n,DDK) = ML(0:ND1,j,n,DDK) &
                                    + tau1(0:ND1,j,Edge_Num,DDK,l) &
                                    * BITF1(j,Edge_Num,n,DDK)
                  ML(ND1,0:ND2,n,DDK) = ML(ND1,0:ND2,n,DDK) &
                                      + tau3(j,0:ND2,Edge_Num,DDK,l) &
                                      * BITF1(j,Edge_Num,n,DDK)
              
               enddo
               ML(ND1,0:ND2,n,DDK) = ML(ND1,0:ND2,n,DDK) &
                                   + tau2(0:ND2,Edge_Num,DDK,l) &
                                   * BITF1(0:ND2,Edge_Num,n,DDK)
!     Add interface penalty : second set of dirichlet
               ML(ND1,0:ND2,n,DDK) = ML(ND1,0:ND2,n,DDK) &
                                   + Sigma_tild(0:ND2,Edge_Num,DDK,l) &
                                   * BITF1(0:ND2,Edge_Num,n,DDK)
!     Add interface penalty : neumann
               ML(ND1,0:ND2,n,DDK) = ML(ND1,0:ND2,n,DDK) &
                                   + SigmaHat(0:ND2,Edge_Num,DDK,l) &
                                   * BITF2(0:ND2,Edge_Num,n,DDK)
      
            end select
      
!     Side 1
            Edge_Num=1
            select case (BC_Type(Edge_Num,DDK))
               case(1)
               do i=0,ND1
                  ML(i,0:ND2,n,DDK)  = ML(i,0:ND2,n,DDK) & 
                                     + (    tauD(i,0:ND2,Edge_Num,DDK,l) &
                                     - tau_tilde(i,0:ND2,Edge_Num,DDK,l)) &
                                     *  BvEdge(i,Edge_Num,n,DDK)
      
               enddo
               case(2,3)
               do i=0,ND1
                  ML(i,0:ND2,n,DDK) = ML(i,0:ND2,n,DDK) &
                                    +    tauND1(i,0:ND2,Edge_Num,DDK,l) &
                                    *  BvEdge(i,Edge_Num,n,DDK)
      
               enddo
               case(0)
      
!     Add interface penalty : first set of dirichlet
               do i=0,ND1
                  ML(i,0:ND2,n,DDK) = ML(i,0:ND2,n,DDK) &
                                    +     tau1(i,0:ND2,Edge_Num,DDK,l) &
                                    *  BITF1(i,Edge_Num,n,DDK)
                  ML(0:ND1,0,n,DDK) = ML(0:ND1,0,n,DDK) &
                                    +     tau3(i,0:ND1,Edge_Num,DDK,l) &
                                    *  BITF1(i,Edge_Num,n,DDK)        
               enddo
               ML(0:ND1,0,n,DDK) = ML(0:ND1,0,n,DDK) &
                                 +     tau2(0:ND1,Edge_Num,DDK,l) &
                                 *  BITF1(0:ND1,Edge_Num,n,DDK)
!     Add interface penalty : second set of dirichlet
               ML(0:ND1,0,n,DDK) = ML(0:ND1,0,n,DDK) &
                                 + Sigma_tild(0:ND1,Edge_Num,DDK,l) &
                                 *    BITF1(0:ND1,Edge_Num,n,DDK)
!     Add interface penalty : neumann
               ML(0:ND1,0,n,DDK) = ML(0:ND1,0,n,DDK) &
                                 + SigmaHat(0:ND1,Edge_Num,DDK,l) &
                                 *  BITF2(0:ND1,Edge_Num,n,DDK)
            end select
      
      !open(unit=19,file="tau.text")
      !Edge_Num=1
      !do j=0,ND2
      !do i=0,ND1
      !write(unit=19,fmt="(2I3,1x,3E15.8e3)") i,j,tauD(i,j,Edge_Num,DDK),tau_tilde(i,j,Edge_Num,DDK),&
      !(tauD(i,j,Edge_Num,DDK) - tau_tilde(i,j,Edge_Num,DDK))
      !
      !enddo
      !enddo
      !close(unit=19)
      !open(unit=55,file="b_cgv1.text")
      !do j=0,ND2
      !do i=0,ND1
      !
      !write(unit=55,fmt="(2I3,1X,E15.8e3)") i,j,b_cg(i,j,DDK)
      !enddo
      !enddo
      !close(unit=55)
      !
!     Side 3
            Edge_Num=3
            select case (BC_Type(Edge_Num,DDK))
               case(1)
               do i=0,ND1
                  ML(i,0:ND2,n,DDK) = ML(i,0:ND2,n,DDK) &
                                    + (    tauD(i,0:ND2,Edge_Num,DDK,l) &
                                    - tau_tilde(i,0:ND2,Edge_Num,DDK,l)) &
                                    *   BvEdge(i,Edge_Num,n,DDK)
               enddo
      
               case(2,3)
               do i=0,ND1
                  ML(i,0:ND2,n,DDK) = ML(i,0:ND2,n,DDK) &
                                    +    tauND1(i,0:ND2,Edge_Num,DDK,l) &
                                    *  BvEdge(i,Edge_Num,n,DDK)
               enddo
      !write(*,*)'PBC3',PBC(0:ND1,Edge_Num,DDK)
      
      ! write(*,*)'tau3',tauND(0:ND1,Edge_Num,DDK)
      
      !open(unit=56,file="b_cgv3.text")
      !    do j=0,ND2
      !        do i=0,ND1
      !
      !        write(unit=56,fmt="(2I3,1X,E15.8e3)") i,j,b_cg(i,j,DDK)
      !        enddo
      !    enddo
      !close(unit=56)
      
               case(0)
      
!     Add interface penalty : first set of dirichlet
               do i=0,ND1
                  ML(i,0:ND2,n,DDK) = ML(i,0:ND2,n,DDK) &
                                    +     tau1(i,0:ND2,Edge_Num,DDK,l) &
                                    *  BITF1(i,Edge_Num,n,DDK)
                  ML(0:ND1,ND2,n,DDK) = ML(0:ND1,ND2,n,DDK) &
                                      +    tau3(i,0:ND1,Edge_Num,DDK,l) &
                                      * BITF1(i,Edge_Num,n,DDK)          
               enddo
               ML(0:ND1,ND2,n,DDK) = ML(0:ND1,ND2,n,DDK) &
                                   +     tau2(0:ND1,Edge_Num,DDK,l) &
                                   *  BITF1(0:ND1,Edge_Num,n,DDK)
!     Add interface penalty : second set of dirichlet
               ML(0:ND1,ND2,n,DDK) = ML(0:ND1,ND2,n,DDK) &
                                   + Sigma_tild(0:ND1,Edge_Num,DDK,l) &
                                   *    BITF1(0:ND1,Edge_Num,n,DDK)
!     Add interface penalty : neumann
               ML(0:ND1,ND2,n,DDK) = ML(0:ND1,ND2,n,DDK) &
                                   + SigmaHat(0:ND1,Edge_Num,DDK,l) &
                                   *  BITF2(0:ND1,Edge_Num,n,DDK)
            end select
         enddo !n
      enddo !DDK
      
!     Multiply mass matrix in x,y direction
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do n = 1, (ND1p*ND2p)
            do j = 0, ND2
               do i = 0, ND1
                  ML(i,j,n,DDK) = ML(i,j,n,DDK) * LGLWeights(j,ND2) & 
                                * LGLWeights(i,ND1)
               enddo
            enddo
         enddo
      enddo
      
!     Check gloabl A matrix times Mass matrix
      open(1238,file='ML.text')
      do DDK = 1, TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do n = 1, ND1p*ND2p
            write(1238,*)n
            do j = 0, ND2
               do i = 0, ND1
                  write(1238,*)i,j,ML(i,j,n,DDK)
               enddo
            enddo
         enddo
      enddo
      
      close(1238)
      
      end subroutine Construct_ML_operator
      
!===================================================================
      subroutine compute_PBC_MG(LD1,LD2,lpts,temp1,temp2,BvEdge,BITF1,BITF2,l)
      use MD2D_Grid
      use Multigrid_Var
      use State_Var

      implicit none
      integer:: i, j, NDFix
      integer:: ND1p, ND2p, lpts, l, LD1,LD2, n
      real(kind=8) :: temp1(0:LD1,0:LD2,1:lpts,1:TotNum_DM), temp2(0:LD1,0:LD2,1:lpts,1:TotNum_DM)
      real(kind=8) :: vEdge(0:LD1,1:4,1:lpts,1:TotNum_DM), dvdnEdge(0:LD1,1:4,1:lpts,1:TotNum_DM)
      real(kind=8) :: BvEdge(0:LD1,1:4,1:lpts,1:TotNum_DM)
      real(kind=8) :: BITF1(0:LD1,1:4,1:lpts,1:TotNum_DM), BITF2(0:LD1,1:4,1:lpts,1:TotNum_DM)
      
!     assign v and dvdn on the edge number
      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,l);ND2 = PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
      
         do n = 1, (ND1p*ND2p)
            do Edge_Num=1,4
               select case (Edge_Num)
                  case (2,4)
      
                     ND=ND2; NDFix=0; if (Edge_Num .eq. 2) NDFix=ND1
      
                     vEdge(0:ND,Edge_Num,n,DDK) = I_U(NDFix,0:ND,n,DDK)
                     dvdnEdge(0:ND,Edge_Num,n,DDK) &
                        = Norvec_x1(0:ND,Edge_Num,DDK,l) * temp1(NDFix,0:ND,n,DDK) &
                        + Norvec_x2(0:ND,Edge_Num,DDK,l) * temp2(NDFix,0:ND,n,DDK)
      
                     BvEdge(0:ND,Edge_Num,n,DDK) &
                        =  bnd_alpha(0:ND,Edge_Num,DDK,l) * vEdge(0:ND,Edge_Num,n,DDK) &
                        +   bnd_beta(0:ND,Edge_Num,DDK,l) * bEdge(0:ND,Edge_Num,DDK,l) &
                        *   dvdnEdge(0:ND,Edge_Num,n,DDK)
      
      
                  case(1,3)
      
                     ND=ND1; NDFix=0; if (Edge_Num .eq. 3) NDFix=ND2
      
                     vEdge(0:ND,Edge_Num,n,DDK) = I_U(0:ND,NDFix,n,DDK)
                     dvdnEdge(0:ND,Edge_Num,n,DDK) &
                        = Norvec_x1(0:ND,Edge_Num,DDK,l) * temp1(0:ND,NDFix,n,DDK) &
                        + Norvec_x2(0:ND,Edge_Num,DDK,l) * temp2(0:ND,NDFix,n,DDK)
      
                     BvEdge(0:ND,Edge_Num,n,DDK) &
                        =  bnd_alpha(0:ND,Edge_Num,DDK,l) * vEdge(0:ND,Edge_Num,n,DDK) &
                        +   bnd_beta(0:ND,Edge_Num,DDK,l) * bEdge(0:ND,Edge_Num,DDK,l) &
                        *   dvdnEdge(0:ND,Edge_Num,n,DDK)
      
               end select
            enddo !Edge_Num
         enddo !index
      enddo ! DDK
      
      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1;ND2p = ND2 + 1
      
         do n = 1, (ND1p*ND2p)
      
            do Edge_Num=1,4
               select case (Edge_Num)
                  case (2,4)
                     ND=ND2; NDFix=0; if (Edge_Num .eq. 2) NDFix=ND1
      
                     select case (BC_Type(Edge_Num,DDK))
                        case(0) ! internal BC
                           DDK_Connect   = DM_Connect(1,Edge_Num,DDK)
                           Edge_Connect  = DM_Connect(2,Edge_Num,DDK)
                           Patch_Type    = DM_Connect(3,Edge_Num,DDK)

                           select case (Patch_Type)
                              case(1)
                                 do j=0,ND
                                    BITF1(j,Edge_Num,n,DDK) = &
                                       vEdge(j,Edge_Num,n,DDK) &
                                    -  vEdge(j,Edge_Connect,n,DDK_Connect)
      
                                    BITF2(j,Edge_Num,n,DDK) = &
                                       a(DDK) * dvdnEdge(j,Edge_Num,n,DDK) &
                                     + a(DDK_Connect) &
                                     * dvdnEdge(j,Edge_Connect,n,DDK_Connect)
                                 enddo
                              case (-1) ! Reverse Patching
                                 do j=0,ND
                                    BITF1(j,Edge_Num,n,DDK) = &
                                       vEdge(j,Edge_Num,n,DDK) &
                                     - vEdge(ND-j,Edge_Connect,n,DDK_Connect)
      
                                    BITF2(j,Edge_Num,n,DDK) = &
                                       a(DDK) * dvdnEdge(j,Edge_Num,n,DDK) &
                                     + a(DDK_Connect) &
                                     * dvdnEdge(ND-j,Edge_Connect,n,DDK_Connect)
                                 enddo
                           end select ! Patch_type
                     end select ! bc_type
                  case(1,3)
                     ND=ND1; NDFix=0; if (Edge_Num .eq. 3) NDFix=ND2
      
                     select case (BC_Type(Edge_Num,DDK))
                        case (0)
                           DDK_Connect  = DM_Connect(1,Edge_Num,DDK)
                           Edge_Connect = DM_Connect(2,Edge_Num,DDK)
                           Patch_Type   = DM_Connect(3,Edge_Num,DDK)

                           select case (Patch_Type)
                              case(1)
                                 do i=0,ND
                                    BITF1(i,Edge_Num,n,DDK) = &
                                       vEdge(i,Edge_Num,n,DDK) &
                                     - vEdge(i,Edge_Connect,n,DDK_Connect)
      
                                    BITF2(i,Edge_Num,n,DDK) = &
                                       a(DDK) * dvdnEdge(i,Edge_Num,n,DDK) &
                                     + a(DDK_Connect) &
                                     * dvdnEdge(i,Edge_Connect,n,DDK_Connect)
                                 enddo
                              case (-1) ! Reverse Patching
                                 do i=0,ND
                                    BITF1(i,Edge_Num,n,DDK) = &
                                       vEdge(i,Edge_Num,n,DDK) &
                                     - vEdge(ND-i,Edge_Connect,n,DDK_Connect)
      
                                    BITF2(i,Edge_Num,n,DDK) = &
                                       a(DDK) * dvdnEdge(i,Edge_Num,n,DDK) &
                                     + a(DDK_Connect) &
                                     * dvdnEdge(ND-i,Edge_Connect,n,DDK_Connect)
                                 enddo ! i
                           end select ! Patch_Type
                     end select ! BC_Type
               end select !Edge_Num
            enddo !Edge_Num
         enddo ! index
      enddo ! DDK
      
      end subroutine
!---------------------------------------------------------------
      subroutine Smoothing_Pack(LD1,LD2,l)

      use constants
      use Legendre
      use MD2D_Grid
      use State_Var
      use CG_Var
      use Multigrid_Var

      implicit none
      integer      :: N_vcycle, m_smooth, LD1, LD2
      integer      :: vcycle, iterNumc, method, ind_JS
      integer      :: k, i, j, n, ND1p, ND2p,l, ND1c, ND2c, Nx, Ny
      real(kind=8) :: smoothpar
      real(kind=8) :: x_vc_ini(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: error_vc(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp(0:LD1,0:LD2,1:TotNum_DM),tmp1(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp2(0:LD1,0:LD2,1:TotNum_DM),tmp3(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp4(0:LD1,0:LD2,1:TotNum_DM),tmp5(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp6(0:LD1,0:LD2,1:TotNum_DM)

      Nx = PolyDegN_DM(1,1,l); Ny = PolyDegN_DM(2,1,l)

      call alloc_mem_jacobismooth_var(PolyDegN_DM(1,1,l),TotNum_DM)

      N_vcycle  = 1
      m_smooth  = 1
      smoothpar = 2.d0/3.d0

!--------------------------------------------------------------------------------
!     Construct Jacobi-smoothing matrix
!      open(347,file='M.text')
      
      do DDK = 1, TotNum_DM
         ND1=PolyDegN_DM(1,DDK,l); ND2=PolyDegN_DM(2,DDK,l)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do j=0,ND2
            do i=0,ND1
               n = i + j*(ND2+1) + 1
               M(n) = 1.d0/ML(i,j,n,DDK)
      !      write(347,*)n,M(n)
            enddo
         enddo
      enddo
      !  close(347)
      
      call copy(r_smooth,rhs,Nx,Ny,l)
      call copy(b_smooth,r_smooth,Nx,Ny,l)
      
      x_vc_ini = 0.d0 ! initial gues

      do vcycle = 1, N_vcycle
         
         xc_in=0.0

         call copy(x_vc,x_vc_ini,Nx,Ny,l)
      
         call axhelm2(Nx,Ny,x_vc,Ax_vc,l)

         call add3s2(r_smooth,b_smooth,Ax_vc,1.0,-1.0,Nx,Ny,l)
         
         call chk_amax('xvc',x_vc,Nx,Ny,l)
         call chk_amax('bJr',r_smooth,Nx,Ny,l)

!        Jacobi smoothing
         do k = 1,m_smooth

            do DDK = 1, TotNum_DM
               ND1=PolyDegN_DM(1,DDK,l)
               ND2=PolyDegN_DM(2,DDK,l)
               do j=0,ND2
                  do i=0,ND1
                     n = i + j*(ND2+1) + 1
                     Mr(i,j,DDK) = M(n) * r_smooth(i,j,DDK)
                  enddo
               enddo
               x_vc(0:ND1,0:ND2,DDK) = x_vc(0:ND1,0:ND2,DDK) &
                                     + smoothpar * Mr(0:ND1,0:ND2,DDK)
            enddo !DDK
            call chk_amax('xas',x_vc,Nx,Ny,l)

            call axhelm2(Nx,Ny,x_vc,Ax_vc,l)
            call add3s2(r_smooth,b_smooth,Ax_vc,1.0,-1.0,Nx,Ny,l)

         enddo ! smooth
      
!        Coarse-grid restriction x <--- x + e, where e is approximated on coarse grid
      
         do DDK = 1 ,TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
            ND1c = PolyDegN_DM(1,DDK,l-1); ND2c = PolyDegN_DM(2,DDK,l-1)
      
            tmp(0:ND1,0:ND2c,DDK) = &
               Matmul(r_smooth(0:ND1,0:ND2,DDK), Ihy(0:ND2,0:ND2c,DDK))
      
            rc_smooth(0:ND1c,0:ND2c,DDK) = &
               Matmul(Ihx_transpose(0:ND1c,0:ND1,DDK),tmp(0:ND1,0:ND2c,DDK))
      
         enddo
         
         call chk_amax('rcs',rc_smooth(0:ND1c,0:ND2c,1:TotNum_DM),ND1c,ND2c,l-1)
      
         iterNumc = 0
         do DDK=1,TotNum_DM
            iterNumc = iterNumc + (PolyDegN_DM(1,DDK,l-1)+1) &
                                   * (PolyDegN_DM(2,DDK,l-1)+1)
         enddo

         call chk_amax('xci',xc_in,PolyDegN_DM(1,1,l-1),PolyDegN_DM(2,1,l-1),l-1)   

         call CG(ec(0:ND1c,0:ND2c,1:TotNum_DM),xc_in(0:ND1c,0:ND2c,1:TotNum_DM),&
         rc_smooth(0:ND1c,0:ND2c,1:TotNum_DM),1,iterNumc,1e-20)

         call chk_amax('ecg',ec,ND1c,ND2c,l-1)
         call chk_amax('rcc',rc_smooth,ND1c,ND2c,l-1)
      
         do DDK = 1 ,TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
            ND1c = PolyDegN_DM(1,DDK,l-1); ND2c = PolyDegN_DM(2,DDK,l-1)
      
            tmp1(0:ND1c,0:ND2,DDK) = &
               Matmul(ec(0:ND1c,0:ND2c,DDK),Ihy_transpose(0:ND2c,0:ND2,DDK))

            ef(0:ND1,0:ND2,DDK) = &
               Matmul(Ihx(0:ND1,0:ND1c,DDK),tmp1(0:ND1c,0:ND2,DDK))
         enddo
      
         call add2s2(x_vc,ef,1.0,Nx,Ny,l)

         call add3s2(error_vc,v,x_vc,1.0,-1.0,Nx,Ny,l)
     
         call chk_amax('ers',error_vc,Nx,Ny,l)
         call copy(x_vc_ini,x_vc,Nx,Ny,l)

!         do DDK =1, TotNum_DM
!         do j=0,Ny
!         do i=0,Nx
!         write(*,*)'x_vc',i,j,DDK,x_vc(i,j,DDK)
!         enddo
!         enddo
!         enddo
      
      enddo ! vcycle
      write(10,*)'End cycle'
      
      end subroutine Smoothing_Pack
      !===============================================================================
      subroutine Smoothing_Pack_Overlapping(LD1,LD2,l)
      use Legendre
      use MD2D_Grid
      use State_Var
      use CG_Var
      use Multigrid_Var

      implicit none
      integer :: N_vcycle, m_smooth, LD1, LD2
      integer :: vcycle, iterNumc, method, ind_JS, iterNum
      integer :: k, i, j, index, ND1p, ND2p,l, ND1c, ND2c, Nx, Ny
      real(kind=8) :: smoothpar
      real(kind=8) :: x_vc_ini(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: error_vc(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp(0:LD1,0:LD2,1:TotNum_DM),tmp1(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp2(0:LD1,0:LD2,1:TotNum_DM),tmp3(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp4(0:LD1,0:LD2,1:TotNum_DM),tmp5(0:LD1,0:LD2,1:TotNum_DM)
      real(kind=8) :: tmp6(0:LD1,0:LD2,1:TotNum_DM)

      Nx = PolyDegN_DM(1,1,l); Ny = PolyDegN_DM(2,1,l)

      N_vcycle = 1
      m_smooth = 1
      smoothpar = 1.0 !2.d0/3.d0
      

      call copy(r_smooth,rhs,Nx,Ny,l)
      call copy(b_smooth,r_smooth,Nx,Ny,l)
      
      x_vc_ini = 0.d0 ! initial gues
      
      do vcycle = 1, N_vcycle

         call copy(x_vc,x_vc_ini,Nx,Ny,l)

         call axhelm2(Nx,Ny,x_vc,Ax_vc,l)

         call add3s2(r_smooth,b_smooth,Ax_vc,1.0,-1.0,Nx,Ny,l)

         call chk_amax('xvc',x_vc,Nx,Ny,l)
         call chk_amax('bJr',r_smooth,Nx,Ny,l)
      
!     Additive Jacobi Schwartz
         do k = 1,m_smooth

            z_ov = 0*x_vc
      
      !-In Ind_JS, we are Multiplying the inverse of the BlockL operator.
            do ind_JS = 1,1
               do DDK = 1, TotNum_DM
                  ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
      
                  tmp2(0:ND1,0:ND2,DDK) = &
                Matmul(r_smooth(0:ND1,0:ND2,DDK), Sy_norm(0:ND2,0:ND2,DDK))
      
                  tmp3(0:ND1,0:ND2,DDK) = &
                Matmul(Sx_t_norm(0:ND1,0:ND1,DDK), tmp2(0:ND1,0:ND2,DDK))
      
                  do j=0,ND2
                     do i=0,ND1
                        index = i + j*(ND2+1) +1
                        tmp4(i,j,DDK) = Diagonal(index,DDK) * tmp3(i,j,DDK)
                     enddo
                  enddo

                  tmp5(0:ND1,0:ND2,DDK) = &
                    Matmul(tmp4(0:ND1,0:ND2,DDK), Sy_t_norm(0:ND2,0:ND2,DDK))
      
                  z_ov_sol(0:ND1,0:ND2,DDK) = &
                    Matmul(Sx_norm(0:ND1,0:ND1,DDK), tmp5(0:ND1,0:ND2,DDK))
      
                z_ov(0:ND1,0:ND2,DDK) = z_ov(0:ND1,0:ND2,DDK) + z_ov_sol(0:ND1,0:ND2,DDK)
      
               enddo! DDK
            enddo !ind_JS
      
            call add2s2(x_vc,z_ov,smoothpar,Nx,Ny,l)      

            call axhelm2(Nx,Ny,x_vc,Ax_vc,l)
            call add3s2(r_smooth,b_smooth,Ax_vc,1.0,-1.0,Nx,Ny,l)

         enddo ! smooth
      
            !-compute the error of x_vc in smoothing part and exact solution
      
!     Coarse-grid restriction x <--- x + e, where e is approximated on coarse grid
      
         do DDK = 1 ,TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,l); ND2 = PolyDegN_DM(2,DDK,l)
            ND1c = PolyDegN_DM(1,DDK,l-1); ND2c = PolyDegN_DM(2,DDK,l-1)
      
            tmp(0:ND1,0:ND2c,DDK) = &
               Matmul(r_smooth(0:ND1,0:ND2,DDK), Ihy(0:ND2,0:ND2c,DDK))
      
            rc_smooth(0:ND1c,0:ND2c,DDK) = &
               Matmul(Ihx_transpose(0:ND1c,0:ND1,DDK),tmp(0:ND1,0:ND2c,DDK))

         enddo
      
         iterNumc = 0; iterNum = 0
         do DDK=1,TotNum_DM
            iterNumc = iterNumc + (PolyDegN_DM(1,DDK,l-1)+1) &
                                * (PolyDegN_DM(2,DDK,l-1)+1)
            iterNum = iterNum + (PolyDegN_DM(1,DDK,l)+1) &
                              * (PolyDegN_DM(2,DDK,l)+1)
         enddo
      
      
!     Calling CG to solve for ec
         call CG(ec,xc_in,rc_smooth,l-1,iterNumc,1e-20)
      
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
         call add2s2(x_vc,ef,1.0,Nx,Ny,l)

         call add3s2(error_vc,v,x_vc,1.0,-1.0,Nx,Ny,l)

         call chk_amax('ers',error_vc,Nx,Ny,l)
         call copy(x_vc_ini,x_vc,Nx,Ny,l)
      
         enddo ! vcycle
          write(10,*)'End cycle'
      
      endsubroutine Smoothing_Pack_Overlapping
    
      !---------------------------------------------------------------------------
      subroutine fd_weights_full(xx,x,n,mm,c)
        use Multigrid_Var
        implicit none
      !
      !     This routine evaluates the derivative based on all points
      !     in the stencils.  It is more memory efficient than "fd_weights"
      !
      !     This set of routines comes from the appendix of
      !     A Practical Guide to Pseudospectral Methods, B. Fornberg
      !     Cambridge Univ. Press, 1996.   (pff)
      !
      !     Input parameters:
      !       xx -- point at wich the approximations are to be accurate
      !       x  -- array of x-ordinates:   x(0:n)
      !       n  -- polynomial degree of interpolant (# of points := n+1)
      !       m  -- highest order of derivative to be approxxmated at xi
      !
      !     Output:
      !       c  -- set of coefficients c(0:n,0:m).
      !             c(j,k) is to be applied at x(j) when
      !             the kth derivative is approxxmated by a
      !             stencil extending over x(0),x(1),...x(n).
      !
      !
        integer:: n, mm, i, j, k, mn
        real(kind=8):: xx
        real(kind=8):: x(0:n)
        real(kind=8):: c(0:n,0:mm)
        real(kind=8):: c1, c2, c3, c4, c5
      !
        c1       = 1.
        c4       = x(0) - xx
      !
        do k=0,mm
          do j=0,n
            c(j,k) = 0.
          enddo
        enddo
        c(0,0) = 1.
      !
        do i=1,n
          mn = min(i,mm)
          c2 = 1.
          c5 = c4
          c4 = x(i)-xx
          do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
              do k=mn,1,-1
              c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
              enddo
              c(i,0) = -c1*c5*c(i-1,0)/c2
                do k=mn,1,-1
                  c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
                enddo
                c(j,0) = c4*c(j,0)/c3
          enddo
          c1 = c2
        enddo
      
        return
      end subroutine
      !!-----------------------------------------------------------------------
      subroutine Interp_mat(l)
        use Legendre
        use MD2D_Grid
        use Multigrid_Var
        implicit none
      
        integer:: l, i, j, ND1c, ND2c 
      
        do DDK = 1, TotNum_DM
          ND1 = PolyDegN_DM(1,DDK,l)
          ND2 = PolyDegN_DM(2,DDK,l)
          xo(0:ND1,DDK) = LGLCoord(0:ND1,ND1)
          yo(0:ND2,DDK) = LGLCoord(0:ND2,ND2)
        enddo
      
        
        do DDK = 1, TotNum_DM
          ND1 = PolyDegN_DM(1,DDK,l-1)
          ND2 = PolyDegN_DM(2,DDK,l-1)
          xi(0:ND1,DDK) = LGLCoord(0:ND1,ND1)
          yi(0:ND2,DDK) = LGLCoord(0:ND2,ND2)
        enddo
      
      
        do DDK = 1, TotNum_DM
          ND1 = PolyDegN_DM(1,DDK,l)
          ND1c = PolyDegN_DM(1,DDK,l-1)
          ND2 = PolyDegN_DM(2,DDK,l)
          ND2c = PolyDegN_DM(2,DDK,l-1)
      
      !-Constructing the interpolation matrix J
          do i = 0, ND1
          call fd_weights_full(xo(i,DDK),xi(:,DDK),ND1c,1,wx(0:ND1c,1:2,DDK))
          Ihx(0:ND1c,i,DDK) = wx(0:ND1c,1,DDK)
      
          enddo
      
      !-Storing the interpolation matrix J and J'
          Ihx_transpose(:,:,DDK) = Ihx(:,:,DDK)
          Ihx(:,:,DDK) = transpose(Ihx(:,:,DDK))
      
      !-Constructing the interpolation matrix J
          do j = 0, ND2
          call fd_weights_full(yo(j,DDK),yi(:,DDK),ND2c,1,wy(0:ND2c,1:2,DDK))
          Ihy(0:ND2c,j,DDK) = wy(0:ND2c,1,DDK)
          enddo
      
      !-Storing the interpolation matrix J and J'
          Ihy_transpose(:,:,DDK) = Ihy(:,:,DDK)
          Ihy(:,:,DDK) = transpose(Ihy(:,:,DDK))
      
        enddo ! DDK
     
!        open(319,file='Jhx.text')
!        do DDK = 1, TotNum_DM
!          ND1 = PolyDegN_DM(1,DDK,l)
!          ND1c = PolyDegN_DM(1,DDK,l-1)
!          do j = 0,ND1c
!            do i =0,ND1
!            write(319,*)i,j,Ihx(i,j,DDK)
!            enddo
!          enddo
!        enddo
!     
!        open(320,file='Jhy.text')
!        do DDK = 1, TotNum_DM
!          ND2 = PolyDegN_DM(2,DDK,l)
!          ND2c = PolyDegN_DM(2,DDK,l-1)
!          do j = 0,ND2c
!            do i =0,ND2
!            write(320,*)i,j,Ihy(i,j,DDK)
!            enddo
!          enddo
!        enddo
!      
!        open(537,file='Jhx_transpose.text')
!        do DDK = 1, TotNum_DM
!          ND1 = PolyDegN_DM(1,DDK,1)
!          ND1c = PolyDegN_DM(1,DDK,2)
!          do j = 0,ND1
!            do i =0,ND1c
!            write(537,*)i,j,Ihx_transpose(i,j,DDK)
!            enddo
!          enddo
!        enddo
!      
!        open(668,file='Jhy_transpose.text')
!        do DDK = 1, TotNum_DM
!          ND2 = PolyDegN_DM(2,DDK,1)
!          ND2c = PolyDegN_DM(2,DDK,2)
!          do j = 0,ND2
!            do i =0,ND2c
!            write(668,*)i,j,Ihy_transpose(i,j,DDK)
!            enddo
!          enddo
!        enddo
      !
      
      end subroutine
      !================================================================================
      !-This subroutine is to construct the Lx and Ly operators which are in one
      ! dimensional. This idea is from the paper "Hybrid Multigrid/Schwarz
      ! Algorithms for the Spectral Element Method
      subroutine Construct_Lx_Ly_operator(level)
      use Multigrid_Var
      use State_Var
      use constants
      use Legendre
      use MD2D_Grid
      implicit none
      integer:: ND1p, ND2p, index, lbw, info
      integer:: i, j, level
      
!     Initialize variables for Lx and Ly operator
!      call alloc_mem_Lx_Ly_var(PolyDegN_DM(1,1,1),TotNum_DM)
      
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
         ND1p = ND1 + 1; ND2p = ND2 + 1

         dudx_tmp(0:ND1,0:ND1,DDK) = &
                   (Jacobin(0:ND1,0:ND2,DDK,level) &
                *     a_pts(0:ND1,0:ND2,DDK,level) ) &
                * (dxi1_dx1(0:ND1,0:ND2,DDK,level)**2.0) &
!                +  dxi1_dx2(0:ND1,0:ND2,DDK,level)**2.0 ) &
                *  Diff_xi1(0:ND1,0:ND1,ND1)
!         dudx_tmp(0:ND1,0:ND1,DDK) = &
!                     a_pts(0:ND1,0:ND2,DDK,level)  &
!                *  Diff_xi1(0:ND1,0:ND1,ND1) &
!                * (dxi1_dx1(0:ND1,0:ND2,DDK,level)) 
      
         dudy_tmp(0:ND2,0:ND2,DDK) = &
                   (Jacobin(0:ND1,0:ND2,DDK,level) &
                *     a_pts(0:ND1,0:ND2,DDK,level) ) &
!                * (dxi2_dx1(0:ND1,0:ND2,DDK,level)**2.0 &
                *  (dxi2_dx2(0:ND1,0:ND2,DDK,level)**2.0 ) &
                *  Diff_xi1(0:ND1,0:ND1,ND1)

!         dudy_tmp(0:ND2,0:ND2,DDK) = &
!                     a_pts(0:ND1,0:ND2,DDK,level)  &
!                *  Diff_xi1(0:ND1,0:ND1,ND1) &
!                *  (dxi2_dx2(0:ND1,0:ND2,DDK,level) )
!      
         Lx(0:ND1,0:ND1,DDK) = - Matmul(Diff_xi1(0:ND1,0:ND1,ND1), &
                                           dudx_tmp(0:ND1,0:ND1,DDK) )
      
         Ly(0:ND2,0:ND2,DDK) = - Matmul(Diff_xi1(0:ND1,0:ND1,ND1), &
                                           dudy_tmp(0:ND2,0:ND2,DDK) )
      enddo

!      write(10,*)'L operator without penalty'
!      do DDK=1,TotNum_DM
!         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
!         do j=0,ND2
!            do i=0,ND1
!               write(10,*)i,j,DDK,Lx(i,j,DDK),Ly(i,j,DDK)
!            enddo
!         enddo
!      enddo
      
!     Start adding penalty for every DDK
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
         ND1p = ND1 + 1; ND2p = ND2 + 1
!     side 4
         Edge_Num=4
         select case (BC_Type(Edge_Num,DDK))
            case(1)
!               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
!               do j=0,ND2
!                  do i=0,ND1
!                     write(10,*)i,j,tauD(i,j,Edge_Num,DDK,level),tau_tilde(i,j,Edge_Num,DDK,level)
!                  enddo
!               enddo
               Lx(0:ND1,0,DDK) = Lx(0:ND1,0,DDK) &
                               + (tauD(0:ND1,0,Edge_Num,DDK,level)   &
                               -  tau_tilde(0:ND1,1,Edge_Num,DDK,level))
            case(2,3)
               do j = 0,ND2
                  Lx(0:ND1,j,DDK) = Lx(0:ND1,j,DDK) &
                                  +    tauND1(0:ND1,j,Edge_Num,DDK,level) !&
               enddo
            case(0)
               write(10,*)'Edge:',Edge_Num,'interface penalty'
!     add interface penalty : first set of dirichlet
               Lx(0:ND1,0,DDK) = Lx(0:ND1,0,DDK) &
                               + tau1(0:ND1,0,Edge_Num,DDK,level) 
               Lx(0,0,DDK) = Lx(0,0,DDK) &
                           +  tau2(1,Edge_Num,DDK,level) !&
!     add interface penalty : second set of dirichlet
               Lx(0,0,DDK) = Lx(0,0,DDK) &
                           + Sigma_tild(1,Edge_Num,DDK,level)  !&
!     add interface penalty : neumann
               Lx(0,0:ND2,DDK) = Lx(0,0:ND2,DDK) &
                               - SigmaHat(0,Edge_Num,DDK,level) *2&
                               * Diff_xi1(0,0:ND1,ND1)
         end select
!     side 2
         Edge_Num=2
         select case (BC_Type(Edge_Num,DDK))
            case(1)
!               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
!               do j=0,ND2
!                  do i=0,ND1
!                     write(10,*)i,j,tauD(i,j,Edge_Num,DDK,level),tau_tilde(i,j,Edge_Num,DDK,level)
!                  enddo
!               enddo
               Lx(0:ND1,ND2,DDK) = Lx(0:ND1,ND2,DDK) &
                                 +(    tauD(0:ND1,ND2,Edge_Num,DDK,level) &
                                 - tau_tilde(0:ND1,ND2-1,Edge_Num,DDK,level))
            case(2,3)
               do j=0,ND2
                  Lx(0:ND1,j,DDK) = Lx(0:ND1,j,DDK) &
                                  +    tauND1(0:ND1,j,Edge_Num,DDK,level) !&
               enddo
            case(0)
!               write(10,*)'Edge:',Edge_Num,'interface penalty'
!               do j=0,ND2
!                  do i=0,ND1
!                     write(10,*)i,j,tau1(i,j,Edge_Num,DDK,level)
!                  enddo
!               enddo
!               do i=0,ND1
!                  write(10,*)i,Sigma_tild(i,Edge_Num,DDK,level),tau2(i,Edge_Num,DDK,level),SigmaHat(i,Edge_Num,DDK,level)
!               enddo
!     add interface penalty : first set of dirichlet
               Lx(0:ND1,ND1,DDK) = Lx(0:ND1,ND1,DDK) &
                                 + tau1(0:ND1,ND1,Edge_Num,DDK,level) !
               Lx(ND1,ND1,DDK) = Lx(ND1,ND1,DDK) &
                               +     tau2(ND1-1,Edge_Num,DDK,level)  !&
!     add interface penalty : second set of dirichlet
               Lx(ND1,ND1,DDK) = Lx(ND1,ND1,DDK) &
                                  + Sigma_tild(ND1-1,Edge_Num,DDK,level)  !&
!     add interface penalty : neumann
               Lx(ND1,0:ND2,DDK) = Lx(ND1,0:ND2,DDK) &
                                    + SigmaHat(ND2,Edge_Num,DDK,level) *2 &
                                    * Diff_xi1(ND1,0:ND1,ND1)
!            do i=0,ND1
!               write(10,*)SigmaHat(ND2,Edge_Num,DDK,level) * Diff_xi1(ND1,i,ND1)
!            enddo
         end select

            write(10,*)'Lx operator after adding penalty'
!            do j=0,ND1
!               do i =0,ND1
!                  write(10,*)i,j,Lx(i,j,DDK)
!               enddo
!            enddo
!     side 1
         Edge_Num=1
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
!               do j=0,ND2
!                  do i=0,ND1
!                     write(10,*)i,j,tauD(i,j,Edge_Num,DDK,level),tau_tilde(i,j,Edge_Num,DDK,level)
!                  enddo
!               enddo
               Ly(0:ND1,0,DDK) = Ly(0:ND1,0,DDK) &
                               + (    tauD(0,0:ND2,Edge_Num,DDK,level) &
                               - tau_tilde(1,0:ND2,Edge_Num,DDK,level))
            case(2,3)
               do i=0,ND1
                  Ly(i,0:ND2,DDK) = Ly(i,0:ND2,DDK) &
                                  +    tauND1(i,0:ND2,Edge_Num,DDK,level) !&
               enddo
            case(0)
!     Add interface penalty : first set of dirichlet
               Ly(0:ND1,0,DDK) = Ly(0:ND1,0,DDK) &
                               + tau1(0,0:ND2,Edge_Num,DDK,level) 
               Ly(0,0,DDK) = Ly(0,0,DDK) &
                           + tau2(1,Edge_Num,DDK,level) !&
!     Add interface penalty : second set of dirichlet
               Ly(0,0,DDK) = Ly(0,0,DDK) &
                           + Sigma_tild(1,Edge_Num,DDK,level)  !&
!     Add interface penalty : neumann
               Ly(0,0:ND2,DDK) = Ly(0,0:ND2,DDK) &
                               - SigmaHat(0,Edge_Num,DDK,level) *2 &
                               * Diff_xi1(0,0:ND2,ND2)
         end select
!     side 3
         Edge_Num=3
         select case (BC_Type(Edge_Num,DDK))
            case(1)
               write(10,*)'Edge:',Edge_Num,'dirichlet penalty'
!               do j=0,ND2
!                  do i=0,ND1
!                     write(10,*)i,j,tauD(i,j,Edge_Num,DDK,level),tau_tilde(i,j,Edge_Num,DDK,level)
!                  enddo
!               enddo
               Ly(0:ND1,ND2,DDK) = Ly(0:ND1,ND2,DDK) &
                                 + (    tauD(ND1,0:ND2,Edge_Num,DDK,level)  &
                                 - tau_tilde(ND1-1,0:ND2,Edge_Num,DDK,level))  !&
            case(2,3)
               do i=0,ND1
                  Ly(i,0:ND2,DDK) = Ly(i,0:ND2,DDK) &
                                  +    tauND1(i,0:ND2,Edge_Num,DDK,level) !&
               enddo
            case(0)
!     Add interface penalty : first set of dirichlet
               Ly(0:ND1,ND2,DDK) = Ly(0:ND1,ND2,DDK) &
                                 + tau1(ND1,0:ND2,Edge_Num,DDK,level) 
               Ly(ND1,ND2,DDK) = Ly(ND1,ND2,DDK) &
                               + tau2(ND1-1,Edge_Num,DDK,level) !&
!     Add interface penalty : second set of dirichlet
               Ly(ND1,ND2,DDK) = Ly(ND1,ND2,DDK) &
                               + Sigma_tild(ND1-1,Edge_Num,DDK,level) !&
!     Add interface penalty : neumann
               Ly(ND1,0:ND2,DDK) = Ly(ND1,0:ND2,DDK) &
                                 + SigmaHat(ND2,Edge_Num,DDK,level) *2 &
                                 * Diff_xi1(ND2,0:ND2,ND2)
         end select
      enddo !DDK

!      do DDK=1,TotNum_DM
!         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
!         ND1p = ND1 + 1; ND2p = ND2 + 1
!         do j = 0, ND2
!            do i = 0, ND1
!               Lx(i,j,DDK) = Lx(i,j,DDK) !* dxi1_dx1(i,j,DDK,level) 
!               Ly(i,j,DDK) = Ly(i,j,DDK) !* dxi2_dx2(i,j,DDK,level) 
!            enddo
!         enddo
!      enddo
!      do DDK=1,TotNum_DM
!         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
!         ND1p = ND1 + 1; ND2p = ND2 + 1
!         do j = 0, ND2
!            do i = 0, ND1
!               Ly(i,j,DDK) = Lx(i,j,DDK) 
!            enddo
!         enddo
!      enddo

!     Multiply mass matrix in x,y direction
      do DDK=1,TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
         ND1p = ND1 + 1; ND2p = ND2 + 1
         do j = 0, ND2
            do i = 0, ND1
               Lx(i,j,DDK) = Lx(i,j,DDK) * LGLWeights(i,ND1)
               Ly(i,j,DDK) = Ly(i,j,DDK) * LGLWeights(i,ND2)
            enddo
         enddo
      enddo

!     Construct Matrix Bx and By which will be used in solving the generalized eigenvalue problem Ax = lambda B x
      do DDK = 1, TotNum_DM
         ND1=PolyDegN_DM(1,DDK,level); ND2=PolyDegN_DM(2,DDK,level)
         do i = 0, ND1
            Bx(i,i,DDK) = LGLWeights(i,ND1)
         enddo
         do j= 0, ND2
            By(j,j,DDK) = LGLWeights(j,ND2)
         enddo
      enddo
      
!     Checking whether Lx operator is SPD or not
      open(637,file='Lx.text')
         do DDK = 1, TotNum_DM
            ND1 = PolyDegN_DM(1,DDK,level)
            do j=0,ND1
               do i =0,ND1
                  write(637,*)i,j,Lx(i,j,DDK)
               enddo
            enddo
         enddo
      close(637)

!     Checking whether Ly operator is SPD or not
      open(638,file='Ly.text')
         do DDK = 1, TotNum_DM
            ND2 = PolyDegN_DM(2,DDK,level)
            do j=0,ND2
               do i =0,ND2
                  write(638,*)i,j,Ly(i,j,DDK)
               enddo
            enddo
         enddo
      close(638)

      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level)
         lbw =  (ND1+1) * (ND1+1)
         call dsygv(1,'V','U',ND1+1,Lx(0:ND1,0:ND1,DDK),&
               ND1+1,Bx(0:ND1,0:ND1,DDK),ND1+1,lamx(0:ND1,DDK),bwx,lbw,info)
      enddo
      write(10,*)'Complete diagonalize Lx operator'
      
!     Storing the eigenvectors that are not normalized yet
      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level)
         Sx_t(0:ND1,0:ND1,DDK) = transpose(Lx(0:ND1,0:ND1,DDK))
         Sx(0:ND1,0:ND1,DDK) = Lx(0:ND1,0:ND1,DDK)
      enddo
!     Start normalizing the eigenvectors of Lx
      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level)
         do j=0,ND1
            do i=0,ND1
               BSx(i,j,DDK) = LGLWeights(i,ND1) * Sx(i,j,DDK)
            enddo
         enddo
!     We need the diagonal part of Sx' B Sx to normalize the eigenvector
         SxBSx(0:ND1,0:ND1,DDK) = Matmul(Sx_t(0:ND1,0:ND1,DDK), &
                                          BSx(0:ND1,0:ND1,DDK))
         do j=0,ND1
            do i =0,ND1
            Sx_norm(i,j,DDK) = Sx(i,j,DDK) / sqrt(SxBSx(j,j,DDK))
!            write(10,*)sqrt(SxBSx(i,j,DDK)),BSx(i,j,DDK)
            enddo
         enddo
         Sx_t_norm(0:ND1,0:ND1,DDK) = transpose(Sx_norm(0:ND1,0:ND1,DDK))
      enddo
      write(10,*)'Complete normalize the eigenvector of Lx operator'
      
!     Start Diagonalize Ly operator
      do DDK = 1, TotNum_DM
         ND2 = PolyDegN_DM(2,DDK,level)
         lbw = 4 * (ND2+1) * (ND2+1)
         call dsygv(1,'V','U',ND2+1,Ly(0:ND2,0:ND2,DDK),&
               ND2+1,By(0:ND2,0:ND2,DDK),ND2+1,lamy(0:ND2,DDK),bwy,lbw,info)
      enddo
      write(10,*)'Complete diagonalize Ly operator'
      
!     Storing the eigenvectors that are not normalized yet
      do DDK = 1, TotNum_DM
         ND2 = PolyDegN_DM(2,DDK,level)
         Sy_t(0:ND2,0:ND2,DDK) = transpose(Ly(0:ND2,0:ND2,DDK))
         Sy(0:ND2,0:ND2,DDK) = Ly(0:ND2,0:ND2,DDK)
      enddo

!     Start normalizing the eigenvectors of Lx
      do DDK = 1, TotNum_DM
         ND2 = PolyDegN_DM(2,DDK,level)
         do j=0,ND2
            do i=0,ND2
               BSy(i,j,DDK) = LGLWeights(i,ND2) * Sy(i,j,DDK)
            enddo
         enddo
!     We need the diagonal part of Sx' B Sx to normalize the eigenvector
         SyBSy(0:ND2,0:ND2,DDK) = Matmul(Sy_t(0:ND2,0:ND2,DDK), &
                                          BSy(0:ND2,0:ND2,DDK))
         do j=0,ND2
            do i =0,ND2
               Sy_norm(i,j,DDK) = Sy(i,j,DDK) / sqrt(SyBSy(j,j,DDK))
            enddo
         enddo
         Sy_t_norm(0:ND1,0:ND1,DDK) = transpose(Sy_norm(0:ND1,0:ND1,DDK))
      enddo
      write(10,*)'Complete normalize the eigenvector of Ly operator'

!      write(*,*)'Output the normalizd eigenvectors of L_x and L_y'
!      do DDK = 1, TotNum_DM
!         ND1 = PolyDegN_DM(1,DDK,level)
!         do j=0,ND1
!            do i =0,ND1
!               write(10,*)i,j,Sx_norm(i,j,DDK),Sy_norm(i,j,DDK)&
!               ,Sx(i,j,DDK),Sy(i,j,DDK),Lx(i,j,DDK)
!            enddo
!         enddo
!      enddo
!      do DDK = 1, TotNum_DM
!         ND1 = PolyDegN_DM(1,DDK,level)
!         do j=0,ND1
!            write(10,*)j,sqrt(SxBSx(j,j,DDK)),sqrt(SyBSy(j,j,DDK))
!         enddo
!      enddo

      open(unit=50,file='Diagonal.text')
      do DDK = 1, TotNum_DM
         ND1 = PolyDegN_DM(1,DDK,level); ND2 = PolyDegN_DM(2,DDK,level)
         do j =0,ND2
            do i = 0,ND1
               index = i + j*(ND2+1) + 1
               Diagonal(index,DDK) = 1.0 / ( lamx(i,DDK) + lamy(j,DDK))
               write(50,*)lamx(i,DDK),lamy(j,DDK)
               write(50,*)i,j,DDK,index,Diagonal(index,DDK)
            enddo
         enddo
      enddo
      close(50)

      !!-Checking......
      !open(1054,file='Compare.text')
      !do DDK = 1, TotNum_DM
      !ND1 = PolyDegN_DM(1,DDK,level)
      !ND2 = PolyDegN_DM(2,DDK,level)
      !
      !tmp2(0:ND1,0:ND2,DDK) = &
      !Matmul(F(0:ND1,0:ND2,DDK), Sy(0:ND2,0:ND2,DDK))
      !
      !tmp3(0:ND1,0:ND2,DDK) = &
      !Matmul(Sx_t(0:ND1,0:ND1,DDK), tmp2(0:ND1,0:ND2,DDK))
      !
      !do j=0,ND2
      !do i=0,ND1
      !index = i + j*(ND2+1) +1
      !tmp4(i,j,DDK) = Diagonal(index,DDK) * tmp3(i,j,DDK)
      !enddo
      !enddo
      !tmp5(0:ND1,0:ND2,DDK) = &
      !Matmul(tmp4(0:ND1,0:ND2,DDK), Sy_t(0:ND2,0:ND2,DDK))
      !
      !test(0:ND1,0:ND2,DDK) = &
      !Matmul(Sx(0:ND1,0:ND1,DDK), tmp5(0:ND1,0:ND2,DDK))
      !
      !do j=0,ND2
      !do i=0,ND1
      !write(1054,*)i,j,test(i,j,DDK),v(i,j,DDK),test(i,j,DDK)-v(i,j,DDK)
      !enddo
      !enddo
      !enddo
      
      
      
      
      end subroutine
