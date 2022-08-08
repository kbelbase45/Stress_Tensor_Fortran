SUBROUTINE FINL_ELECT_STR(lm,lmmax,coulomb_str_int,E_xc_m_V_xc,exc_vxc) 

!!===================================================================================================
!! This program was written by Kamal Belbase to complete PhD project and tested extensively.         |
!! This subroutine computes the electrostatic stress tensor. The electrostatic stress tensor         |
!! is the sum of the Coulomb stress, the Madelung stress. In addition, the exchange correlation      |
!! stress is also calculated here. For details about them, please read the following manuscript      |
!!                                                                                                   |
!! Kamal Belbase, Andreas Tröster, and Peter Blaha                                                   |
!! Phys. Rev. B 104, 174113 – Published 29 November 2021                                             |
!! In the following if equation numbers are given then these equations can be found in the above     |
!! paper                                                                                             |
!!                                                                                                   |
!! The variables are                                                                                 |
!! lm                L angular and M magnetic qauntum number of charge density and potential         |
!! lmmax             total number of lm component for a given system                                 |
!!                                                                                                   |
!! coulomb_str_int   The columb stress tensor in the interestitial region, the first line of Eq.(G5) |
!!                                                                                                   |
!! E_xc_m_V_xc       The exchange correlation stress tensor calculated in the interestitial region   |
!!                   interestitial part of Eq. (G7), calculate using the fast Fourier Transform      |
!!                                                                                                   |                     
!! exc_vxc           The exchange correlation stress tensor inside the atomic region                 |
!!                   Eq.(G7) inside an atomic sphere => \int_0^{R_a}.                                                |
!!                                                                                                   |
!! In the following                                                                                  |
!! clm2  = Total electronic charge density, rho in formula and read from .clmsum extension file      |
!! vpot  = Total potential, sum of spherical (lm = 00) and non spherical potential (higher lm). They |
!!         are read from .vsp and .vns extension file.                                               |
!!====================================================================================================


!Get the necessary variables from the following modules
   USE densit ,    ONLY : lmmaxx,clm
   USE work   ,    ONLY : V, exc, vxc
   USE struct ,    ONLY : zz, mult, jri,RMT,dx,iatnr,nat,r0,vol,irel
   USE parallel    
   USE gga_str_kb, ONLY : gga_fin_str
   
   IMPLICIT NONE
   INCLUDE 'param.inc'
  
!Take the required variables for the message passing parallelization run  
#ifdef Parallel
   INCLUDE 'mpif.h'
#endif

   INTEGER,INTENT(IN)   :: lm(2,ncom+3,nat),lmmax(nat)
   REAL*8, INTENT(IN)   :: E_xc_m_V_xc
   
!! Two EXTERNAL functions are used to compute gaunt numbers i.e the surface 
!! integral of product of three spherical harmonics

   INTEGER, EXTERNAL    :: NOTRI
   REAL*8 , EXTERNAL    :: GAUNT
   LOGICAL              :: HEXA, CUBIC
   
   COMPLEX*16           :: coulomb_str_int(1:9), fc(ncom+3,nat), clm2(1:nrad,ncom+3,nat),                        &
                           vpot(1:nrad,ncom+3,nat), value_c(1:nrad), cabt(3,-1:1), ctest,                        &
                           zeroc, out_c
                           
   REAL*8               :: exc_vxc(1:9), exc_vxc_buf(1:9), VALUE(1:nrad),r(1:nrad),col_buf_pot(1:9),              &
                           col_pot_1(1:9), madelung_1(1:9),madelung_3(1:9), madelung_1_buf(nat,1:9),              &
                           madelung_2_buf(nat,1:9),madelung_3_buf(nat,1:9), madelung_final(1:9),                  &                           
                           madelung_final_sum(1:9), tot_stress(1:3,1:3) 
   
   REAL*8               :: out1,out2,gaunt1, gaunt2, test_col_sum, zero, pi, RMT1,RMT_i, sqrt_4pi, sqrt_4pii,     &
                           fopi_i, dx1, pone                        
   
   !Most of the single-letter integers (bad practice :-) are the same as those used in formulas in my PhD thesis.
   INTEGER              :: lmtot(2,ncom+3), lmmax_2, lm1, lm1p, m0, l0, l0p, m0p, t, tp, v1, s, lmtot1, alpha,    & 
                           beta, j, index2,jatom,mu,imax, lmmax1(nat), index_a                                                    


   !Define constants
   zeroc     = (0.d0,0.d0)
   zero      = 0.d0
   pi        = 4.0d+0*ATAN(1.d0)
   fopi_i    = 1.d0/( 4.d0*pi )
   sqrt_4pi  = SQRT( 4.d0*pi )
   sqrt_4pii = 1.d0/SQRT( 4.d0*pi )
   
   !Intialized the variables
   
   !For the Coulomb part, the first term of the second line of Ed.(G5)
   col_pot_1      = zero   
   
   !Madelung part of electrostatic, the second term of the second line of Ed.(G5)
   madelung_1     = zero   
   madelung_3     = zero   
   madelung_1_buf = zero
   madelung_2_buf = zero
   madelung_3_buf = zero
   madelung_final = zero
   
   !Calculate the c-matrix given in equations (C2) and (C6) in the KB paper
   CALL c_alpha_m(cabt)                     
   
   !Combine the index for the same and a different type of atom
   index_a = 0

   DO 29 jatom = nstart(myid), nstop(myid)                                               !! Loop over the different types of atoms
   DO 30 mu    = 1, mult(jatom)                                                          !! Loop over the same kind of atoms
   
         index_a   = index_a + 1
         RMT1      = RMT(jatom)                                                          !! Sphere boundary depending on different atoms, R_a in Eq.(G5)
         RMT_i     = 1.d0/RMT1   
         dx1       = dx(jatom)                                                           !! Radial mesh increment
         imax      = jri(jatom)                                                          !! Maximum points in a radial mesh
         HEXA      = .FALSE.
         CUBIC     = .FALSE.
         
         !For cubic or non-cubic lattice iatnr(jatom) is specified differently in the WIEN2k input file
         IF(iatnr(jatom) .LT. 0) HEXA  = .TRUE.
         IF(iatnr(jatom) .GT. 0) CUBIC = .TRUE.
         
         !Create logarithmic radial grid points within an atomic sphere
         DO J = 1, imax
            r(j) = r0(jatom)*exp( dx(jatom)*dble(j-1) )
         ENDDO
   
   !The electronic charge density and the potential are stored as a complex number
   !For lm=(0,0), normalization constant sqrt( 4.d0*pi ) is handled differently
   
   DO lm1 = 1, lmmax(jatom)     
      clm2(1:imax,lm1,jatom)                = clm(1:imax,lm1,jatom)
      vpot(1:imax,lm1,jatom)                = V(1:imax,lm1,jatom)
      IF(lm1 .EQ. 1) clm2(1:imax,lm1,jatom) = clm2(1:imax,lm1,jatom)*sqrt_4pii            
   ENDDO
   
   !Preprocessing for the non-cubic unit cell system
   IF(iatnr(jatom) .LT. 0) THEN         
     lmmax_2 = lmmax(jatom)
     
     DO lm1 = 1, lmmax(jatom)
        M0 = lm(2,lm1,jatom)
        
        IF(M0 .NE. 0) THEN  
        !rho/pot_{1,1} = rho/pot_{1,-1} and so on for m .not. 0
          lmmax_2 = lmmax_2 + 1          
          IF(lmmax_2 .GT. ncom+3) STOP  'ERROR IN "finl_elect_str.f" FOR LM COMBINATION'
          clm2(1:imax,lmmax_2,jatom) = clm2(1:imax,lm1,jatom)
          vpot(1:imax,lmmax_2,jatom) = vpot(1:imax,lm1,jatom)
          
        ENDIF
     ENDDO
    
    !Calculate multiplicative constants when converting complex spherical harmonics to real values
    CALL MULTFC(fc,jatom,lmmax_2,lm(1,1,jatom),iatnr(jatom))    
    
    DO lm1 = 1, lmmax_2       
       clm2(1:imax,lm1,jatom) = clm2(1:imax,lm1,jatom)*fc(lm1,jatom)      
       vpot(1:imax,lm1,jatom) = vpot(1:imax,lm1,jatom)*fc(lm1,jatom)      
    ENDDO
    
    !Combine -l and +l of the same m in the real spherical harmonics
    lmtot1 = 0
    DO 39 lm1 = 1, lmmax_2
          lmtot1 = lmtot1 + 1
          lmtot(1,lmtot1) = lm(1,lm1,jatom)
          lmtot(2,lmtot1) = lm(2,lm1,jatom)
          clm2(1:imax,lmtot1,jatom) = clm2(1:imax,lm1,jatom)
          vpot(1:imax,lmtot1,jatom) = vpot(1:imax,lm1,jatom)
          
          IF(lmtot1 .EQ. 1) GOTO 39
          IF( ( IABS( lmtot(1,lmtot1)) .EQ. IABS( lmtot(1,lmtot1-1))  ) .AND.   &
              (  lmtot(2,lmtot1) .EQ.  lmtot(2,lmtot1-1)  )  ) THEN
              
              clm2(1:imax,lmtot1-1,jatom) = clm2(1:imax,lmtot1-1,jatom) + clm2(1:imax,lmtot1,jatom)            
              vpot(1:imax,lmtot1-1,jatom) = vpot(1:imax,lmtot1-1,jatom) + vpot(1:imax,lmtot1,jatom)
              lmtot1 = lmtot1 - 1
              
          ENDIF    
          
39  CONTINUE  
 
    ELSE     
       !Preprocessing for the cubic unit cell system       
      
       CALL LM_COMBINE(imax,lmmax(jatom),clm2(1,1,jatom),lm(1,1,jatom),jatom)
       CALL LM_COMBINE(imax,lmmax(jatom),vpot(1,1,jatom),lm(1,1,jatom),jatom)              
       
       lmtot1 = lmmax(jatom)
       DO lm1 = 1, lmmax(jatom)          
          lmtot(1,lm1) =  lm(1,lm1,jatom)
          lmtot(2,lm1) =  lm(2,lm1,jatom)                                        
          IF( lm(2,lm1,jatom) .ne. 0) THEN
            lmtot1 = lmtot1 + 1
            IF( lmtot1 .gt. ncom+3 ) STOP 'ERROR IN finl_elect_str.F '
            lmtot(1,lmtot1) =  lm(1,lm1,jatom)
            lmtot(2,lmtot1) = -lm(2,lm1,jatom)
         
            clm2(1:imax,lmtot1,jatom) = clm2(1:imax,lm1,jatom)           !uncommenting this one leads to the perfect unstymmetry
            vpot(1:imax,lmtot1,jatom) = vpot(1:imax,lm1,jatom)
          ENDIF
       ENDDO
                             
    ENDIF    
    
    ! Calculate the Coulomb stress tensor within the atomic sphere, the first term in the second line of Eq.(G5)
    ! The Gaunt number comes after expanding the charge density and potential with respect to spherical harmonics.
    ! The product of two gaunt numbers are due to the surface integral of the four spherical harmonics after 
    ! expanding the charge density and the potential.
    
     col_buf_pot = zero     
     ctest       = zero
     
     DO lm1 = 1, lmtot1             
        l0 = lmtot(1,lm1)
        m0 = lmtot(2,lm1)
        pone = (-1.0)**m0
        
        DO lm1p = 1, lmtot1
           l0p = lmtot(1,lm1p)
           m0p = lmtot(2,lm1p)
           
           !Calculate integrand of the first term of the second line of Eq.(G5) 
           DO j = 1, imax
             value_c(j) = (clm2(j,lm1,jatom))*vpot(j,lm1p,jatom)                          
           ENDDO
           
           !Compute the integration of the real part of the above integrand 
           VALUE(:) = REAL( value_c(:) )
           CALL charge(r,dx1,VALUE,1,imax,out1)
           
           !Compute the integration of the complex part of the above integrand 
           VALUE(:) = AIMAG( value_c(:) )
           CALL charge(r,dx1,VALUE,1,imax,out2)
           
           !Combine the real and complex part
           out_c = dcmplx(out1,out2)
           
           !Check the condition that gives the non-zero contribution
           DO s = 0,2
              !Check the triangle inequality conditon 
              IF( NOTRI(s,1,1) .lt. 0    ) CYCLE
              IF( NOTRI(l0,l0p,s) .lt. 0 ) CYCLE
              
              DO v1 = -s,s
                 IF( (m0+m0p+v1) .ne. 0 ) CYCLE                 
                 DO t = -1,1
                    DO tp = -1,1
                       IF( (-v1+t+tp) .ne. 0 ) CYCLE
                       gaunt1 = GAUNT(s,1,1,v1,t,tp)
                       gaunt2 = GAUNT(l0,l0p,s,-m0,m0p,v1)                                               
                       
                       index2 = 0
                       DO alpha = 1, 3
                          DO beta = 1, 3
                             index2 = index2 + 1
                             col_buf_pot(index2) = col_buf_pot(index2) - pone*    &
                                                   cabt(alpha,t)*cabt(beta,tp)*   &                                                    
                                                   gaunt1*gaunt2*out_c              
                          ENDDO
                       ENDDO
                       
                    ENDDO  !tp
                    
                 ENDDO     !t
                 
              ENDDO        !v1
              
           ENDDO           !s
           
        ENDDO              !lm1p
        
     ENDDO                 !lm1        
     
     ! This line is for debugging purpose
     ! WRITE(6,*) 'col,ctest',col_buf_pot(2),ctest                
     
     !Symmetrization for the equivalent atoms using the rotloc and rotij 
     CALL SYMMETRY_STRESS_ROTIJ(index_a,jatom,col_buf_pot)
          col_pot_1 = col_pot_1 + col_buf_pot                  
     
    ! Calculate Madelung Stress Tensor, the last term of Eq(G5)
    ! When the expression of VM is put in the last term of Eq(G5), there will be three
    ! different contributions that are calcualted here.
    
     madelung_1 = zero
     DO lm1 = 1, lmtot1        
        l0 = lmtot(1,lm1)
        m0 = lmtot(2,lm1)
        pone = (-1)**m0
        
          DO J = 1, imax
             value_c(j) = clm2(j,lm1,jatom)*(1.d0/r(j)-RMT_i) 
          ENDDO
          
          !Be careful, dx from module is only defined for the first atom, so use dx1
          VALUE(:) = REAL( value_c(:) )          
          CALL charge(r,dx1,VALUE,1,imax,out1)
          
          VALUE(:) = AIMAG( value_c(:) )
          CALL charge(r,dx1,VALUE,1,imax,out2)   
          
          out_c = dcmplx(out1,out2)                        
              
        IF(NOTRI(l0,1,1) .lt. 0) CYCLE
        DO t = -1,1
           DO tp = -1,1
              IF( (m0+t+tp) .NE. 0) CYCLE                          
              gaunt1 = GAUNT(l0,1,1,-m0,t,tp)                           
                           
              index2 = 0
              DO alpha = 1,3
                 DO beta = 1,3
                    index2 = index2 + 1                                        
                    
                   madelung_1(index2) = madelung_1(index2) +                       &
                                        pone*cabt(alpha,t)*cabt(beta,tp)*          &
                                        gaunt1*out_c  
                   
                 ENDDO
              ENDDO
              
           ENDDO      !tp
        ENDDO         !t
        
     ENDDO            !lm1
     
     
                    
     DO t = -1,1
        DO tp = -1,1
           !Condition check for the Kronecker delta function
           IF( -t .NE. tp) CYCLE
           index2 = 0
           
           DO alpha = 1,3
              DO beta = 1,3
                 index2 = index2 + 1
                                  
                 madelung_2_buf(jatom,index2) = madelung_2_buf(jatom,index2) +        &
                                               (-1)**t*cabt(alpha,t)*cabt(beta,tp)*   &
                                               ZZ(jatom)*fopi_i*RMT_i
              ENDDO
           ENDDO
           
        ENDDO
     ENDDO
     
     
        
     madelung_3 = zero
     DO lm1 = 1, lmtot1          
        l0   = lmtot(1,lm1)
        m0   = lmtot(2,lm1)        
        pone = (-1)**m0
        IF( NOTRI(l0,1,1) .lt. 0 ) CYCLE
        DO t = -1,1
           DO tp = -1,1
              IF( (m0+t+tp) .ne. 0 ) CYCLE
              gaunt1 = GAUNT(l0,1,1,-m0,t,tp)

              index2 = 0
              DO alpha = 1,3
                 DO beta = 1,3
                    index2 = index2 + 1
                    madelung_3(index2) = madelung_3(index2) +                       &
                                         0.5d0*pone*cabt(alpha,t)*cabt(beta,tp)*    &
                                         vpot(imax,lm1,jatom)*gaunt1*fopi_i                                        
                 ENDDO
              ENDDO

           ENDDO
        ENDDO
     ENDDO   
     
     
   
     CALL SYMMETRY_STRESS_ROTIJ(index_a,jatom,madelung_1)
          madelung_1_buf(jatom,:) = madelung_1_buf(jatom,:) + madelung_1(:)

     CALL SYMMETRY_STRESS_ROTIJ(index_a,jatom,madelung_3)
          madelung_3_buf(jatom,:) = madelung_3_buf(jatom,:) + madelung_3(:)           

          
30  CONTINUE
29  CONTINUE
   
!This line was added when an error was found in the calculation of NaCl, MgS, GaAs
   DO jatom = nstart(myid), nstop(myid)         
      DO index2 = 1, 9                  
      
        madelung_final(index2) = madelung_final(index2)                   +     &
                                 ZZ(jatom)*( madelung_1_buf(jatom,index2) +     &
                                 madelung_2_buf(jatom,index2)             +     &
                                 madelung_3_buf(jatom,index2) )                  
      ENDDO
   ENDDO 


!! mpi_reduce combines, specified by MPI_SUM, the same contributions that are calculated 
!! by different processor for the same atom   

#ifdef Parallel
   CALL mpi_reduce(col_pot_1,col_buf_pot,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   CALL mpi_reduce(madelung_final,madelung_final_sum,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   CALL mpi_reduce(exc_vxc,exc_vxc_buf,9,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   
   IF(myid.eq.0) THEN
       madelung_final = madelung_final_sum
       col_pot_1      = col_buf_pot
       exc_vxc        = exc_vxc_buf
   ENDIF        
#endif
   
!!    Later, this subroutine is called outside of this subroutine
!!    CALL EXCHANGE_STR(lm,lmmax,exc_vxc)
   
!! Only the diagonal position of the exchange correlation stress tensor contributes to the final result
   DO index2 = 1, 9
      IF(index2 .EQ. 1 .OR. index2 .EQ. 5 .OR. index2 .EQ. 9) THEN
        exc_vxc(index2) = exc_vxc(index2) + E_xc_m_V_xc      
      ENDIF            
   ENDDO   
   

   
   WRITE(6,*) '======================= TOTAL STRESS DECOMPOSITION ==============================='
   WRITE(6,'("  i  j     coul_pot_sp      madelung     Coul_pot_int     exc_vxc       GGA_cont")')
   
   index2 = 0
   DO beta = 1, 3
      DO alpha = 1, 3
         index2 = index2 + 1
         tot_stress(alpha,beta) =  0.5*REAL(col_pot_1(index2))  + &
                                   madelung_final(index2)       + &
                                   coulomb_str_int(index2)      + &
                                   exc_vxc(index2)              - &
                                   gga_fin_str(index2)   
                                   
        WRITE(6,'(2i3,5f15.6)') alpha,beta, col_pot_1(index2)  ,     &
                                madelung_final(index2)         ,     &
                                REAL(coulomb_str_int(index2))  ,     &
                                exc_vxc(index2)                ,     &
                                gga_fin_str(index2)
      ENDDO
   ENDDO
   
   !Stress are written as Ry/vol
   tot_stress = tot_stress /vol
   
   !Only the master processor allows to WRITE in the output file
   !Print ONLY for the non-relativistic case
   
   IF(myid.eq.0) THEN
     
     IF(irel.eq.'NREL') THEN             
     WRITE(21,*)
     WRITE(21,80)    
     DO alpha = 1, 3
        WRITE(21,77) alpha, (tot_stress(alpha,beta),beta = 1, 3)
     ENDDO
     WRITE(21,*)
     ENDIF  
     
   ENDIF  

      
   
77 FORMAT(':STR_DENS',i3.3,':',2x,3f20.12)       
80 FORMAT(':STRESS (Ry/Bohr^3)')
   
RETURN
END SUBROUTINE FINL_ELECT_STR


SUBROUTINE LM_COMBINE(imax,lmmax,clm2,lm1,jatom)
!This combine the different lm component of cubic spherical harmonics
   USE struct, ONLY : nat
   
   IMPLICIT NONE
   INCLUDE 'param.inc'
   
   INTEGER, INTENT(IN)  :: imax,lmmax,lm1(2,ncom+3), jatom
   COMPLEX*16           :: clm2(1:nrad,ncom+3)
   REAL*8               :: c_kub(0:10,0:10),sq1,sqrt2,c1,c2,c3
   INTEGER              :: i,j   
   
!Coefficient for REAL spherical harmonics in cubic system
  c_kub=0.0d0
  c_kub(0,0)=1.d0
  c_kub(3,2)=1.d0
  c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
  c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
  c_kub(6,0)=.5d0*SQRT(.5d0)
  c_kub(6,2)=.25d0*SQRT(11.d0)
  c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
  c_kub(6,6)=-.25d0*SQRT(5.d0)
  c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
  c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
  c_kub(8,0)=.125d0*SQRT(33.d0)
  c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
  c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
  c_kub(9,2)=.25d0*SQRT(3.d0)
  c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
  c_kub(9,6)=-.25d0*SQRT(13.d0)
  c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
  c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
  c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
  c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
  c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
  c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
  c_kub(10,10)=-.0625d0*SQRT(85.d0)

  sqrt2 = sqrt(2.d0)
  
  !In the following different lm components of the electronic charge density combine differently
  
  i = 1
5 CONTINUE
     IF( i .gt. lmmax ) GOTO 6
     IF( ( lm1(1,i) .eq. 0 ) .and. ( lm1(2,i) .eq. 0 ) ) THEN
     
     i = i + 1
     ELSEIF ( ( lm1(1,i) .eq. -3 ) .and. ( lm1(2,i) .eq. 2 )) THEN
     
        DO j = 1, imax
           clm2(j,i) = clm2(j,i)/sqrt2
        ENDDO
        
     i = i + 1
     ELSEIF (( lm1(1,i) .eq. 4 ) .OR. ( lm1(1,i) .eq. 6 ) .OR. &
            ( lm1(1,i) .eq. -7 ) .OR. ( lm1(1,i) .eq. -9 )) THEN
            
            IF( lm1(2,i) .eq. 0 ) THEN
               sq1 = 1.d0
            ELSE
               sq1 = sqrt2
            ENDIF
            
            c1 = c_kub(IABS(lm1(1,i)),lm1(2,i) )
            c2 = c_kub(IABS(lm1(1,i)),lm1(2,i)+4 )
            
            DO J = 1, imax
               clm2(j,i)   = clm2(j,i)*c1 + clm2(j,i+1)*c2
               clm2(j,i+1) = clm2(j,i)*c2/sqrt2
               clm2(j,i)   = clm2(j,i)*c1/sq1               
            ENDDO
            
     i = i + 2
     ELSEIF( ( lm1(1,i) .eq. 8 ) .OR.( lm1(1,i) .eq. 10 ) ) THEN   
     
           IF( lm1(2,i) .eq. 0 ) THEN
              sq1 = 1.d0
           ELSE
              sq1 = sqrt2
           ENDIF
           
           c1 = c_kub(IABS(lm1(1,i)),lm1(2,i) )
           c2 = c_kub(IABS(lm1(1,i)),lm1(2,i)+4 )
           c1 = c_kub(IABS(lm1(1,i)),lm1(2,i)+8 )
           
           DO j = 1, imax
              clm2(j,i)   = clm2(j,i)*c1 + clm2(j,i+1)*c2 + clm2(j,i+2)*c3
              clm2(j,i+1) = clm2(j,i)*c2/sqrt2
              clm2(j,i+2) = clm2(j,i)*c3/sqrt2
              clm2(j,i)   = clm2(j,i)*c1/sq1              
           ENDDO
           
     ELSE 
           WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
           stop '=============incorrect LM in LM_COMBINE========= '
     ENDIF      
  GOTO 5
6 CONTINUE  
  
  
RETURN
END SUBROUTINE LM_COMBINE
           
SUBROUTINE EXCHANGE_STR(lm,lmmax,exc_vxc)
!! This subroutine calculates the exchange correlation stress tensor as given in Eq.(G7) 
!! inside atomic sphere region. In the following, different cubic l and m are handled
!! differently

   USE densit  ,   ONLY: clmsp
   USE work    ,   ONLY: vxc,exc
   USE struct  ,   ONLY: mult,jri,r0,dx,rmt,iatnr,nat
   USE parallel,   ONLY: nstart,nstop,myid
   USE norm_kub,   ONLY: c_kub
   
   IMPLICIT NONE
   INCLUDE 'param.inc'   
   
   INTEGER, INTENT(in)  :: lm(2,ncom+3,nat), lmmax(nat)
   LOGICAL              :: HEXA, CUBIC
   INTEGER              :: jatom, lm1, mu, l0, m0, ispin, J, alpha,beta,              &
                           index2
                           
   REAL*8               :: r(nrad), rho,rho1, sqrt_4pi,sqrt_4pii, pi, out1,           &
                           VALUE(nrad), exc_vxc(1:9), exc_vxc_buf1(ncom,nat),         &
                           buf_exc_vxc, exc_vxc_buf
                               
  
   pi        = 4.0d+0*ATAN(1.d0)
   sqrt_4pi  = sqrt( 4.D0*pi)
   sqrt_4pii = 1.d0/sqrt_4pi
   
   exc_vxc_buf = 0.d0   
   
   DO 20 jatom = nstart(myid), nstop(myid)
   DO 20 mu = 1, mult(jatom)
   
         HEXA  = .FALSE.
         CUBIC = .FALSE.
         IF(iatnr(jatom) .lt. 0) HEXA  = .TRUE.
         IF(iatnr(jatom) .gt. 0) CUBIC = .TRUE.
         
         DO J = 1, jri(jatom)
            r(j) = r0(jatom)*exp( dx(jatom)*(j-1) )
         ENDDO
         
         DO 40 lm1 = 1, lmmax(jatom)
            l0 = IABS( lm(1,lm1,jatom) )
            m0 =       lm(2,lm1,jatom)
            
            exc_vxc_buf1(lm1,jatom) = 0.d0
         DO 40 ispin = 1, 2
            
               IF(lm1 .eq. 1) THEN
               
                 DO j = 1, jri(jatom)
                    RHO = clmsp(j,lm1,jatom,ispin)*sqrt_4pii
                    VALUE(j) = rho*( exc(j,lm1,jatom,ispin)-vxc(j,lm1,jatom,ispin) )*0.5d0
                 ENDDO
                 
               ELSEIF(CUBIC .and. (l0 .ne. 9 .and. m0 .eq. 4)) THEN
                  GOTO 40
               ELSEIF(CUBIC.and. (m0.eq.6 .or. m0.eq.8 .or. m0.eq. 10 ) ) THEN
                  GOTO 40
               ELSEIF(CUBIC.and. (l0.eq.4 .or. l0.eq.6 .or. l0.eq.7 .or. l0.eq.8 .or.l0.eq.9 .or.l0.eq.10 ) ) THEN
                 DO j = 1, jri(jatom)
                    rho = clmsp(j,lm1,jatom,ispin)*c_kub(l0,m0) + &
                          clmsp(j,lm1+1,jatom,ispin)*c_kub(l0,m0+4)
                    IF(l0 .eq. 8 .or. l0 .eq. 10) rho = rho+ clmsp(j,lm1+2,jatom,ispin)*c_kub(l0,m0+8)     
                    rho1 = ( Exc(j,lm1,jatom,ispin)-Vxc(j,lm1,jatom,ispin) )*c_kub(l0,m0) + &
                           ( Exc(j,lm1+1,jatom,ispin)-Vxc(j,lm1+1,jatom,ispin) )*c_kub(l0,m0+4)  
                    IF(l0 .eq. 8 .or. l0 .eq. 10) rho1 = rho1 + &
                        ( Exc(j,lm1+2,jatom,ispin)-Vxc(j,lm1+2,jatom,ispin) )*c_kub(l0,m0+8)
                        
                    VALUE(j) = rho*rho1*0.5d0    
                 ENDDO
                 
                 ELSE
                 
                   DO J = 1, jri(jatom)
                      rho = clmsp(j,lm1,jatom,ispin)
                      VALUE(j) = rho*(Exc(j,lm1,jatom,ispin)-Vxc(j,lm1,jatom,ispin))*0.5d0
                   ENDDO
                   
               ENDIF
               
               call charge(r,dx(jatom),VALUE,1,jri(jatom),out1)
               
               !sum over different spin states
               exc_vxc_buf1(lm1,jatom) = exc_vxc_buf1(lm1,jatom) + out1  !sum over all spin states              
    40 CONTINUE                             
    20 CONTINUE                                  
       
       DO 50 jatom = nstart(myid),nstop(myid)
          DO 50 mu = 1, mult(jatom)
             HEXA  = .FALSE.
             CUBIC = .FALSE.
             IF(iatnr(jatom) .lt. 0) HEXA  = .TRUE.
             IF(iatnr(jatom) .gt. 0) CUBIC = .TRUE.
             buf_exc_vxc = 0.d0
             DO 55 lm1 = 1, lmmax(jatom)
                l0 = IABS( lm(1,lm1,jatom) )
                m0 =       lm(2,lm1,jatom)                
                IF(CUBIC) THEN
                   IF( m0.eq.4 .or. m0.eq.6 .or. m0 .eq.8 .or. m0 .eq. 10 ) GOTO 55
                ENDIF
                buf_exc_vxc = buf_exc_vxc + exc_vxc_buf1(lm1,jatom)     !sum over all lm states
       55 CONTINUE     
             exc_vxc_buf = exc_vxc_buf + buf_exc_vxc                    !sum over equivalent and non-equivalent atoms
  50   CONTINUE     
       
     
       index2 = 0
       DO alpha = 1, 3
          DO beta = 1, 3
             index2 = index2 + 1
             IF(alpha .eq. beta) exc_vxc(index2) = exc_vxc_buf*2.d0                           
          ENDDO
       ENDDO
       
RETURN
END SUBROUTINE EXCHANGE_STR


SUBROUTINE SYMMETRY_STRESS_ROTIJ(index_a,jatom,ele_tensor)

  USE struct,  ONLY : rotij, rotloc, lattic, ndif

  IMPLICIT NONE
    INTEGER, INTENT(IN)    :: index_a,jatom
    REAL*8 , INTENT(INOUT) :: ele_tensor(1:9) 
    
    REAL*8                 :: buf_prod1(1:3,1:3), buf_prod2(1:3,1:3),             &
                              detinv,B1(3,3), A1(3,3), buf_prod3(1:3,1:3),        &
                              final_stress(3,3), rotij_st(3,3,ndif), BR4(3,3)
                                                             
    INTEGER                :: alpha, beta, index2, ii
    
    rotij_st = rotij        
    
    IF(LATTIC(1:1) .EQ. 'H') THEN
      BR4(1,1)=SQRT(3.d0)/2.d0                                              
      BR4(1,2)=-.5d0     
      BR4(1,3)=0.0d0                                                      
      BR4(2,1)=0.0d0                                                     
      BR4(2,2)=1.0d0                                                      
      BR4(2,3)=0.0d0                                                      
      BR4(3,1)=0.0d0                                                      
      BR4(3,2)=0.0d0                                                      
      BR4(3,3)=1.d0                                   

    ELSEIF(LATTIC(1:1) .EQ. 'R') THEN
      BR4(1,1)=1/2.d0/sqrt(3.d0)
      BR4(1,2)=-1/2.d0                                                     
      BR4(1,3)=1/3.d0                                                      
      BR4(2,1)=1/2.d0/SQRT(3.d0)                                          
      BR4(2,2)=1*0.5d0                                                
      BR4(2,3)=1/3.d0                                                      
      BR4(3,1)=-1/SQRT(3.d0)                                         
      BR4(3,2)=0.d0                                                
      BR4(3,3)=1/3.d0            
    ENDIF


    !This if condition only applicable for Hexagonal and Rhombohedral lattice
    IF( (LATTIC(1:1) .eq. 'H' ) .OR. (LATTIC(1:1) .eq. 'R' ) ) THEN   
    
      A1(:,:) = BR4(:,:)
      
      !Calculate the determinant of 3x3 matrix
      detinv = 1/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
              
      IF(detinv .eq. 0.0) STOP 'inverse of sym_struct does not exist'    
      
      !IF the determinant is not zero, compute the inverse of the matrix
      B1(1,1) = +detinv * (A1(2,2)*A1(3,3) - A1(2,3)*A1(3,2))
      B1(2,1) = -detinv * (A1(2,1)*A1(3,3) - A1(2,3)*A1(3,1))
      B1(3,1) = +detinv * (A1(2,1)*A1(3,2) - A1(2,2)*A1(3,1))
      B1(1,2) = -detinv * (A1(1,2)*A1(3,3) - A1(1,3)*A1(3,2))
      B1(2,2) = +detinv * (A1(1,1)*A1(3,3) - A1(1,3)*A1(3,1))
      B1(3,2) = -detinv * (A1(1,1)*A1(3,2) - A1(1,2)*A1(3,1))
      B1(1,3) = +detinv * (A1(1,2)*A1(2,3) - A1(1,3)*A1(2,2))
      B1(2,3) = -detinv * (A1(1,1)*A1(2,3) - A1(1,3)*A1(2,1))
      B1(3,3) = +detinv * (A1(1,1)*A1(2,2) - A1(1,2)*A1(2,1))
   
      buf_prod2(:,:) = matmul(B1(:,:),rotij(:,:,index_a))
      buf_prod1(:,:) = matmul(buf_prod2(:,:),A1(:,:))
      
      rotij_st(:,:,index_a) = buf_prod1(:,:)

    ENDIF
   
    index2 = 0
    DO alpha = 1, 3
       DO beta = 1,3
          index2 = index2 + 1
          buf_prod1(alpha,beta) = ele_tensor(index2)
       ENDDO
    ENDDO
    
    buf_prod3 = buf_prod1
    
    !! This is for rotloc matrix that converts the coordinate from local to
    !! the global coordinate system
    A1(:,:) = rotloc(:,:,jatom)
    detinv = 1/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
              
    IF(detinv .eq. 0.0) STOP 'INVERSE OF ROTLOC DOES NOT EXIST -- EC stress'    
    
    B1(1,1) = +detinv * (A1(2,2)*A1(3,3) - A1(2,3)*A1(3,2))
    B1(2,1) = -detinv * (A1(2,1)*A1(3,3) - A1(2,3)*A1(3,1))
    B1(3,1) = +detinv * (A1(2,1)*A1(3,2) - A1(2,2)*A1(3,1))
    B1(1,2) = -detinv * (A1(1,2)*A1(3,3) - A1(1,3)*A1(3,2))
    B1(2,2) = +detinv * (A1(1,1)*A1(3,3) - A1(1,3)*A1(3,1))
    B1(3,2) = -detinv * (A1(1,1)*A1(3,2) - A1(1,2)*A1(3,1))
    B1(1,3) = +detinv * (A1(1,2)*A1(2,3) - A1(1,3)*A1(2,2))
    B1(2,3) = -detinv * (A1(1,1)*A1(2,3) - A1(1,3)*A1(2,1))
    B1(3,3) = +detinv * (A1(1,1)*A1(2,2) - A1(1,2)*A1(2,1))
    
    buf_prod2(:,:) = matmul( B1(:,:),buf_prod1(:,:) )
    buf_prod1(:,:) = matmul( buf_prod2(:,:),A1(:,:) )
    
    !! This is for rotij matrix that converts the coordinate system of one atom 
    !! into another atoms.
    A1(:,:) = rotij_st(:,:,index_a)
        
    detinv = 1/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
              
    IF(detinv .eq. 0.0) STOP 'INVERSE OF ROTIJ DOES NOT EXIST -- -- EC stress'
    
    B1(1,1) = +detinv * (A1(2,2)*A1(3,3) - A1(2,3)*A1(3,2))
    B1(2,1) = -detinv * (A1(2,1)*A1(3,3) - A1(2,3)*A1(3,1))
    B1(3,1) = +detinv * (A1(2,1)*A1(3,2) - A1(2,2)*A1(3,1))
    B1(1,2) = -detinv * (A1(1,2)*A1(3,3) - A1(1,3)*A1(3,2))
    B1(2,2) = +detinv * (A1(1,1)*A1(3,3) - A1(1,3)*A1(3,1))
    B1(3,2) = -detinv * (A1(1,1)*A1(3,2) - A1(1,2)*A1(3,1))
    B1(1,3) = +detinv * (A1(1,2)*A1(2,3) - A1(1,3)*A1(2,2))
    B1(2,3) = -detinv * (A1(1,1)*A1(2,3) - A1(1,3)*A1(2,1))
    B1(3,3) = +detinv * (A1(1,1)*A1(2,2) - A1(1,2)*A1(2,1))    
    
    buf_prod2(:,:) = matmul(B1(:,:),buf_prod1(:,:))
    buf_prod1(:,:) = matmul(buf_prod2(:,:),A1(:,:))
    

   index2 = 0
   DO alpha = 1, 3
      DO beta = 1,3
         index2 = index2+1
         ele_tensor(index2) = buf_prod1(alpha,beta)
      ENDDO
   ENDDO

   
RETURN
END SUBROUTINE SYMMETRY_STRESS_ROTIJ

