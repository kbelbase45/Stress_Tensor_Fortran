! This subroutine calculates the trace of core correction stress tensor.  
! It is called each time from within the loop through a different atom of the system.

  SUBROUTINE core_corr(natm, nrp, r1, rhoc, jatom, index_a) 
 !natm different atoms in the unit cell
 !nrp  total radial points upto the atomic boundary
 !r1   radial points on the radial grid (1 to nrp)
 !rhoc r*r*spherical core charge density 
 !dv1p spherical potential (L=0,M=0)
       USE TestCases
       USE struct,  only: mult, dx
       USE core_sp, only: dv1p
       IMPLICIT REAL*8 (A-H,O-Z)
       INCLUDE 'param.inc'             
       !NRAD specified in param.inc, which is 881       
                  
       INTEGER, intent(in)    :: natm, nrp, jatom, index_a
       REAL*8,  intent(in)    :: r1(nrad+141), rhoc(nrad+141)
       INTEGER                :: ir,jatm,imu, counter, ii
       REAL*8                 :: dx1, dpot(nrad+141),value(nrad+141),V1(NRAD+141)
       REAL*8                 :: SQ4pi,one,two,TC,Pi,value2(nrad+141)       
       REAL*8                 :: sum1,sum_p_last,tot_char,int_test                                       
       REAL*8                 :: tot_charge_781,tot_charge_922
        
       PI    = 4.0D+0*ATAN(1.0D+0)             
       sq4pi = 1.d0/sqrt(4.d0*pi) 
       one   = 1.d0
       two   = 2.d0
       !increment of the logarithmic mesh of different atoms
       dx1 = dx(jatom)
        
       sum1       = 0.d0
       sum_p_last = 0.d0           
       
       ! The radial mesh used in the core calculation program has 781 radial points. 
       ! This is 141 additional point then specified in structure file, nrp = 781.
       
       counter = nrp + 141
          ! imu is loop over the same type of atoms in the unit cell
          ! for silicon, jatom = 1 and MULT(jatom) = 2.
          DO imu  = 1, MULT(jatom)
               !loop over the radial point inside an atomic sphere
               DO ir = 1, counter                      
                     V1(ir) = two*DV1P(ir)             ! H to Ry      
               ENDDO               
               
               !calculate the derivate of the potential
               call dergl(counter,r1,v1,dpot,.false.,g0,.false.)          
               
               !calculate r^3*dv and r^2*v on each radial points 
               DO ir = 1, counter                          
                  value(ir) = rhoc(ir)*dpot(ir)*r1(ir)
                  value2(ir)= rhoc(ir)*V1(ir)
               enddo
               
               !integrate int_0^R(r^3*dv); R = r1|_(ir=781)
               !R atomic sphere size
               call chargel2(r1,1.d0,one,value(1),dx1,nrp,TC) 
               
               !sum for all equivalent and non-equivalent atoms
               sum1 = sum1 + TC               
               
               !Integration int_0^R(r^3*dv); R = r1|_(ir=781+141)
               !R larger than atomic sphere size
               call chargel2(r1,1.d0,one,value(1),dx1,counter,TC)                     
               sum_p_last = sum_p_last + TC               
               
               !Integration int_0^R(r^2*core_charge_density)
               !compute total core charge in the whole region 
               call chargel2(r1,1.d0,one,rhoc,dx1,counter,tot_char)                              
               
               !Integration int_0^R(r^2*v); R = r1|_(ir=781+141)
               call chargel2(r1,1.d0,one,value2(1),dx1,counter,int_test)
                                         
          ENDDO
    
       
       !Calculate all nine different components of the core correction stress tensor
       call core_corr_stress(natm,nrp,r1,rhoc,dpot,v1,jatom, index_a) 
             
    RETURN
    END  
    
    SUBROUTINE core_corr_stress(natm,nrp,r1,rhoc,dpot,v1,jatom, index_a)
    ! This is the subroutine that computes Eq (6.48) of Core_Correction_stress.pdf.
    ! Structural parameters are accessed through struct module.
    ! The different l,m components (lm11_st) of the non-spherical potential are 
    ! supplied through the module core_nsp. This is not useful for the total energy
    ! calculation because only the spherically average (l=0,m=0) potential contributes.
      USE TestCases
      USE struct  , only : nat, iatnr, mult      
      USE core_nsp, only : lm11_st, lmmax_st, vns_st, sum_tot
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'param.inc'
      
      INTEGER, INTENT(IN) :: natm,nrp,jatom, index_a
      
      REAL*8,  INTENT(IN) :: dpot(nrad+141),r1(nrad+141),rhoc(nrad+141),V1(nrad+141)                                                                
                             
      REAL*8              :: one,two, pi, dx,  TC,TC1, sq4pi, sqrt2, sq4_pi 
                                  
      REAL*8              :: cor_corr_sph(1:9),cor_corr_ns1(1:9),cor_corr_ns2(1:9),              &
                             cor_corr_tot(1:9),cor_corr_ns1_buf(1:9),cor_corr_ns2_buf(1:9)
                             
                             
      INTEGER             :: t,tp,s, alpha,beta,index2, lm1p,lm1, lmmax_p, ll_p, mm_p,            &
                             counter,ri,insv, lmmax2,lmmax3, lmtot(2,ncom+3), lmtot1
                             
                             
      COMPLEX*16          :: cabt(3,-1:1),ca_st_lm(3,-1:1,-1:1), zeroc, zero, imag1, imag2,        &
                             fac_st(ncom+3,nat),fc(ncom+3,nat),value_c(nrad+141),buf_sum(1:9),     &   
                             out_c,int_ns1(ncom+3),int_ns2(ncom+3),vtmp_cmplx(nrad+141,ncom+3,nat)                               
      
      REAL*8              :: dvtmp(nrad+141),g0, gaunt1, int_out, vtmp(nrad+141), value1(nrad+141)
            
      INTEGER             :: lm(2,ncom+3), lmmax22, lm11(2,ncom+3),jatm,mu,  ir
      
      REAL*8              :: val_real(nrad+141),val_cmplx(nrad+141),value(nrad+141), &
                             dval_real(nrad+141),dval_cmplx(nrad+141), out1, out2 
      
      ! The following two functions are used to compute the gaunt numbers (G in .pdf).                        
      REAL*8,  external   :: GAUNT
      INTEGER, external   :: NOTRI   
      
     !=========here what we have===========
     !v1       = spherical potential
     !dpot     = derivative of v1 
     !rhoc     = spherical density          
     !cabt     = c-matrix coefficients associated with the unit vector expansion with respect to spherical harmonics
     !ca_st_lm = c-matrix coefficients associated with the derivative of the spherical harmonics
            
      one     = 1.d0
      two     = 2.d0
      PI      = 4.0D+0*ATAN(1.0D+0) 
      dx      = log(r1(nrp)/r1(1))/(nrp-1)     
      sq4pi   = 1.d0/sqrt(4.d0*pi)
      sq4_pi  = sqrt(4.d0*pi)
      zeroc   = (0.d0,0.d0) 
      zero    = 0.d0
      sqrt2   = sqrt(2.D0)
      imag2   = (0.d0,1.d0)
      counter = nrp+141
      
      call c_alpha_m(cabt)      !c-matrix of unit vector expansion
                           
      cor_corr_sph  = zero!zeroc      !spherical part of the core correction stress tensor
      cor_corr_ns1  = zero!zeroc      !non-spherical part1         
      cor_corr_ns2  = zero!zeroc      !non-spherical part2                    
      cor_corr_ns1_buf = zero!zeroc
      cor_corr_ns2_buf = zero!zeroc            
      
      lmmax_p = lmmax_st(jatom)        
!============reading of non-spherical potential==========================      
      DO lm1p = 1, lmmax_p
         lm11(1,lm1p) = lm11_st(1,lm1p,jatom)
         lm11(2,lm1p) = lm11_st(2,lm1p,jatom)                  
         
         DO ri = 1, nrp
             vtmp_cmplx(ri,lm1p,jatom) = vns_st(ri,lm1p,jatom)
             val_real(ri) = vtmp_cmplx(ri,lm1p,jatom)
         ENDDO
                          
         call dergl(counter,r1,val_real,dval_real,.false.,g0,.false.)         
         
         DO ri = 1, counter
           value(ri) = r1(ri)*rhoc(ri)*dval_real(ri) 
         ENDDO          
         
         call chargel2(r1,1.d0,one,value(1),dx,counter,out1)
         int_ns1(lm1p) = out1
           
         DO ri = 1, counter  
           value(ri) = rhoc(ri)*val_real(ri)
         ENDDO           
         call chargel2(r1,1.d0,one,value(1),dx,counter,out2)    
         int_ns2(lm1p) = out2
         
      ENDDO      
!================end reading non-spherical potential===================        

!===================New modification====================2021     
        lmmax22  = lmmax_p
        DO lm1p = 1, lmmax_p           
           mm_p = lm11(2,lm1p)
           IF(mm_p .ne. 0) then
             lmmax22 = lmmax22 + 1 
              lm11(1,lmmax22) =  lm11(1,lm1p)
              lm11(2,lmmax22) = -lm11(2,lm1p)              
           endif
        ENDDO              
!                
      IF(iatnr(jatom) .lt. 0) THEN         
         lmmax_new  = lmmax_p
         DO lm1 = 1, lmmax_p            
            mm_p = lm11(2,lm1)
            if(mm_p .ne. 0) then
              lmmax_new = lmmax_new + 1
              if(lmmax_new .gt. ncom+3) STOP 'error in "core_corr.f" in lm list'
              vns_st(1:nrp,lmmax_new,jatom) = vns_st(1:nrp,lm1,jatom)
              vtmp_cmplx(1:counter,lmmax_new,jatom) =  vtmp_cmplx(1:counter,lm1,jatom)              
              int_ns1(lmmax_new) = int_ns1(lm1)
              int_ns2(lmmax_new) = int_ns2(lm1)                                                       
            endif                                    
         ENDDO
                  
         call multfc_core( fc,jatom,lmmax_new,lm11(1,1) )            
                  
         DO lm1 = 1, lmmax_new                      
!             vns_st(1:nrp,lm1,jatom) = vns_st(1:nrp,lm1,jatom)*fc(lm1,jatom)
            vtmp_cmplx(1:counter,lm1,jatom) = vtmp_cmplx(1:counter,lm1,jatom)*fc(lm1,jatom)
                int_ns1(lm1) = int_ns1(lm1)*fc(lm1,jatom)
                int_ns2(lm1) = int_ns2(lm1)*fc(lm1,jatom)    
         ENDDO               
      !end recently added lines
      
         lmtot1 = 0
         DO 39 lm1 = 1, lmmax_new
            lmtot1 = lmtot1 + 1
            lmtot(1,lmtot1) = lm11(1,lm1)
            lmtot(2,lmtot1) = lm11(2,lm1)
            
             vns_st(1:nrp,lmtot1,jatom) = vns_st(1:nrp,lm1,jatom)
              vtmp_cmplx(1:counter,lmtot1,jatom) = vtmp_cmplx(1:counter,lm1,jatom)
              int_ns1(lmtot1) = int_ns1(lm1)
              int_ns2(lmtot1) = int_ns2(lm1)                                            
            
            IF(lmtot1 .eq. 1) GOTO 39
            IF( ( IABS( lmtot(1,lmtot1)) .eq. IABS( lmtot(1,lmtot1-1))  ) .AND.   &
              (  lmtot(2,lmtot1) .eq.  lmtot(2,lmtot1-1)  )  ) THEN
!               vns_st(1:nrp,lmtot1-1,jatom) = vns_st(1:nrp,lmtot1-1,jatom) + vns_st(1:nrp,lmtot1,jatom)
              vtmp_cmplx(1:counter,lmtot1-1,jatom) = vtmp_cmplx(1:counter,lmtot1-1,jatom) + vtmp_cmplx(1:counter,lmtot1,jatom)
                  int_ns1(lmtot1-1) = int_ns1(lmtot1-1) + int_ns1(lmtot1)
                  int_ns2(lmtot1-1) = int_ns2(lmtot1-1) + int_ns2(lmtot1)                  
              lmtot1 = lmtot1 - 1
            ENDIF                                      
39       CONTINUE         

      ELSE
!          call lm_combine(counter,lmmax_p,vns_st(1,1,jatom),lm11(1,1),jatom)
         call lm_combine(counter,lmmax_p,vtmp_cmplx(1,1,jatom),lm11(1,1),jatom)
            call multsu(int_ns1(1),lmmax_p,lm11(1,1))
            call multsu(int_ns2(1),lmmax_p,lm11(1,1))
         
         lmtot1 = lmmax_p
         DO lm1 = 1, lmmax_p
            lmtot(1,lm1) = lm11(1,lm1)
            lmtot(2,lm1) = lm11(2,lm1)
            IF( lm11(2,lm1) .ne. 0 ) THEN
               lmtot1 = lmtot1 + 1
               IF(lmtot1 .gt. ncom + 3) STOP 'error in subroutine "core_corr_stress" '
               lmtot(1,lmtot1) =  lm11(1,lm1)
               lmtot(2,lmtot1) = -lm11(2,lm1)
               
!                vns_st(1:nrp,lmtot1,jatom) = vns_st(1:nrp,lm1,jatom)
               vtmp_cmplx(1:counter,lmtot1,jatom) = vtmp_cmplx(1:counter,lm1,jatom)
               int_ns1(lmtot1) = int_ns1(lm1)
               int_ns2(lmtot1) = int_ns2(lm1)
            ENDIF
         ENDDO
         
      ENDIF
!===================end new modification================         
      
   DO mu = 1, MULT(jatom)               !loop over the same type of atoms in the unit cell
      index_a = index_a + 1
      
            DO ir = 1, counter
               value1(ir) = r1(ir)*rhoc(ir)*dpot(ir)               
            ENDDO
            call chargel2(r1,1.d0,one,value1(1),dx,counter,TC)           
            
            DO t = -1,1
               DO tp = -1,1
                  IF(-t .ne. tp) CYCLE
                  index2 = 0
                  DO alpha = 1, 3
                     DO beta = 1, 3
                        index2 = index2 + 1
                        cor_corr_sph(index2) = cor_corr_sph(index2) + &
                        (-1)**t*cabt(alpha,t)*cabt(beta,tp)*TC                                                
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO                                                                                         
            
      !first part of the non-spherical core correction contribution      
      cor_corr_ns1_buf = 0.0
      buf_sum  = (0.0d0, 0.d0)
      DO lm1p = 1, lmtot1
         ll_p = lmtot(1,lm1p)
         mm_p = lmtot(2,lm1p)                                                       
                  
         DO ri = 1, counter
            val_real(ri)  = REAL( vtmp_cmplx(ri,lm1p,jatom) )  
            val_cmplx(ri) = AIMAG( vtmp_cmplx(ri,lm1p,jatom) )
            
!             val_real(ri)  = REAL( vtmp_cmplx(ri,lm1p,1) )  
!             val_cmplx(ri) = AIMAG( vtmp_cmplx(ri,lm1p,1) )            
            
!             val_real(ri)  = REAL( vns_st(ri,lm1p,jatom) )  
!             val_cmplx(ri) = AIMAG( vns_st(ri,lm1p,jatom) )            
         ENDDO              
         
         call dergl(counter,r1,val_real,dval_real,.false.,g0,.false.)
         call dergl(counter,r1,val_cmplx,dval_cmplx,.false.,g0,.false.)
         
         DO ri = 1, counter
              val_real(ri)  = r1(ri)*rhoc(ri)*dval_real(ri)
              val_cmplx(ri) = r1(ri)*rhoc(ri)*dval_cmplx(ri)
         ENDDO

         call chargel2(r1,1.d0,one,val_real(1),dx,counter,out1)
         call chargel2(r1,1.d0,one,val_cmplx(1),dx,counter,out2)                  
         
         out_c = dcmplx(out1,out2)                      
         
         IF( NOTRI(ll_p,1,1) .lt. 0) CYCLE           
           DO t = -1,1
              DO tp = -1,1
                 IF( (mm_p+t+tp) .ne. 0 ) CYCLE                      
                  gaunt1 = GAUNT(ll_p,1,1,-mm_p,t,tp)                                       
                                      
                  index2 = 0
                  DO alpha = 1, 3
                     DO beta = 1, 3
                        index2 = index2 + 1                         
                        
                        cor_corr_ns1_buf(index2) = cor_corr_ns1_buf(index2) + &
                        cabt(alpha,t)*cabt(beta,tp)*gaunt1*(-1)**mm_p*out_c                                                    
                        
!                         cor_corr_ns1_buf(index2) = cor_corr_ns1_buf(index2) + &
!                         cabt(alpha,t)*cabt(beta,tp)*gaunt1*(-1)**mm_p*int_ns1(lm1p)                          
                          
                     ENDDO
                  ENDDO                                         
              ENDDO
           ENDDO           
      ENDDO                            
      
      call symmetry_stress_rotij(index_a,jatom,cor_corr_ns1_buf)
      cor_corr_ns1 = cor_corr_ns1 + cor_corr_ns1_buf                                    
      
      !second part of the non-spherical core correction contribution
      cor_corr_ns2_buf = 0.0
      DO lm1p = 1, lmtot1
         ll_p = lmtot(1,lm1p)
         mm_p = lmtot(2,lm1p)
         
         DO ri = 1, counter
            value_c(ri) = rhoc(ri)*vtmp_cmplx(ri,lm1p,jatom) 
            
!             value_c(ri) = rhoc(ri)*vtmp_cmplx(ri,lm1p,1)             
!             value_c(ri) = rhoc(ri)*vns_st(ri,lm1p,jatom) 
         ENDDO
         val_real  = REAL( value_c )
         val_cmplx = AIMAG( value_c )
         call chargel2(r1,1.d0,one,val_real(1),dx,counter,out1)
         call chargel2(r1,1.d0,one,val_cmplx(1),dx,counter,out2)
         
         out_c = dcmplx(out1,out2)                  
         
         DO s = -1,1,2
            IF( (ll_p+s) .ne. 1) CYCLE
            DO t = -1,1
               DO tp = -1,1
                  IF(-tp .ne. (mm_p+t) ) CYCLE
                  CALL C_FACTOR(ll_p,mm_p,ca_st_lm)                     
                   index2 = 0
                   DO alpha = 1, 3
                      DO beta = 1,3
                         index2 = index2 + 1                                                   
                         cor_corr_ns2_buf(index2) = cor_corr_ns2_buf(index2) + &
                         0.5d0*( ca_st_lm(alpha,s,t)*cabt(beta,tp) + &
                           ca_st_lm(beta,s,t)*cabt(alpha,tp) )*(-1)**tp*out_c
                     !or      
!                          cor_corr_ns2_buf(index2) = cor_corr_ns2_buf(index2) + &
!                          0.5d0*( ca_st_lm(alpha,s,t)*cabt(beta,tp) + &
!                            ca_st_lm(beta,s,t)*cabt(alpha,tp) )*(-1)**tp*real(int_ns2(lm1p)  )                                                     
                      ENDDO
                   ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO             
      
      call symmetry_stress_rotij(index_a,jatom,cor_corr_ns2_buf)
      cor_corr_ns2 = cor_corr_ns2 + cor_corr_ns2_buf            
  ENDDO   !loop over multiplicity  
     
     !===================look here==================
      DO index2 = 1, 9                    
         cor_corr_tot(index2) = ( cor_corr_sph(index2)*sq4pi*sq4pi + cor_corr_ns1(index2)*sq4pi + &
                                cor_corr_ns2(index2)*sq4pi )!*MULT(jatom)                                           
         sum_tot(index2,jatom) = sum_tot(index2,jatom) + REAL(cor_corr_sph(index2))*sq4pi*sq4pi + REAL(cor_corr_ns1(index2))*sq4pi &
                                    + REAL(cor_corr_ns2(index2))*sq4pi*sq4pi*2.d0                                                            
         WRITE(21,77) index2, REAL(cor_corr_sph(index2))*sq4pi*sq4pi,  REAL(cor_corr_ns1(index2))*sq4pi*sq4pi, &
                                    REAL(cor_corr_ns2(index2))*sq4pi*sq4pi         
      ENDDO
     !===============================================         
      
77 FORMAT(':COR__CORR__',i3.3,':',2x,3ES21.12)   
    RETURN
    END SUBROUTINE core_corr_stress

!=======factors required to combine -l and +l
  SUBROUTINE multfc_core(fc,jatom,lmmax3,lm)
    use struct, only : nat
    implicit none
    include 'param.inc'
    integer, intent(in)   :: lmmax3, jatom
    integer, intent(in)   :: lm(2,ncom+3)
    integer               :: lm1p,mm,minu
    complex*16            :: fc(ncom+3,nat), imag1, imag
    
    imag  = (0.d0,1.d0)
    DO 1 lm1p = 1, lmmax3
         mm = lm(2,lm1p)
         IF( mm .eq. 0) THEN
            fc(lm1p,jatom) = (1.d0,0.d0)
            goto 1
         ENDIF
       imag1 = (1.d0,0.d0)
       minu  = 1
         IF(lm(1,lm1p) .lt. 0) THEN
           imag1 = -imag
           minu  = -1
         ENDIF
         IF(MOD(mm,2) .eq. 1) then
           imag1 = -imag1
           minu  = -minu
         endif
         IF(mm .gt. 0) minu = 1
         fc(lm1p,jatom) = imag1/sqrt(2.d0)*minu
1   CONTINUE    
  RETURN
  END SUBROUTINE multfc_core  
  
subroutine symmetry_stress_rotij(index_a,jatom,inout_tensor)
!This subroutine symmetrizes the core correction stress tensor using rotational matrices of the system
!Both the rotij and rotloc are read from the case.struct file
use struct, only : rotij, rotloc, lattic, ndif
implicit none
    include 'param.inc'
    integer, intent(in)    :: index_a,jatom
    real*8 , intent(inout) :: inout_tensor(1:9) 
    real*8              :: buf_prod1(1:3,1:3),         &
                           buf_prod2(1:3,1:3), A1(3,3),&
                           buf_prod3(1:3,1:3), &
                           detinv,B1(3,3), final_stress(3,3), &
                           rotij_st(3,3,ndif), BR4(3,3)
    integer             :: alpha, beta, index2, ii
    
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



    IF( (LATTIC(1:1) .eq. 'H' ) .OR. (LATTIC(1:1) .eq. 'R' ) ) THEN   
      A1(:,:) = BR4(:,:)
      detinv = 1.0/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
      IF(detinv .eq. 0.0) STOP 'inverse of sym_struct does not exist'    
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
          buf_prod1(alpha,beta) = inout_tensor(index2)
       ENDDO
    ENDDO    
    
    
!       buf_prod2(:,:) = matmul( transpose(rotloc(:,:,jatom)),buf_prod1(:,:) )
!       buf_prod1(:,:) = matmul( buf_prod2(:,:),rotloc(:,:,jatom) )
! 
!       buf_prod2(:,:) = matmul(transpose(rotij_st(:,:,index_a)),buf_prod1(:,:))
!       buf_prod1(:,:) = matmul(buf_prod2(:,:),rotij_st(:,:,index_a))
     
! !-------Below modification is done because of the problem in Rutile 
! A^inv Tensor A gives the same result as A^T Tensor A

      A1(:,:) =  rotloc(:,:,jatom) 
      
      detinv = 1.0/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
      IF(detinv .eq. 0.0) STOP 'inverse of rotloc does not exist -- core stress'    
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
        
      A1(:,:) = rotij_st(:,:,index_a)           
      detinv = 1.0/(A1(1,1)*A1(2,2)*A1(3,3) - A1(1,1)*A1(2,3)*A1(3,2)&
              - A1(1,2)*A1(2,1)*A1(3,3) + A1(1,2)*A1(2,3)*A1(3,1)&
              + A1(1,3)*A1(2,1)*A1(3,2) - A1(1,3)*A1(2,2)*A1(3,1))
      IF(detinv .eq. 0.0) STOP 'inverse of rotij does not exist -- -- core stress'    
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
         inout_tensor(index2) = buf_prod1(alpha,beta)
      ENDDO
   ENDDO   
   
RETURN
end subroutine symmetry_stress_rotij
  
subroutine lm_combine(imax,lmmax,clm2,lm1,jatom)   
   use norm_kub, only: c_kub
   implicit none
   include 'param.inc'
   
   integer, intent(in)  :: imax,lmmax,lm1(2,ncom+3), jatom
   complex*16           :: clm2(nrad+141,ncom+3)!,vpot(1:nrad,ncom+3)
!    real*8               :: c_kub(0:10,0:10)
   integer              :: i,j
   real*8               :: sq1,sqrt2,c1,c2,c3
!=======coefficient for real spherical harmonics in cubic system========   
!   c_kub=0.0d0
!   c_kub(0,0)=1.d0
!   c_kub(3,2)=1.d0
!   c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
!   c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
!   c_kub(6,0)=.5d0*SQRT(.5d0)
!   c_kub(6,2)=.25d0*SQRT(11.d0)
!   c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
!   c_kub(6,6)=-.25d0*SQRT(5.d0)
!   c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
!   c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
!   c_kub(8,0)=.125d0*SQRT(33.d0)
!   c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
!   c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
!   c_kub(9,2)=.25d0*SQRT(3.d0)
!   c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
!   c_kub(9,6)=-.25d0*SQRT(13.d0)
!   c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
!   c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
!   c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
!   c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
!   c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
!   c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
!   c_kub(10,10)=-.0625d0*SQRT(85.d0)
!===========================================================================   
  sqrt2 = sqrt(2.d0)
  
  i = 1
5 CONTINUE
     IF( i .gt. lmmax ) GOTO 6
     IF( ( lm1(1,i) .eq. 0 ) .and. ( lm1(2,i) .eq. 0 ) ) THEN
        i = i + 1
     ELSEIF ( ( lm1(1,i) .eq. -3 ) .and. ( lm1(2,i) .eq. 2 )) THEN
        DO j = 1, imax
           clm2(j,i) = clm2(j,i)/sqrt2
!            vpot(j,i) = vpot(j,i)/sqrt2
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
               
!                vpot(j,i)   = vpot(j,i)*c1 + vpot(j,i+1)*c2
!                vpot(j,i+1) = vpot(j,i)*c2/sqrt2
!                vpot(j,i)   = vpot(j,i)*c1/sq1
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
              
!               vpot(j,i)   = vpot(j,i)*c1 + vpot(j,i+1)*c2 + vpot(j,i+2)*c3
!               vpot(j,i+1) = vpot(j,i)*c2/sqrt2
!               vpot(j,i+2) = vpot(j,i)*c3/sqrt2
!               vpot(j,i)   = vpot(j,i)*c1/sq1
           ENDDO
     ELSE 
           WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
           stop '=============incorrect LM in lm_combine========= '
     ENDIF      
  GOTO 5
6 CONTINUE  
  
  
RETURN
end subroutine lm_combine



SUBROUTINE multsu(qq,lmmax,lm)
!                                                                       
!                                
!     

!   USE struct, ONLY: nat
!   use parallel, only: abort_parallel
  use norm_kub, only: c_kub
  IMPLICIT NONE
  INCLUDE 'param.inc'
  !
  COMPLEX*16   :: qq(ncom+3),imag,a,b,c, ak1,bk1,ck1
  INTEGER      :: lm(2,ncom+3),lmmax,llmm,i
  REAL*8       :: sqrt2,sq1,c1,c2,c3
!   REAL*8       :: c_kub(0:10,0:10)  
!----------------------------------------------------------------------    
!
!   c_kub=0.0d0
!   c_kub(0,0)=1.d0
!   c_kub(3,2)=1.d0
!   c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
!   c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
!   c_kub(6,0)=.5d0*SQRT(.5d0)
!   c_kub(6,2)=.25d0*SQRT(11.d0)
!   c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
!   c_kub(6,6)=-.25d0*SQRT(5.d0)
!   c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
!   c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
!   c_kub(8,0)=.125d0*SQRT(33.d0)
!   c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
!   c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
!   c_kub(9,2)=.25d0*SQRT(3.d0)
!   c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
!   c_kub(9,6)=-.25d0*SQRT(13.d0)
!   c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
!   c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
!   c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
!   c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
!   c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
!   c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
!   c_kub(10,10)=-.0625d0*SQRT(85.d0)

  sqrt2=SQRT(2.d0)
  !     
     !        CUBIC CASE
        
     i=1
     DO         
        IMAG=CMPLX(1.D0,0.D0)
        IF(lm(1,i).LT.0) IMAG=CMPLX(0.D0,-1.D0)
        
        IF(i.gt.lmmax) EXIT
        IF(lm(1,i).EQ.0.AND.lm(2,i).EQ.0) THEN
           i=i+1
        ELSEIF (lm(1,i).EQ.-3.AND.lm(2,i).EQ.2) THEN  
           qq(i)=qq(i)*imag/sqrt2                      
           i=i+1
        ELSEIF (lm(1,i).EQ.4.OR.lm(1,i).EQ.6.OR. &
             lm(1,i).EQ.-7.OR.lm(1,i).EQ.-9) THEN
           IF (lm(2,i).EQ.0) THEN
              sq1=1.d0
           ELSE
              sq1=sqrt2
           ENDIF
           c1=c_kub(ABS(lm(1,i)),lm(2,i))
           c2=c_kub(ABS(lm(1,i)),lm(2,i)+4)
           a=qq(i)*imag
           b=qq(i+1)*imag
           qq(i)=a*c1*c1/sq1 + b*c1*c2/sq1
           qq(i+1)=a*c1*c2/sqrt2 + b*c2*c2/sqrt2                       
           i=i+2
        ELSEIF (lm(1,i).EQ.8.OR.lm(1,i).EQ.10) THEN
           IF (lm(2,i).EQ.0) THEN
              sq1=1.d0
           ELSE
              sq1=sqrt2
           ENDIF
           a=qq(i)
           b=qq(i+1)
           c=qq(i+2)
           c1=c_kub(ABS(lm(1,i)),lm(2,i))
           c2=c_kub(ABS(lm(1,i)),lm(2,i)+4)
           c3=c_kub(ABS(lm(1,i)),lm(2,i)+8)
           qq(i)=a*c1*c1/sq1 + b*c1*c2/sq1 + &
                c*c1*c3/sq1
           qq(i+1)=a*c1*c2/sqrt2 + b*c2*c2/sqrt2  &
                + c*c2*c3/sqrt2
           qq(i+2)=a*c1*c3/sqrt2 + b*c2*c3/sqrt2  &
                + c*c3*c3/sqrt2                           
           i=i+3
           
        ELSE
           WRITE(6,*) 'UNCORRECT LM LIST FOR CUBIC STRUCTURE'
           WRITE(6,*) 'MULTSU.F'
           WRITE(6,'(a2,i4,a3,i4)') 'L=',lm(1,i),' M=', &
                lm(2,i)
!            call abort_parallel
           STOP
        ENDIF        
     END DO
  RETURN                                                            
END SUBROUTINE MULTsu


subroutine read_vns()
!This subroutine reads the non-spherical potential to calculate the non-spherical part of the core correction stress.
!This is being called from outside of jatom loop in hfsd.f
   use struct  , only: nat, jri
   use core_nsp, only: vns_st, lmmax_st, lm11_st
   implicit none
   include 'param.inc'
   integer    :: insv, ll_p, mm_p, lm(2, ncom+3), ri, nrp, lm1p, ii, &
                 jatom
   real*8     :: vtmp(nrad), val_real(nrad), dval_real(nrad), r1(nrad)  , &
                 out1, out2, one
   
      allocate( vns_st(nrad,ncom+3,nat) )
      allocate( lmmax_st(nat) ); allocate( lm11_st(2,ncom+3,nat) ) 
            
      vtmp       = 0.d0
      vns_st = (0.d0,0.d0)
      
!       rewind(19)
      read(19,'(//)',iostat=insv)
      IF(insv .ne. 0) STOP 'Error in reading vns file "core.corr.f" '
   DO jatom = 1, nat   
      nrp = jri(jatom)
      read(19,'(3x)') 
      read(19,'(15x,i3//)') lmmax_st(jatom)                                           
      
      DO lm1p = 1, lmmax_st(jatom)
         read(19,'(15x,i3,5x,i2/)') ll_p, mm_p         
         lm(1,lm1p) = ll_p
         lm(2,lm1p) = mm_p      
         lm11_st(1,lm1p,jatom) = ll_p
         lm11_st(2,lm1p,jatom) = mm_p
         
         read(19,'(3x,4ES19.12)') (vtmp(ri), ri = 1, nrp)
                 
         DO ri = 1, nrp
            vns_st(ri,lm1p,jatom) = vtmp(ri)
         ENDDO                  
         
         read(19,'(/)')
      ENDDO
      read(19,'(///)')

   ENDDO   
   
return
end subroutine read_vns
