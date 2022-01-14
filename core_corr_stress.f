    SUBROUTINE core_corr(natm, nrp, r1, rhoc, jatom, index_a)!,cmb_index,rot_loc11,rot_ij11)
 !this subroutine calculates the core correction according to formula 6.49 with lm = 0
 !natm is non_equ_atm
 !nrp total number of the radial points upto Rmt
 !r1 radial points at each points (1 to nrp)
 !rhoc is the charge density same as one printed in case.clmcor
 !dv1p is the spherical potential
       use TestCases
       use struct, only:mult, dx
       use core_sp, only: dv1p
       implicit REAL*8 (A-H,O-Z)
       include 'param.inc'             
!NRAD specified in param.inc, which is 881       
                  
       integer, intent(in)    :: natm, nrp, jatom, index_a
       real*8, intent(in)     :: r1(nrad+141), rhoc(nrad+141)
       integer                :: ir,jatm,imu
       real*8                 :: dx1, dpot(nrad+141),value(nrad+141),V1(NRAD+141)
       real*8                 :: SQ4pi,one,two,TC,Pi,value2(nrad+141)
       integer                :: counter, ii
       real*8                 :: sum1,sum_p_last,tot_char,int_test                                       
       real*8                 :: tot_charge_781,tot_charge_922
        
       PI = 4.0D+0*ATAN(1.0D+0)             
       sq4pi=1.d0/sqrt(4.d0*pi) 
       one=1.d0
       two=2.d0
!        dx = log(r1(nrp)/r1(1))/(nrp-1)                     
       dx1 = dx(jatom)
        
       sum1 = 0.d0
       sum_p_last= 0.d0           
            
       counter = nrp + 141                             !141 extra points than the information given in the input file               
    !   DO jatm = 1, natm                               !no of diff atoms, when this subroutine is being called it is inside jatom loop
          DO imu  = 1, MULT(jatom)                             !diff no of same atom
               DO ir = 1, counter                      
                     V1(ir) = two*DV1P(ir)             ! H to Ry      
               ENDDO               
               
               call dergl(counter,r1,v1,dpot,.false.,g0,.false.)          
               
               DO ir = 1, counter                          
                  value(ir) = rhoc(ir)*dpot(ir)*r1(ir)
                  value2(ir)= rhoc(ir)*V1(ir)
               enddo
               
               call chargel2(r1,1.d0,one,value(1),dx1,nrp,TC)                     
               sum1 = sum1 + TC               
               
               call chargel2(r1,1.d0,one,value(1),dx1,counter,TC)                     
               sum_p_last = sum_p_last + TC               

               call chargel2(r1,1.d0,one,rhoc,dx1,counter,tot_char)                              
               
               call chargel2(r1,1.d0,one,value2(1),dx1,counter,int_test)
                              
           !    write(21,'(4es19.12)') (DV1P(ir) , ir = 1, nrp)
           !    write(21,'(4es19.12)') (rhoc(ir) , ir = 1, nrp)
          ENDDO
    !   enddo
       
       !The following subroutine contains the explicit calculation of core correction stress
       call core_corr_stress(natm,nrp,r1,rhoc,dpot,v1,jatom, index_a)    !this is for full core correction stress
!        WRITE(21,*)'-------------from core_corr', nrp       
!        write(21,2021) nrp,r1(nrp),sum1/2.d0,value(nrp),DV1P(nrp)
       
!        write(21,2022) counter,tot_char,sum_p_last/2.d0     
       
!        write(21,2025) counter,int_test                   !to calculate int (r°2*rhoc*Veff)
!        write(21,'("r at 781 and r at 781+141",2x,2f10.5)') r1(nrp),r1(counter)
!        WRITE(21,*)'--------------------------'
      
       call chargel2(r1,1.d0,one,rhoc,dx,nrp,tot_charge_781)
       call chargel2(r1,1.d0,one,rhoc,dx,counter,tot_charge_922) 
       
!        WRITE(21,'(a,f10.5,a,f10.5)') ':tot_charge_781',tot_charge_781,':tot_charge_922',tot_charge_922
       
!sum1       :-is the core kinetic energy at the muffin-tin sphere boundary at 781
!sum_p_last :- is the core kinetic energy at the last points of the radial mesh, 
!which is extended well beyond the muffin-tin sphere of a host atom. 921
!tot_char   :- is the core density when integrate upto 921
!int_test   :- integral which appeared as rho*v_{eff} upto 921

 2032 FORMAT(50X,I2,//)       
 1980 FORMAT(3X)
 2000 FORMAT(16X,I2,//)
 2031 FORMAT(/)
 2021 FORMAT(':core_kin_at  ',I8, 4(3X,ES19.12))
 2022 FORMAT(':total_core charge and cor_kin at ',I6,2(3x,ES19.12))
 2025 FORMAT(':core_int_at  ',I8, 3X,1ES19.12)
 2033 FORMAT(///)
 2040 FORMAT(///////////////////)
 2041 FORMAT(31X,f8.5)
    RETURN
    END  
    
    SUBROUTINE core_corr_stress(natm,nrp,r1,rhoc,dpot,v1,jatom, index_a)
      use TestCases
      use struct, only : nat, iatnr, mult
      use core_nsp, only: lm11_st, lmmax_st, vns_st, sum_tot
      implicit REAL*8 (A-H,O-Z)
      include 'param.inc'
      integer             :: natm,nrp,jatm,mu,  ir, jatom, index_a
      REAL*8              :: dpot(nrad+141),r1(nrad+141),rhoc(nrad+141),          &
                             V1(nrad+141),value1(nrad+141)                                   
      REAL*8              :: one,two, pi, dx,  TC,TC1, sq4pi, sqrt2, sq4_pi 
                                  
      real*8              :: cor_corr_sph(1:9), cor_corr_ns1(1:9), &
                             cor_corr_ns2(1:9), cor_corr_tot(1:9), &
                             cor_corr_ns1_buf(1:9), cor_corr_ns2_buf(1:9)
                             
      integer             :: t,tp,s, alpha,beta,index2, lm1p,lm1, lmmax_p, ll_p, mm_p, counter,&
                             ri,insv, lmmax2,lmmax3, lmtot(2,ncom+3), lmtot1
      real*8,  external   :: GAUNT
      integer, external   :: NOTRI
      complex*16          :: cabt(3,-1:1),ca_st_lm(3,-1:1,-1:1), zeroc, zero, imag1, imag2, fac_st(ncom+3,nat)
      real*8              :: dvtmp(nrad+141),g0, gaunt1, int_out, vtmp(nrad+141)            
!       real(kind=16)        ::                                     
      complex*16           :: fc(ncom+3,nat)
      integer              :: lm(2,ncom+3), lmmax22, lm11(2,ncom+3)
      real*8               :: val_real(nrad+141),val_cmplx(nrad+141),value(nrad+141), &
                              dval_real(nrad+141),dval_cmplx(nrad+141), out1, out2 
      complex*16           :: value_c(nrad+141), out_c,int_ns1(ncom+3),int_ns2(ncom+3), &
                              vtmp_cmplx(nrad+141,ncom+3,nat), buf_sum(1:9)
     !=========here what we have===========
     !v1   = spherical potential
     !dpot = derivative of v1 
     !rhoc = spherical density          
     !cabt     = c-matrix coefficient during unit vectors expansion in terms of spherical harmonics
     !ca_st_lm = c-matrix coefficient the derivative of the spherical harmoincs
            
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
                           
      cor_corr_sph  = zero!zeroc      !spherical part
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
                 
         !Uncomment this if int_ns1 and int_ns2 is used during the calculation
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
      
   DO mu = 1, MULT(jatom)               !loop over same type of atoms, we don't loop over but multiply at the end      
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
            
      !first part of non-spherical core contribution      
      cor_corr_ns1_buf = 0.0
      buf_sum  = (0.0d0, 0.d0)
      DO lm1p = 1, lmtot1
         ll_p = lmtot(1,lm1p)
         mm_p = lmtot(2,lm1p)                                                       
         
         !if we use commented line then replace counter by nrp
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
         
!          IF( mm_p .ne. 0) CYCLE
         
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
      
      !second part of non-spherical core contribution
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
         
!          IF( mm_p .ne. 0) CYCLE
         
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
     
! !-------Below modification is done because of problem in Rutile if A^inv T A but gives the same result 

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
!     MULTSU COMBINES THE MULTIPOLMOMENTS QQ                            
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
!The subroutine reads the non-spherical potential to calculate the non-spherical part of the core correction stress.
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
