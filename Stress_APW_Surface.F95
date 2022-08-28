SUBROUTINE stress_APW_surface(nemin,nemax,jatom,rmt,ALM,BLM,CLM,A_surf2,buf_APW)
   
   !===================================================================================
   ! This subroutine resides within the k=point, atomic index loop. This subroutine    !
   ! computes the surface term at the sphere boundary when APW+lo method is used.      !
   ! Wave functions are involved to cimputes this contribution and thereby the         !
   ! computational cost is little bit higher then those computed with charge density.  !
   ! An expression will be available in our second paper.                              !
   !===================================================================================
   
   ! Variables
   ! nemin         = lowest band index for a given k-point
   ! nemax         = maximum band index for a given k-point
   ! jatom         = atom index
   ! rmt           = rmt sphere for given jatom
   ! ALM, BLM, CLM = matching coefficient at the sphere boundary
   ! A_surf2       = surface term at the sphere boundary from the interstitial side
   ! buf_APW       = output integration (from_inside_sphere-from_outside_sphere)|_{at boundry}
   
   ! psi_is        = the part of wave function in interstitial 
   !                 this expand interms of spherical harmonics and bessel function using
   !                 Rayleigh expansion
   ! psi_ra        = the part of wave function in an atomic sphere
   ! Isurf_int ds  = surface integral
   ! buf_APW       = Isurf_int [psi_ra*(d psi_ra/dr_a) - psi_is*(d psi_is /dr_a) ] ds                   
   
   USE param   
   USE atspdt, ONLY: P,DP,PE,DPE
   USE lo    , ONLY: loor,ilo,plo,dplo
   USE xa    , ONLY: weight
   
   IMPLICIT NONE
   
   REAL*8 ,    INTENT(IN)   :: rmt
   INTEGER,    INTENT(IN)   :: nemin, nemax,jatom
   COMPLEX*16, INTENT(IN)   :: CLM((lmax2+1)*(lmax2+1),nume,nloat),          &
                               ALM((lmax2+1)*(lmax2+1),nume),                &
                               BLM((lmax2+1)*(lmax2+1),nume),                &                               
                               A_surf2((lmax2+1)*(lmax2+1),nume)
                               
   COMPLEX*16, INTENT(OUT)  :: buf_APW(3,3,nume)   
   
   COMPLEX*16               :: c_am(3,-1:1), buf_temp1(3,3,nume),                     &
                               buf_temp2(3,3,nume), buf_temp3(3,3,nume),              &
                               buf_temp4(3,3,nume), zeroc
                          
   INTEGER                  :: l0, lp, s0, v0, t0, tp, alpha, beta, m0, mp, num,       &
                               index1,index2,jlo,jlop, index, index3,lmda,
                               lmdap 
                           
   REAL*8                   :: gaunt1,gaunt2,rmt3,rmt2,tw_by_r                        
   INTEGER, EXTERNAL        :: NOTRI
   REAL*8 , EXTERNAL        :: GAUNT
   
   
   zeroc = (0.d0,0.0d0)
   call c_alpha_m(c_am)
   
   buf_temp1 = zeroc
   buf_temp2 = zeroc
   buf_temp3 = zeroc
   buf_temp4 = zeroc
   
   buf_APW   = zeroc        
   
   rmt3 = rmt*rmt*rmt
   rmt2 = rmt*rmt
   tw_by_r = 2.d0/rmt
   
   index1 = 0
   index2 = 0
   DO s0 = 0,2
      ! l character of spherical harmonics can't be negative
      IF(NOTRI(s0,1,1) .lt. 0) CYCLE  
      DO v0 = -s0,s0
         DO t0 = -1,1
            DO tp = -1,1
               IF(-v0+t0+tp .ne. 0) CYCLE
               index1 = 0
               DO l0 = 0,lmax2
                  DO m0 = -l0,l0
                     index1 = index1 + 1                     
                     index2 = 0
                     DO lp = 0,lmax2
                        DO mp = -lp,lp
                           index2 = index2+1
                           IF(NOTRI(l0,lp,s0) .lt. 0) CYCLE
                           IF(-m0+mp+v0 .ne. 0) CYCLE
                           
                           !Because of the surface integral of four spherical harmonics, there are
                           !the product of two gaunt number
                           
                           gaunt1 = GAUNT(s0,1,1,v0,t0,tp)
                           gaunt2 = GAUNT(l0,lp,s0,m0,mp,v0)
                           
                           DO alpha = 1, 3
                              DO beta = 1,3
                                 DO num = nemin,nemax
                                    buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                    c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                    (conjg(ALM(index1,num))*P(l0)+conjg(BLM(index1,num))*PE(l0))*&
                                    (ALM(index2,num)*DP(l0)+BLM(index2,num)*DPE(l0)-A_surf2(index2,num))
                                 ENDDO
                                 
                                 DO num = nemin,nemax
                                    DO jlo = 1,ilo(l0)
                                       IF(loor(jlo,l0)) THEN
                                         buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                         c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                         conjg(clm(index1,num,jlo))*plo(jlo,l0)*   &
                                         (ALM(index2,num)*DP(l0)+BLM(index2,num)*DPE(l0)-A_surf2(index2,num))                                                                                  
                                       ENDIF
                                    ENDDO
                                 ENDDO
                                 
                                 DO num = nemin,nemax
                                    DO jlop = 1,ilo(lp)
                                       IF(loor(jlop,lp)) THEN
                                          buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                          c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                          clm(index2,num,jlop)*dplo(jlop,lp)*&
                                          (conjg(ALM(index1,num))*P(l0)+conjg(BLM(index1,num))*PE(l0))
                                       ENDIF
                                    ENDDO
                                 ENDDO
                                 
                                 DO num = nemin,nemax
                                    DO jlo = 1,ilo(l0)
                                       IF(loor(jlo,l0))THEN
                                          DO jlop = 1,ilo(lp)
                                             IF(loor(jlop,lp)) THEN
                                                buf_APW(alpha,beta,num) = buf_APW(alpha,beta,num) + &
                                                c_am(alpha,t0)*c_am(beta,tp)*gaunt1*gaunt2*rmt2*&
                                                conjg(clm(index1,num,jlo))*plo(jlo,l0)*clm(index2,num,jlop)*dplo(jlop,lp)
                                             ENDIF
                                          ENDDO
                                       ENDIF
                                    ENDDO
                                 ENDDO
                                 
                              ENDDO
                           ENDDO                                                                                 
                           
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   
   
RETURN
END SUBROUTINE
