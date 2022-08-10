!!----------------------------------------------------------------------------------------------------
!! The coefficients calculated in this subprogram will be used in our following publication          |
!! They are used wherever the unit vector is expanded in term of spherical harmonics and the         |
!! derivative of the spherical harmonics need to be computed.                                        |
!!                                                                                                   |
!! Kamal Belbase, Andreas Tröster, and Peter Blaha                                                   |
!! Phys. Rev. B 104, 174113 – Published 29 November 2021                                             |
!!                                                                                                   |
!! c_matrix is equivalent to c_{\alpha t} of Eq.(C6)                                                 |
!! c_ast_lm is equivalent to c_{\alpha}^{s t}(l,m) of Eq.(C4)                                        |
!!----------------------------------------------------------------------------------------------------

  SUBROUTINE c_alpha_m(c_matrix)
     !This c matrix comes into the formalism when a radial unit vector is expanded interm of spherical harmonics
     
     IMPLICIT NONE
     COMPLEX*16            :: c_matrix(3,-1:1)
     COMPLEX*16, PARAMETER :: imag  = (0.D0,1.d0)
     COMPLEX*16, PARAMETER :: zeroc = (0.d0,0.d0)
     REAL*8, PARAMETER     :: pi    = acos(-1.d0)
     REAL*8, PARAMETER     :: two   = 2.D0
     REAL*8, PARAMETER     :: three = 3.d0
     REAL*8, PARAMETER     :: twopi = 2.D0*pi
    
       c_matrix(1,-1) = sqrt(twopi/three)
       c_matrix(1,0)  = zeroc
       c_matrix(1,1)  = -sqrt(twopi/three)

       c_matrix(2,-1) = imag*sqrt(twopi/three)    !our work doesn't have minus, kuppel does
       c_matrix(2,0)  = zeroc
       c_matrix(2,1)  = imag*sqrt(twopi/three)    !our work doesn't have minus, kuppel does

       c_matrix(3,-1) = zeroc
       c_matrix(3,0)  = sqrt(two*twopi/three)
       c_matrix(3,1)  = zeroc
  
  RETURN
  END SUBROUTINE c_alpha_m
  
      
  SUBROUTINE C_factor(ll,mm, c_ast_lm)
     !SUBROUTINE to calculate c_ast(l,m) which is being used to calculate the derivative of the spherical hermonics
     !This subprogram is called inside the loop of ll and mm. ll and mm are the l and m characters of the charge density
     ! or the potential
     
     IMPLICIT NONE
     INTEGER,    intent(in)   :: ll,mm     
     COMPLEX*16, intent(out)  :: c_ast_lm(3,-1:1,-1:1)
     COMPLEX*16, PARAMETER    :: zeroc = (0.d0,0.d0)
     COMPLEX*16, PARAMETER    :: imag  = (0.D0,1.D0)
     INTEGER,    PARAMETER    :: one = 1.D0, two = 2.D0, three = 3.D0, f4r = 4.D0      
     INTEGER                  :: index1, ss, tt, index, alpha, beta
     REAL*8                   :: delta(3,3), C1, C2, FF(6)

                    
     C1  = (ll-mm+one)*(ll-mm+two)
     C2  = f4r*(two*ll+one)*(two*ll+three)
     FF(1) = sqrt(C1/C2)

     C1 = (ll+mm-one)*(ll+mm)
     C2 = f4r*(two*ll-one)*(two*ll+one)
     FF(2) = sqrt(C1/C2)

     C1 = (ll+mm+one)*(ll+mm+two)
     C2 = f4r*(two*ll+one)*(two*ll+three)
     FF(3) = sqrt(C1/C2)

     C1 = (ll-mm-one)*(ll-mm)
     C2 = f4r*(two*ll-one)*(two*ll+one)
     FF(4) = sqrt(C1/C2)

     C1 = (ll+one)*(ll+one) - mm*mm
     C2 = (two*ll+one)*(two*ll+three)
     FF(5) = sqrt(C1/C2)

     C1 = ll*ll - mm*mm
     C2 = (two*ll-one)*(two*ll+one)
     FF(6) = sqrt(C1/C2)
     
                                                  
     delta  = 0
     DO alpha = 1, 3
        DO beta = 1, 3
           IF(alpha .EQ. beta) THEN
               delta(alpha,beta) = 1.d0
           ELSE
               delta(alpha,beta) = 0.d0
           ENDIF
        ENDDO
     ENDDO
     
     c_ast_lm(:,0,:) = zeroc
         
     DO alpha = 1, 3
           
        c_ast_lm(alpha,-1,-1) = -(ll+1)*(delta(alpha,1)+imag*delta(alpha,2))*FF(2)
        c_ast_lm(alpha,-1, 0) =  (ll+1)*delta(alpha,3)*FF(6)
        c_ast_lm(alpha,-1, 1) =  (ll+1)*(delta(alpha,1)-imag*delta(alpha,2))*FF(4)
              
        c_ast_lm(alpha, 1,-1) = -ll*(delta(alpha,1)+imag*delta(alpha,2))*FF(1)
        c_ast_lm(alpha, 1, 0) = -ll*delta(alpha,3)*FF(5)
        c_ast_lm(alpha, 1, 1) =  ll*(delta(alpha,1) - imag*delta(alpha,2))*FF(3)
                           
     ENDDO     
     
  RETURN
  END SUBROUTINE C_factor


