!****************************************************************************    
! *   FILE         = MULTIPLE FILES                                         *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *  
!**************************************************************************** 
      PROGRAM BLASISUS_EQ
      USE GLOBAL
      IMPLICIT NONE


       DY= 0.25

!      BLASIUS EQUATION 2F''' + FF" = 0

!      BC AT ETA = 0, F = 0; F' = 0
       F(1) = 0.0D0   ! F
       G(1) = 0.0D0   ! F'


!      INITIAL GUESS1 F"(0)
       P=0.3359375D0
       H(1)=P    ! F"
       CALL RK
       PE1= G(N)-1.0D0

!       G(N) IS BOUNDARY CONDITION AT ETA --> INFINTY  ==> F'(INFINITY) = 1

      
!      INITIAL GUESS2 F"(0)
       Q = 0.010D0
       H(1) = Q    ! F"
       CALL RK
       PE2 = G(N)-1.0D0

!      CHOOSE INTIAL GUSESSES SUCH THAT PE1*PE2 = -VE TO SHOOT USING BISECTION METHOD

       PRINT*, PE1, PE2 

       IF (PE1*PE2 .LT. 0) THEN

           DO KM = 1, 100
             R = (P+Q)/2.0D0
             H(1) = R
             CALL RK

             PRINT*, 'ITERATION-', KM, ',INITIAL GUESS F"(0)-', R
            

             IF ( PE1 .GT. 0) THEN
               IF ((G(N)-1.0D0) .LT. 0) THEN
                  Q=R
               ELSE IF  ((G(N)-1.0D0) .GT. 0) THEN
                  P=R
               ELSE IF ((G(N)-1.0D0) .EQ. 0)  THEN
                  EXIT
               END IF
             ELSE
               IF ((G(N)-1.0D0) .LT. 0) THEN
                  P=R
               ELSE IF  ((G(N)-1.0D0) .GT. 0) THEN
                  Q=R
               ELSE IF ((G(N)-1.0D0) .EQ. 0)  THEN
                  EXIT
               END IF
             END IF

           END DO
           
           PRINT*, 'SOLUTION CONVERGED'
           PRINT*, 'ETA', '  AND', ' F''= U/V'

           DO I =1, N
              PRINT*, I, G(I)
           END DO

        ELSE
           PRINT*, 'BISECTION METHOD FAILED:MODIFY INITIAL GUESS'
           PRINT*, PE1, PE2
        END IF


      END PROGRAM BLASISUS_EQ


      SUBROUTINE RK
      USE GLOBAL
      DOUBLE PRECISION :: K1(3), K2(3), K3(3), K4(3)
      DOUBLE PRECISION :: A, B, C, F1


      DO I= 1, N-1
 
        A = F(I)
        B = G(I)
        C = H(I)

        DO J=1, 3
         K1(J) = F1(A, B, C, J)
        END DO

        A = F(I)+0.5D0*K1(1)*DY
        B = G(I)+0.5D0*K1(2)*DY
        C = H(I)+0.5D0*K1(3)*DY

        DO J =1,3
           K2(J) = F1(A, B, C, J)
        END DO

        A = F(I)+0.5D0*K2(1)*DY
        B = G(I)+0.5D0*K2(2)*DY
        C = H(I)+0.5D0*K2(3)*DY

        DO J =1,3
           K3(J) = F1(A, B, C, J)
        END DO

        A = F(I)+K3(1)*DY
        B = G(I)+K3(2)*DY
        C = H(I)+K3(3)*DY

        DO J =1,3
           K4(J) = F1(A, B, C, J)
        END DO

       
         F(I+1) = F(I) +DY*( K1(1) + 2.0D0*K2(1) + 2.0D0*K3(1) + K4(1))/6.0D0
         G(I+1) = G(I) +DY*( K1(2) + 2.0D0*K2(2) + 2.0D0*K3(2) + K4(2))/6.0D0
         H(I+1) = H(I) +DY*( K1(3) + 2.0D0*K2(3) + 2.0D0*K3(3) + K4(3))/6.0D0

      END DO

      RETURN
      END SUBROUTINE RK


      FUNCTION F1(Y1,Y2,Y3, W)
      USE GLOBAL
      DOUBLE PRECISION :: Y1,Y2,Y3, F1
      INTEGER :: W
      
      IF (W .EQ. 1) THEN
      F1 = Y2
      ELSE IF (W .EQ. 2) THEN
      F1 = Y3  
      ELSE 
      F1 = -Y1*Y3/2.0D0
      END IF
 
      RETURN
      END FUNCTION
 
       

