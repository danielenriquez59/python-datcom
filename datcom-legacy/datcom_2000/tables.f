      SUBROUTINE TABLES(B, MACH, ALPHA, NDM)
C
C***  READ MACH TABLES OF CD EQUATION  REGESSION COEFFICIENTS
C
      DIMENSION B(16), BT(2,2,16), M(14)
      REAL MACH, M
      LOGICAL NDM
C
      DATA M /0.00, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00, 1.10,
     1        1.20, 1.30, 1.40, 1.50, 2.00, 2.50 /
C
      ALP = ABS(ALPHA)
      NDM = .FALSE.
      IF(MACH .LE. 0.9) IAM = 18
      IF(MACH .GT. 0.9 .AND. MACH .LT. 1.0) IAM = 11
      IF(MACH .GE. 1.0 .AND. MACH .LT. 1.1) IAM = 12
      IF(MACH .GE. 1.1) IAM = 15
      IA = ALP
      ALP1 = IAM
      IF(ALP .GT. ALP1) NDM = .TRUE.
      IF(NDM) GO TO 1030
        IF(IA .EQ. IAM) IA = IAM-1
        IA1 = IA+1
        IA2 = IA+2
        ALP1 = IA
        DO 1010 I=2,14
           IF(MACH .GE. M(I-1) .AND. MACH .LE. M(I)) GO TO 1020
 1010   CONTINUE
        NDM = .TRUE.
 1020   CONTINUE
        IM1 = I-1
        IM2 = I
 1030 CONTINUE
      IF(NDM) GO TO 1070
       I = 0
       DO 1050 IM=IM1,IM2
          J = 0
          I = I+1
          DO 1040 IA=IA1,IA2
             J = J+1
             IF(IM .LE. 4)                  CALL TBSUB(BT,IM,IA,I,J)
             IF(IM .GE. 5 .AND. IM .LE. 10) CALL TBTRN(BT,IM,IA,I,J)
             IF(IM .GE. 11)                 CALL TBSUP(BT,IM,IA,I,J)
 1040     CONTINUE
 1050  CONTINUE
C
       DO 1060 I=1,16
          BA = BT(1,1,I)+(BT(1,2,I)-BT(1,1,I))*(ALP-ALP1)
          BB = BT(2,1,I)+(BT(2,2,I)-BT(2,1,I))*(ALP-ALP1)
          B(I) = BA+(BB-BA)*(MACH-M(IM1))/(M(IM2)-M(IM1))
 1060  CONTINUE
 1070 CONTINUE
      RETURN
      END
