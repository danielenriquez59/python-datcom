      SUBROUTINE WBLIFT(A,B,DB,AIW,WB,BW,FACT,BODY,WING,WINGIN)
C
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /CONSNT/ PI, DEG, UNUSED, RAD
      COMMON /FLGTCD/ FLC(160)
      DIMENSION ROUTID(2),Q1210A(3),Q1210B(3),Q212A1(4),Q212A2(4),
     1          Q1212C(3),Q412B(3),Q412C(3)
      DIMENSION A(195),B(49),WB(39),BODY(101),WING(101),WINGIN(4)
      DIMENSION FACT(1),BW(1),ANG(20),CB(6),QCLINT(2)
      DIMENSION X10A(11),Y10A(11),C10A(6)
      DIMENSION X10B(11),Y10B(11),C10B(6)
      DIMENSION X12A1(11),Y12A1(11),C12A1(6)
      DIMENSION X12A2(11),Y12A2(11),C12A2(6)
      DIMENSION X12C(13),Y12C(13),C12C(6)
      DIMENSION XA12(5),XB12(7),Y412C(35),Y412B(35)
      DATA ROUTID /4HWBLI,4HFT  /
      DATA Q1210A /4H4.3.,4H1.2-,4H10-A/,
     1     Q1210B /4H4.3.,4H1.2-,4H10-B/,
     2     Q212A1 /4H4.3.,4H1.2-,4H12-A,4H1   /,
     3     Q212A2 /4H4.3.,4H1.2-,4H12-A,4H2   /,
     4     Q1212C /4H4.3.,4H1.2-,4H12-C/,
     5     Q412B  /4H4.3., 4H1.4-, 4H12B  /,
     6     Q412C  /4H4.3., 4H1.4-, 4H12C  /
C
C*********       KW(B) VS D/B  DATA FIGURE 4.3.1.2-10-A
C
      DATA X10A/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
      DATA Y10A/1.0,1.08,1.16,1.26,1.36,1.46,1.56,1.67,1.78,1.89,2.0/
      DATA I10A/0/
C
C ****         KB(W) VS D/B  DATA  FIGURE 4.3.1.2-10-B
C
      DATA X10B/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
      DATA Y10B/0.0,0.13,0.29,0.45,0.62,0.80,1.0,1.22,1.45,1.70,2.0/
      DATA I10B/0/
C
C ***          FIGURE 4.3.1.2-12-A1
C
      DATA X12A1/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
      DATA Y12A1/1.0,0.97,0.95,4*0.94,0.95,0.96,0.98,0.99/
      DATA I12A1/0/
C
C*********          FIGURE 4.3.1.2-12-A2          **********************
C
      DATA X12A2/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
      DATA Y12A2/0.0,0.11,0.21,0.31,0.41,0.51,0.60,0.70,0.80,0.90,1.0/
      DATA I12A2/0/
C
C*********          FIGURE 4.3.1.2-12-C           **********************
C
      DATA X12C/0.0,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.83,0.9,1.0/
      DATA Y12C/2*1.0,0.999,0.99,0.98,0.965,0.95,0.933,2*0.92,0.928,
     10.95,1.0/
      DATA I12C/0/
C
C*********          FIGURE 4.3.1.4-12             **********************
C
      DATA XA12 /1., 2., 4., 6., 12. /
      DATA XB12 / 0., 0.05, 0.10, 0.15, 0.20, 0.25, 0.30 /
      DATA Y412B /
     1    1.0, 0.998, 0.995, 0.989, 0.982, 0.968, 0.950,
     2    1.0, 0.995, 0.989, 0.978, 0.963, 0.937, 0.917,
     3    1.0, 0.992, 0.982, 0.965, 0.942, 0.900, 0.864,
     4    1.0, 0.967, 0.932, 0.911, 0.907, 0.945, 1.038,
     5    1.0, 0.983, 0.968, 0.962, 0.977, 1.027, 1.134 /
      DATA Y412C /
     1    1.0, 1.013, 1.017, 1.013, 1.000, 0.982, 0.956,
     2    1.0, 1.003, 1.001, 0.992, 0.973, 0.937, 0.894,
     3    1.0, 0.989, 0.973, 0.944, 0.898, 0.822, 0.745,
     4    1.0, 0.980, 0.943, 0.879, 0.758, 0.648, 0.613,
     5    1.0, 0.942, 0.845, 0.688, 0.594, 0.552, 0.532 /
      DATA IY / 0 /, QCLINT / 4HCLIN, 4HT    /
C
C     COMPUTE BODY-DIAMETER-TO-WING-SPAN RATIO
C
      SSPN=WINGIN(4)
      RATIO= DB/(2.0*SSPN)
      AAA=1.0
      DO 1000 J=1,NALPHA
        ANG(J) = FLC(J+22)+AIW-B(49)
 1000 CONTINUE
      DO 1040 J=1,NALPHA
         IF(RATIO.GT.0.80) GO TO 1020
C
C*********          FIGURE 4.3.1.2-10-A
C
         CALL TBFUNX(RATIO,WB(2),DYDX,11,X10A,Y10A,C10A,I10A,MI,NG,
     1                0,0,Q1210A,3,ROUTID)
C
C*********          FIGURE 4.3.1.2-10-B
C
         CALL TBFUNX(RATIO,WB(3),DYDX,11,X10B,Y10B,C10B,I10B,MI,NG,
     1               0,0,Q1210B,3,ROUTID)
         WB(4)= WB(2)*WING(101)*AAA
         WB(5)= WB(3)*WING(101)*AAA
         BW(101)= BODY(101)+WB(4)+WB(5)
         IF(AIW.NE.0.0) GO TO 1010
         GO TO 1030
C
C*********          FIGURE 4.3.1.2-12-A1
C
 1010    CALL TBFUNX(RATIO,WB(7),DYDX,11,X12A1,Y12A1,C12A1,I12A1,MI,NG,
     1               0,0,Q212A1,4,ROUTID)
C
C*********          FIGURE 4.3.1.2-12-A2
C
         CALL TBFUNX(RATIO,WB(8),DYDX,11,X12A2,Y12A2,C12A2,I12A2,MI,NG,
     1               0,0,Q212A2,4,ROUTID)
         WB(9)=WB(7)*WING(101)*AAA
         WB(10)=WB(8)*WING(101)*AAA
         WB(11)=WB(9)+WB(10)
         GO TO 1030
C
C****    CALCULATIONS WITH WING EXTENDING ENTIRE LENGTH OF BODY  *******
C****        (FOR LARGE BODY-DIAMETER-TO-WING-SPAN RATIOS)       *******
C****                                                            *******
C*********          FIGURE 4.3.1.2-12-C
C
 1020    CALL TBFUNX(RATIO,WB12,DYDX,13,X12C,Y12C,C12C,I12C,MI,NG,
     1               0,0,Q1212C,3,ROUTID)
         BW(101)= WB12*WING(101)*AAA
C
C*********          WING BODY LIFT IN NON-LINEAR            ************
C*********          ANGLE-OF-ATTACK RANGE                   ************
C
 1030    CONTINUE
C
C*********** LIFT COEFFICIENT (CL) CALCULATION *************************
C
         ALPEF=FLC(J+22)-B(49)+((WB(7)+WB(8))/(WB(2)+WB(3)))*AIW
         CALL TBFUNX(ALPEF,CLINT,DYDX,NALPHA,ANG,WING(21),
     1              CB,IY,MI,NG,1,1,QCLINT,2,ROUTID)
         BW(J+20)=BODY(J+20)+(WB(2)+WB(3))*CLINT
         BW(J+20)=BW(J+20)+FACT(J+1)*FACT(J+21)*B(J+22)*FACT(1)*
     1            WING(101)
C
C*********** NORMAL FORCE (CN) CALCULATION *****************************
C
         BW(J+60)=BW(J+20)*COS(B(J+22)/RAD)+BW(J)*SIN(B(J+22)/RAD)
C
C     ******* SEE SUBROUTINE WBCM FOR CA CALCULATION *******
C
 1040 CONTINUE
      IF(RATIO .GT. 0.30) GO TO 1050
      CALL TLINEX(XA12,XB12,Y412B,5,7,A(160),RATIO,WB(20),
     1           0,0,0,0,Q412B,3,ROUTID)
      CALL TLINEX(XA12,XB12,Y412C,5,7,A(160),RATIO,WB(21),
     1           0,0,0,0,Q412C,3,ROUTID)
      WB(22) = WB(20)*B(44)
      WB(23) = WB(21)*B(43)
 1050 CONTINUE
      RETURN
      END
