      SUBROUTINE SDDVC(X, ABETA, SWEPTE, SWEPR, TAPR, TAPER, DHB, ICASE)
C
C***  FIND DOWNWARD DISPLACEMENT OF VORTEX CORE
C***  FIGURE 4.4.1-74
C
      DIMENSION X(2), DHB(2), TAPER(4), DDVC(4)
      DIMENSION ROUTID(2), QA(3), QB(3), QC(3), QD(3)
      DIMENSION LGH(2), VAR(2), CC(6)
      DIMENSION XA(7), YA(10), XB(12), YB1(32), YB2(32)
      DIMENSION XC(8), YC1(15), YC2(15), XD(9), YD(20)
      LOGICAL   FLAG
C
      DATA ROUTID / 4HSDDV, 4HC    /
      DATA QA     / 4H4.4., 4H1-74, 4HA    /
      DATA QB     / 4H4.4., 4H1-74, 4HB    /
      DATA QC     / 4H4.4., 4H1-74, 4HC    /
      DATA QD     / 4H4.4., 4H1-74, 4HD    /
C
C***            FIGURE 4.4.1-74A
C
      DATA XA  /
     1          .5, 1.5, 3.0, 5.0, 6.0,
     2          8., 12. /
      DATA YA  /
     1          .5, .81, 1.50, 2.47, 2.93,
     2          .3, .61, 1.21, 2.00, 2.38 /
C
C***            FIGURE 4.4.1-74B
C
      DATA XB  /
     1          .5, .75, 1.50, 2.00, 2.50, 3.00, 3.50, 6.00,
     2          2., 4., 8., 12. /
      DATA YB1 /
     1          1.25, 1.60, 1.60, 2.07, 2.96, 3.81, 4.67, 9.03,
     2          .700, .750, 1.12, 1.49, 1.94, 2.51, 2.51, 2.51,
     3          .410, .460, .700, .900, 1.11, 1.32, 1.52, 2.58,
     4          .310, .330, .500, .630, .800, .960, 1.11, 1.89 /
      DATA YB2 /
     1          1.10, 1.50, 1.50, 1.96, 2.75, 3.56, 4.38, 8.33,
     2          .600, .670, 1.00, 1.30, 1.63, 2.08, 2.52, 2.52,
     3          .310, .340, .540, .700, .860, 1.02, 1.18, 1.94,
     4          .220, .221, .320, .450, .580, .700, .810, 1.40 /
C
C***            FIGURE 4.4.1-74C
C
      DATA XC  /
     1          .5, 1.0, 2.0, 3.5, 6.0,
     2          3.2, 6.4, 12.8 /
      DATA YC1 /
     1          .910, 1.03, 1.46, 2.41, 4.16,
     2          .500, .600, .940, 1.58, 2.66,
     3          .250, .320, .620, 1.08, 1.81 /
      DATA YC2 /
     1          .910, 1.03, 1.40, 2.31, 3.99,
     2         .500, .560, .880, 1.46, 2.41,
     3          .250, .320, .570, .960, 1.59 /
C
C***            FIGURE 4.4.1-74D
C
      DATA XD  /
     1          1.0, 1.5, 2.0, 3.75, 6.0,
     2          2., 4., 8., 12.  /
      DATA YD  /
     1          1.00, 1.41, 1.88, 4.00, 6.77,
     2          .510, .600, .710, 1.37, 2.20,
     3          .290, .300, .370, .660, 1.02,
     4          .130, .120, .200, .400, .690 /
C
      DO 1020 I=1,2
C
C***  FIGURE 4.4.1-74A -- TAPER = 0.0
C
         VAR(1) = X(I)
         VAR(2) = ABETA
         LGH(1) = 5
         LGH(2) = 2
         CALL INTERX(2, XA, VAR, LGH, YA, DDVC(1), 5, 10,
     1               1, 1, 0, 0, 1, 1, 0, 0, QA, 3, ROUTID)
C
C***  FIGURE 4.4.1-74B -- TAPER = 0.25, 0.50 -- SWEEP(C/2) = 0.0
C
         FLAG = (ICASE .EQ. 2) .OR. (ICASE .EQ. 3)
         IF(.NOT. FLAG) GO TO 1000
           LGH(1) = 8
           LGH(2) = 4
           CALL INTERX(2, XB, VAR, LGH, YB1, DDVC(2), 8, 32,
     1                 2, 2, 0, 0, 2, 2, 0, 0, QB, 3, ROUTID)
           CALL INTERX(2, XB, VAR, LGH, YB2, DDVC(3), 8, 32,
     1                 2, 2, 0, 0, 2, 2, 0, 0, QB, 3, ROUTID)
          IF(ICASE .EQ. 3) DDVC2 = DDVC(2)
          IF(ICASE .EQ. 3) DDVC3 = DDVC(3)
 1000    CONTINUE
C
C***  FIGURE 4.4.1-74C -- TAPER = 0.25, 0.50 -- SWEEP(TE) = 0.0
C
         FLAG = ICASE .EQ. 2
         IF(FLAG) GO TO 1010
           LGH(1) = 5
           LGH(2) = 3
           CALL INTERX(2, XC, VAR, LGH, YC1, DDVC(2), 5, 15,
     1                 2, 2, 0, 0, 1, 2, 0, 0, QC, 3, ROUTID)
           CALL INTERX(2, XC, VAR, LGH, YC2, DDVC(3), 5, 15,
     1                 2, 2, 0, 0, 1, 2, 0, 0, QC, 3, ROUTID)
           IF(ICASE .EQ. 3) DDVC(2) = DDVC2+SWEPTE*(DDVC(2)-DDVC2)/SWEPR
           IF(ICASE .EQ. 3) DDVC(3) = DDVC3+SWEPTE*(DDVC(3)-DDVC3)/SWEPR
 1010    CONTINUE
C
C***  FIGURE 4.4.1-74D -- TAPER = 1.0
C
         LGH(1) = 5
         LGH(2) = 4
         CALL INTERX(2, XD, VAR, LGH, YD, DDVC(4), 5, 20,
     1               2, 2, 0, 0, 1, 2, 0, 0, QD, 3, ROUTID)
C
C***  CALCULATE DHB
C
         IN = 0
         CALL TBFUNX(TAPR, DHB(I), DYDX, 4, TAPER, DDVC, CC,
     1               IN, MI, NG, 0, 0, 4HDHB , 1, ROUTID)
 1020 CONTINUE
      RETURN
      END
