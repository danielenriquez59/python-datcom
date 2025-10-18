      SUBROUTINE SDWC(X, Y, Z, AB, SDW)
C
C***  DOWNWASH AT SUPERSONIC SPEEDS
C***  FIGURE 4.4.1-76C (TAPER = 0.25, 0.50)
C
      DIMENSION SDW(2)
      DIMENSION ROUTID(2), QC1(3), QC2(3), QC3(3), QC4(3)
      DIMENSION XC1(11), YC1(28), XC2(12), YC2(12)
      DIMENSION XC3(5), YC3A(5), YC3B(5), XC4(4), YC4(4)
      DIMENSION LGH(2), VAR(2), CC(6)
C
      DATA ROUTID / 4HSDWC, 4H     /
      DATA QC1    / 4H4.4., 4H1-76, 4HC1   /
      DATA QC2    / 4H4.4., 4H1-76, 4HC2   /
      DATA QC3    / 4H4.4., 4H1-76, 4HC3   /
      DATA QC4    / 4H4.4., 4H1-76, 4HC4   /
C
C***           FIGURE 4.4.1-76C1
C
      DATA XC1 /
     1           .900, .950, 1.12, 1.40, 1.80, 2.20, 2.40,
     2           0.00, 0.10, 0.30, 0.50 /
      DATA YC1 /
     1           0.80, 1.00, 1.30, 1.63, 1.96, 2.08, 2.10,
     2           0.00, 0.70, 1.00, 1.39, 1.70, 1.81, 1.80,
     3           0.00, 0.00, 0.33, 0.85, 1.20, 1.40, 1.40,
     4           0.00, 0.00, 0.00, 0.38, 0.79, 0.99, 1.00 /
C
C***           FIGURE 4.4.1-76C2
C
      DATA XC2 /
     1           0.00, 2.50, 4*0.,
     2           2.70, 3.20, 5.40, 6.40, 10.7, 12.8 /
      DATA YC2 /
     1           0.0, 4.60,
     2           0.0, 3.55,
     3           0.0, 2.90,
     4           0.0, 2.50,
     5           0.0, 1.84,
     6           0.0, 1.59 /
C
C***           FIGURE 4.4.1-76C3
C
      DATA XC3  / 0.0, 3.00, 3.50, 3.75, 4.00 /
      DATA YC3A / 0.0, 3.00, 3.50, 3.75, 4.00 /
      DATA YC3B / 0.0, 2.35, 3.00, 3.50, 3.50 /
C
C***           FIGURE 4.4.1-76C4
C
      DATA XC4 /
     1           0.0, 4.00,
     2           0.0, .30 /
      DATA YC4 /
     1           0.0, .780,
     2           0.0, .820 /
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76C1
C
      VAR(1) = X
      VAR(2) = Z
      LGH(1) = 7
      LGH(2) = 4
      CALL INTERX(2, XC1, VAR, LGH, YC1, DUM, 7, 28,
     1            2, 0, 0, 0, 1, 2, 0, 0, QC1, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76C2
C
      VAR(1) = DUM
      VAR(2) = AB
      LGH(1) = 2
      LGH(2) = 6
      CALL INTERX(2, XC2, VAR, LGH, YC2, DUM, 6, 12,
     1            0, 2, 0, 0, 1, 2, 0, 0, QC2, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76C3 AND 76C4 AT TAPER = 0.25
C
      IN = 0
      CALL TBFUNX(DUM, DUM2, DYDX, 5, XC3, YC3A, CC, IN,
     1            MI, NG, 0, 1, QC3, 3, ROUTID)
C
      VAR(1) = DUM2
      VAR(2) = Y
      LGH(1) = 2
      LGH(2) = 2
      CALL INTERX(2, XC4, VAR, LGH, YC4, SDW(1), 2, 4,
     1            0, 1, 0, 0, 1, 1, 0, 0, QC4, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76C3 AND 76C4 AT TAPER = 0.50
C
      IN = 0
      CALL TBFUNX(DUM, DUM2, DYDX, 5, XC3, YC3B, CC, IN,
     1            MI, NG, 0, 1, QC3, 3, ROUTID)
C
      VAR(1) = DUM2
      VAR(2) = Y
      LGH(1) = 2
      LGH(2) = 2
      CALL INTERX(2, XC4, VAR, LGH, YC4, SDW(2), 2, 4,
     1            0, 1, 0, 0, 1, 1, 0, 0, QC4, 3, ROUTID)
      RETURN
      END
