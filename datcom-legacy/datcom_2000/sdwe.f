      SUBROUTINE SDWE(X, Y, Z, AB, TAPER, SDW)
C
C***  DOWNWASH AT SUPERSONIC SPEEDS
C***  FIGURE 4.4.1-76E
C
      DIMENSION ROUTID(2), QE1(3), QE2(3), QE3(3), QE4(3)
      DIMENSION XE1(9), YE1(20), XE2(7), YE2(10)
      DIMENSION XE3(4), YE3(4), XE4(8), YE4(15)
      DIMENSION LGH(2), VAR(2)
C
      DATA ROUTID / 4HSDWE, 4H     /
      DATA QE1    / 4H4.4., 4H1-76, 4HE1   /
      DATA QE2    / 4H4.4., 4H1-76, 4HE2   /
      DATA QE3    / 4H4.4., 4H1-76, 4HE3   /
      DATA QE4    / 4H4.4., 4H1-76, 4HE4   /
C
C***           FIGURE 4.4.1-76E1
C
      DATA XE1 /
     1           1.0, 1.4, 1.8, 2.2, 2.4,
     2           0.0, 0.1, 0.3, 0.5 /
      DATA YE1 /
     1           1.68, 2.77, 3.38, 3.60, 3.62,
     2           1.05, 2.11, 2.78, 3.00, 3.00,
     3           0.00, 1.20, 1.91, 2.23, 2.30,
     4           0.00, 0.20, 1.10, 1.45, 1.51 /
C
C***           FIGURE 4.4.1-76E2
C
      DATA XE2 /
     1           0.0, 1.0, 2.0, 3.0, 4.0,
     2           4.0, 5.0 /
      DATA YE2 /
     1           0.00, 0.92, 1.97, 3.00, 4.10,
     2           0.00, 1.63, 3.36, 5.20, 7.30 /
C
C***           FIGURE 4.4.1-76E3
C
      DATA XE3 /
     1           0.00, 5.50,
     2           0.25, 0.50 /
      DATA YE3 /
     1           0.0, 5.5,
     2           0.0, 5.0 /
C
C***           FIGURE 4.4.1-76E4
C
      DATA XE4 /
     1           0.0, 0.6, 1.0, 1.5, 6.0,
     2           0.0, .15, .30 /
      DATA YE4 /
     1           0.00, .120, .200, .300, 1.20,
     2           .040, .040, .180, .280, 1.14,
     3           .040, .040, .120, .260, .960 /
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76E1
C
      VAR(1) = X
      VAR(2) = Z
      LGH(1) = 5
      LGH(2) = 4
      CALL INTERX(2, XE1, VAR, LGH, YE1, DUM, 5, 20,
     1            2, 0, 0, 0, 1, 2, 0, 0, QE1, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76E2
C
      VAR(1) = DUM
      VAR(2) = AB
      LGH(1) = 5
      LGH(2) = 2
      CALL INTERX(2, XE2, VAR, LGH, YE2, DUM, 5, 10,
     1            0, 2, 0, 0, 1, 2, 0, 0, QE2, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76E3
C
      VAR(1) = DUM
      VAR(2) = TAPER
      LGH(1) = 2
      LGH(2) = 2
      CALL INTERX(2, XE3, VAR, LGH, YE3, DUM, 2, 4,
     1            0, 1, 0, 0, 1, 1, 0, 0, QE3, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76E4
C
      VAR(1) = DUM
      VAR(2) = Y
      LGH(1) = 5
      LGH(2) = 3
      CALL INTERX(2, XE4, VAR, LGH, YE4, SDW, 5, 15,
     1            0, 1, 0, 0, 1, 1, 0, 0, QE4, 3, ROUTID)
      RETURN
      END
