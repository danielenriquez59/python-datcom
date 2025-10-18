      SUBROUTINE SDWD(X, Y, Z, AB, SDW)
C
C***  DOWNWASH AT SUPERSONIC SPEEDS
C***  FIGURE 4.4.1-76D (TAPER = 0.00, 0.25, 0.50)
C
      DIMENSION SDW(3)
      DIMENSION ROUTID(2), QD1(3), QD2(3), QD3(3), QD4(3)
      DIMENSION XD1(9), YD1(20), XD2(9), YD2(20)
      DIMENSION XD3(2), YD3A(2), YD3B(2), YD3C(2), XD4(8), YD4(15)
      DIMENSION LGH(2), VAR(2), CC(6)
C
      DATA ROUTID / 4HSDWD, 4H     /
      DATA QD1    / 4H4.4., 4H1-76, 4HD1   /
      DATA QD2    / 4H4.4., 4H1-76, 4HD2   /
      DATA QD3    / 4H4.4., 4H1-76, 4HD3   /
      DATA QD4    / 4H4.4., 4H1-76, 4HD4   /
C
C***           FIGURE 4.4.1-76D1
C
      DATA XD1 /
     1           1.0, 1.4, 1.8, 2.2, 2.4,
     2           0.0, 0.1, 0.3, 0.5 /
      DATA YD1 /
     1           1.68, 2.77, 3.38, 3.60, 3.62,
     2           1.05, 2.11, 2.78, 3.00, 3.00,
     3           0.00, 1.20, 1.91, 2.23, 2.30,
     4           0.00, 0.20, 1.10, 1.45, 1.51 /
C
C***           FIGURE 4.4.1-76D2
C
      DATA XD2 /
     1           0.0, .2, .3, 1.0, 4.0,
     2           2.0, 4., 8., 12. /
      DATA YD2 /
     1           0.00, .100, .200, 1.00, 7.00,
     2           0.00, .200, .300, 1.00, 4.00,
     3           0.00, .600, .800, 1.20, 2.40,
     4           0.00, .400, .600, 1.00, 1.60 /
C
C***           FIGURE 4.4.1-76D3
C
      DATA XD3  / 0.0, 5.5 /
      DATA YD3A / 0.0, 6.3 /
      DATA YD3B / 0.0, 5.5 /
      DATA YD3C / 0.0, 5.0 /
C
C***           FIGURE 4.4.1-76D4
C
      DATA XD4 /
     1           0.0, 0.6, 1.0, 1.5, 6.0,
     2           0.0, .15, .30 /
      DATA YD4 /
     1           0.00, .120, .200, .300, 1.20,
     2           .040, .040, .180, .280, 1.14,
     3           .040, .040, .120, .260, .960 /
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76D1
C
      VAR(1) = X
      VAR(2) = Z
      LGH(1) = 5
      LGH(2) = 4
      CALL INTERX(2, XD1, VAR, LGH, YD1, DUM, 5, 20,
     1            2, 0, 0, 0, 1, 2, 0, 0, QD1, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76D2
C
      VAR(1) = DUM
      VAR(2) = AB
      LGH(1) = 5
      LGH(2) = 4
      CALL INTERX(2, XD2, VAR, LGH, YD2, DUM, 5, 20,
     1            0, 2, 0, 0, 1, 2, 0, 0, QD2, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76D3 AND 76D4 AT TAPER = 0.00
C
      IN = 0
      CALL TBFUNX(DUM, DUM2, DYDX, 2, XD3, YD3A, CC, IN,
     1            MI, NG, 0, 1, QD3, 3, ROUTID)
C
      VAR(1) = DUM2
      VAR(2) = Y
      LGH(1) = 5
      LGH(2) = 3
      CALL INTERX(2, XD4, VAR, LGH, YD4, SDW(1), 5, 15,
     1            0, 1, 0, 0, 1, 1, 0, 0, QD4, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76D3 AND 76D4 AT TAPER = 0.25
C
      IN = 0
      CALL TBFUNX(DUM, DUM2, DYDX, 2, XD3, YD3B, CC, IN,
     1            MI, NG, 0, 1, QD3, 3, ROUTID)
C
      VAR(1) = DUM2
      VAR(2) = Y
      LGH(1) = 5
      LGH(2) = 3
      CALL INTERX(2, XD4, VAR, LGH, YD4, SDW(2), 5, 15,
     1            0, 1, 0, 0, 1, 1, 0, 0, QD4, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76D3 AND 76D4 AT TAPER = 0.50
C
      IN = 0
      CALL TBFUNX(DUM, DUM2, DYDX, 2, XD3, YD3C, CC, IN,
     1            MI, NG, 0, 1, QD3, 3, ROUTID)
C
      VAR(1) = DUM2
      VAR(2) = Y
      LGH(1) = 5
      LGH(2) = 3
      CALL INTERX(2, XD4, VAR, LGH, YD4, SDW(3), 5, 15,
     1            0, 1, 0, 0, 1, 1, 0, 0, QD4, 3, ROUTID)
      RETURN
      END
