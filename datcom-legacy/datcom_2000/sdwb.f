      SUBROUTINE SDWB(X, Y, Z, AB, SDW)
C
C***  DOWNWASH AT SUPERSONIC SPEEDS
C***  FIGURE 4.4.1-76B (TAPER = 1.00)
C
      DIMENSION ROUTID(2), QB1(3), QB2(3), QB3(3)
      DIMENSION XB1(12), YB1(32), XB2(8), YB2(16), XB3(6), YB3(6)
      DIMENSION LGH(2), VAR(2)
C
      DATA ROUTID / 4HSDWB, 4H     /
      DATA QB1    / 4H4.4., 4H1-76, 4HB1   /
      DATA QB2    / 4H4.4., 4H1-76, 4HB2   /
      DATA QB3    / 4H4.4., 4H1-76, 4HB3   /
C
C***           FIGURE 4.4.1-76B1
C
      DATA XB1 /
     1           1.0, 1.2, 1.4, 1.5, 1.6, 2.0, 2.2, 2.4,
     2           0.0, 0.1, 0.3, 0.5 /
      DATA YB1 /
     1           0.32, 0.74, 1.21, 1.39, 1.44, 1.56, 1.62, 1.77,
     2           0.11, 0.60, 1.10, 1.25, 1.30, 1.42, 1.50, 1.62,
     3           0.00, 0.20, 0.76, 1.00, 1.10, 1.34, 1.42, 1.47,
     4           0.00, 0.00, 0.00, 0.50, 0.75, 1.06, 1.12, 1.20 /
C
C***           FIGURE 4.4.1-76B2
C
      DATA XB2 /
     1           0.00, .8, 1.70, 3.0,
     2           2.00, 4., 8.00, 12. /
      DATA YB2 /
     1           0.00, .60, 4.43, 12.7,
     2           0.00, .80, 1.70, 3.00,
     3           0.00, .40, 0.80, 1.40,
     4           0.00, .27, 0.50, 0.94 /
C
C***           FIGURE 4.4.1-76B3
C
      DATA XB3 /
     1           0.0, 5.0, 0.0,
     2           0.0, .15, .30 /
      DATA YB3 /
     1           0.0, 1.00,
     2           0.0, 1.10,
     3           0.0, 1.20 /
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76B1
C
      VAR(1) = X
      VAR(2) = Z
      LGH(1) = 8
      LGH(2) = 4
      CALL INTERX(2, XB1, VAR, LGH, YB1, DUM, 8, 32,
     1            2, 0, 0, 0, 1, 2, 0, 0, QB1, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76B2
C
      VAR(1) = DUM
      VAR(2) = AB
      LGH(1) = 4
      LGH(2) = 4
      CALL INTERX(2, XB2, VAR, LGH, YB2, DUM, 4, 16,
     1            0, 2, 0, 0, 1, 2, 0, 0, QB2, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76B3
C
      VAR(1) = DUM
      VAR(2) = Y
      LGH(1) = 2
      LGH(2) = 3
      CALL INTERX(2, XB3, VAR, LGH, YB3, SDW, 3, 6,
     1            0, 0, 0, 0, 1, 2, 0, 0, QB3, 3, ROUTID)
      RETURN
      END
