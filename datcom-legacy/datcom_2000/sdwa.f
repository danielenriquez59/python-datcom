      SUBROUTINE SDWA(X, Y, Z, AB, SDW)
C
C***  DOWNWASH AT SUPERSONIC SPEEDS
C***  FIGURE 4.4.1-76A (TAPER = 0.00)
C
      DIMENSION ROUTID(2), QA1(3), QA2(3), QA3(3)
      DIMENSION XA1(10), YA1(24), XA2(6), YA2(6), XA3(10), YA3(10)
      DIMENSION LGH(2), VAR(2)
C
      DATA ROUTID / 4HSDWA, 4H     /
      DATA QA1    / 4H4.4., 4H1-76, 4HA1   /
      DATA QA2    / 4H4.4., 4H1-76, 4HA2   /
      DATA QA3    / 4H4.4., 4H1-76, 4HA3   /
C
C***           FIGURE 4.4.1-76A1
C
      DATA XA1 /
     1            .90, 1.0, 1.2, 1.4, 1.6, 2.4,
     2            0.0, 0.1, 0.3, 0.5  /
      DATA YA1 /
     1            .90, 1.18, 1.65, 1.95, 2.10, 2.34,
     2            .55, 0.84, 1.30, 1.60, 1.74, 2.00,
     3            0.0, 0.00, 0.50, 0.87, 1.10, 1.50,
     4            0.0, 0.00, 0.00, 0.30, 0.61, 1.00 /
C
C***           FIGURE 4.4.1-76A2
C
      DATA XA2 /
     1            0.0, 3.0, 0.0,
     2            4.0, 8.0, 12. /
      DATA YA2 /
     1            0.0, 4.45,
     2            0.0, 3.00,
     3            0.0, 2.50 /
C
C***           FIGURE 4.4.1-76A3
C
      DATA XA3 /
     1            0.0, 3.0, 0.0, 0.0, 0.0,
     2            0.0, 0.3, 0.4, 0.6, 0.8 /
      DATA YA3 /
     1            0.00, 0.60,
     2            -.04, 0.54,
     3            -.09, 0.44,
     4            -.17, 0.35,
     5            -.28, 0.26 /
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76A1
C
      VAR(1) = X
      VAR(2) = Z
      LGH(1) = 6
      LGH(2) = 4
      CALL INTERX(2, XA1, VAR, LGH, YA1, DUM, 6, 24,
     1            2, 0, 0, 0, 1, 2, 0, 0, QA1, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76A2
C
      VAR(1) = DUM
      VAR(2) = AB
      LGH(1) = 2
      LGH(2) = 3
      CALL INTERX(2, XA2, VAR, LGH, YA2, DUM, 3, 6,
     1            0, 1, 0, 0, 1, 1, 0, 0, QA2, 3, ROUTID)
C
C***  TABLE LOOK UP - FIGURE 4.4.1-76A3
C
      VAR(1) = DUM
      VAR(2) = Y
      LGH(1) = 2
      LGH(2) = 5
      CALL INTERX(2, XA3, VAR, LGH, YA3, SDW, 5, 10,
     1            0, 0, 0, 0, 1, 2, 0, 0, QA3, 3, ROUTID)
      RETURN
      END
