      SUBROUTINE FIG60B(BETA, BTANA, CNAA)
C ***
C *** INTERPOLATE DATCOM FIGURE 4.1.3.3-60B FOR FIN CNAA
C ***
C *** REFERENCE   DATCOM SECTION 4.1.3.3
C ***
C *** LIMITATIONS   SUPERSONIC FIN LEADING EDGE
C ***               ATTACHED SHOCK
C ***               STRAIGHT TAPERED FINS
C ***
C *** INPUT   BETA    MACH SIMILARITY PARAMETER, SQRT(MACH**2-1.0)
C ***         BTANA   BETA*TAN(ALPHA), ALPHA = ANGLE OF ATTACK
C ***
C *** OUTPUT  CNAA    CN PER SIN(ALPHA)**2
C ***
C ***
      DIMENSION XA(13), YA(9), Z60B(117), BTA(13)
      DIMENSION Q60B(3), ROUT(2), C(6)
C ***
C *** CONSTANT DATA
C ***
      DATA DUM  / 1.0 /
      DATA Q60B / 4H4.1., 4H3.3-, 4H60B  /
      DATA ROUT / 4HFIG6, 4H0B   /
C ***
C *** FIGURE 4.1.3.3-60B
C ***
      DATA XA / 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
     *          2.0, 2.2, 2.4 /
      DATA YA / 1.25, 1.50, 1.75, 2.00, 2.50, 3.0, 4.0, 5.0, 20.0 /
      DATA Z60B / 0.000,  .047,  .089,  .127,  .160,  .188,  .211,
     *             .231,  .247,  .258,  .266,  .269,  .270,
     *            0.000,  .104,  .192,  .259,  .311,  .357,  .395,
     *             .426,  .452,  .474,  .492,  .505,  .515,
     *            0.000,  .112,  .215,  .308,  .405,  .499,  .572,
     *             .624,  .667,  .704,  .735,  .760,  .780,
     *            0.000,  .132,  .270,  .400,  .521,  .637,  .732,
     *             .820,  .878,  .925,  .962,  .986, 1.000,
     *            0.000,  .183,  .353,  .534,  .700,  .873, 1.040,
     *            1.190, 1.300, 1.388, 1.461, 1.520, 1.563,
     *            0.000,  .223,  .447,  .652,  .865, 1.096, 1.325,
     *            1.525, 1.697, 1.842, 1.978, 2.105, 2.222,
     *            0.000,  .270,  .576,  .845, 1.106, 1.447, 1.803,
     *            2.222, 2.511, 2.707, 2.875, 3.014, 3.125,
     *            0.000,  .270,  .575,  .890, 1.170, 1.563, 1.960,
     *            2.473, 2.921, 3.333, 3.744, 4.148, 4.545,
     *            0.000,  .306,  .629,  .925, 1.232, 1.667, 2.262,
     *            3.004, 4.065, 5.814, 9.259, 20.83, 1000. /
C ***
C *** CONSTRUCT CURVE OF CNAA VS BETA*TAN(ALPHA) AT INPUT BETA
C ***           XA   - CNAA CURVE
C ***           YA   - BETA CURVE
C ***           Z60B - BETA*TAN(ALPHA) VS CNAA AND BETA
C ***
      BTA(1) = 0.0
      DO 1000 I=2,13
        CALL TLINEX(YA, XA, Z60B, 9, 13, BETA, XA(I), BTA(I),
     1              1, 1, 0, 0, Q60B, 3, ROUT)
 1000 CONTINUE
C
C *** LOOK-UP FOR CNAA
C
      BTAN = BTANA
      IF(BTAN .GT. BTA(13)) BTAN = BTA(13)
      CALL TBFUNX(BTAN, CNAA, DYDX, 13, BTA, XA, C, IN,
     1            MI, NG, 0, 0, Q60B, 3, ROUT)
      RETURN
      END
