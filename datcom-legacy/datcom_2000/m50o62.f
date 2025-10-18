      SUBROUTINE M50O62                                                 M50O62
C                                                                       M50O62
C  AIRFOIL SECTION AERODYNAMICS AND GEOMETRIC PROPERTIES                M50O62
C                                                                       M50O62
      COMMON /OVERLY/IJKDUM(8),NOVLY                                    M50O62
      COMMON /IBODY/ PB, NACA(80), BF(232), CBAR                        M50O62
      COMMON /IWING/ PW, X(60)                                          M50O62
      COMMON / IHT / PH, XU(60),XL(60),YU(60),YL(60)                    M50O62
      COMMON / IVT / PV, YUN(60),YLN(60)                                M50O62
      COMMON / IVF / PVF, BVF(380)                                      M50O62
      COMMON / IBW / PBW, BW(380)                                       M50O62
      COMMON / IBH / PBH, BBH(380)                                      M50O62
      COMMON / IBV / PBV, BBV(380)                                      M50O62
      COMMON /IBWH/  PBWH, AI,ALO,CLI,ASEP,CMCO4,CLA0,CLA(20),CLMAX0,   M50O62
     1               CLMAX(20),BH(232)                                  M50O62
      COMMON /IBWHV/ PBWHV,RHO,TMAX,DELTAY,XOVC,TOVC,COVC,CXVC,CAMBER,  M50O62
     1               BV(232)                                            M50O62
      COMMON /WINGI/ WINGIN(100)                                        M50O62
      COMMON /HTI/   HTIN(154)                                          M50O62
      COMMON /VTI/   VTIN(154),TVTIN(8),VFIN(154)                       M50O62
      COMMON /FLGTCD/FLC(93)                                            M50O62
      COMMON /CONSNT/PI,DEG,UNUSED,RAD                                  M50O62
      LOGICAL NACAIN,CAMBER                                             M50O62
      INTEGER A4HN,A4HA,A4HC                                            M50O62
      DIMENSION THN(60),CAM(60)                                         M50O62
      EQUIVALENCE (THN(1),BBH(1)),(CAM(1),BBH(61))                      M50O62
      DATA A4HN,A4HA,A4HC/4HN   ,4HA   ,4HC   /                         M50O62
      NOVLY=50                                                          M50O62
      DO 1000 I=1,232                                                   M50O62
         BH(I)=BBH(I+130)                                               M50O62
         BV(I)=BBV(I+130)                                               M50O62
         BF(I)=BVF(I+130)                                               M50O62
 1000 CONTINUE                                                          M50O62
      IN=0                                                              M50O62
 1010 IN=IN+1                                                           M50O62
      CALL INIZ                                                         M50O62
      NACAIN=.FALSE.                                                    M50O62
      CAMBER=.FALSE.                                                    M50O62
      IF(IN .GT. 4) GO TO 1150                                          M50O62
      GO TO (1020,1030,1040,1050), IN                                   M50O62
 1020 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- HERE FOR WING                                                     M50O62
C                                                                       M50O62
      IF(WINGIN(15) .EQ. UNUSED) GO TO 1010                             M50O62
      CALL SECI(BW(131),CAMBER,ATYPE,WINGIN,L)                          M50O62
      RHO = WINGIN(62)                                                  M50O62
      GO TO 1060                                                        M50O62
 1030 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- HERE FOR HORIZONTAL TAIL                                          M50O62
C                                                                       M50O62
      IF(HTIN(15) .EQ. UNUSED) GO TO 1010                               M50O62
      CALL SECI(BH(1),CAMBER,ATYPE,HTIN,L)                              M50O62
      RHO = HTIN(62)                                                    M50O62
      GO TO 1060                                                        M50O62
 1040 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- HERE FOR VERTICAL TAIL                                            M50O62
C                                                                       M50O62
      IF(VTIN(15) .EQ. UNUSED) GO TO 1010                               M50O62
      CALL SECI(BV(1),CAMBER,ATYPE,VTIN,L)                              M50O62
      GO TO 1060                                                        M50O62
 1050 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- HERE FOR VENTRAL FIN                                              M50O62
C                                                                       M50O62
      IF(VFIN(15) .EQ. UNUSED) GO TO 1010                               M50O62
      CALL SECI(BF(1),CAMBER,ATYPE,VFIN,L)                              M50O62
      RHO = VTIN(62)                                                    M50O62
 1060 CONTINUE                                                          M50O62
      IF(IN .EQ. 1 .AND. ATYPE .LT. 0.) WRITE(6,1160)                   M50O62
      IF(IN .EQ. 2 .AND. ATYPE .LT. 0.) WRITE(6,1170)                   M50O62
      IF(IN .EQ. 3 .AND. ATYPE .LT. 0.) WRITE(6,1180)                   M50O62
      IF(IN .EQ. 4 .AND. ATYPE .LT. 0.) WRITE(6,1190)                   M50O62
      IF(NACA(1).EQ.A4HN .AND. NACA(2).EQ.A4HA .AND.                    M50O62
     1   NACA(3).EQ.A4HC .AND. NACA(4).EQ.A4HA)NACAIN=.TRUE.            M50O62
      IF(NACAIN .AND. ABS(ATYPE) .LE. UNUSED) GO TO 1090                M50O62
      IF(ABS(ATYPE).LE.UNUSED) GO TO 1010                               M50O62
      IVAL=ABS(ATYPE)-0.5                                               M50O62
      CALL XYCORD(IVAL,IN)                                              M50O62
 1070 CONTINUE                                                          M50O62
      DO 1080 M=2,L                                                     M50O62
         IF(THN(M) .GT. THN(M-1)) TOVC=THN(M)                           M50O62
         IF(THN(M) .GT. THN(M-1)) XOVC=X(M)                             M50O62
         IF(CAM(M) .GT. CAM(M-1)) COVC=CAM(M)                           M50O62
         IF(CAM(M) .GT. CAM(M-1)) CXVC=X(M)                             M50O62
 1080 CONTINUE                                                          M50O62
      IF(COVC .GT. UNUSED) CAMBER=.TRUE.                                M50O62
      IF(TMAX .EQ. UNUSED) TMAX=2.*TOVC                                 M50O62
      IF(RHO.EQ.UNUSED) RHO = (X(2)**2+THN(2)**2)/(2.0*X(2))            M50O62
      CALL DELY                                                         M50O62
      CALL THEORY                                                       M50O62
      CALL MAXCL                                                        M50O62
      GO TO 1100                                                        M50O62
 1090 CALL AIRFOL                                                       M50O62
      GO TO 1070                                                        M50O62
 1100 CONTINUE                                                          M50O62
      GO TO (1110,1120,1130,1140), IN                                   M50O62
 1110 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- SET CALCULATED WING PARAMETERS                                    M50O62
C                                                                       M50O62
      CALL SECO(WINGIN(1),CAMBER,ATYPE)                                 M50O62
      CALL CSLOPE                                                       M50O62
      GO TO 1010                                                        M50O62
 1120 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- SET CALCULATED HORIZONTAL TAIL PARAMETERS                         M50O62
C                                                                       M50O62
      CALL SECO(HTIN(1),CAMBER,ATYPE)                                   M50O62
      GO TO 1010                                                        M50O62
 1130 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- SET CALCULATED VERTICAL TAIL PARAMETERS                           M50O62
C                                                                       M50O62
      CALL SECO(VTIN(1),CAMBER,ATYPE)                                   M50O62
      GO TO 1010                                                        M50O62
 1140 CONTINUE                                                          M50O62
C                                                                       M50O62
C---- SET CALCULATED VENTRAL FIN PARAMETERS                             M50O62
C                                                                       M50O62
      CALL SECO(VFIN(1),CAMBER,ATYPE)                                   M50O62
      GO TO 1010                                                        M50O62
 1150 CONTINUE                                                          M50O62
 1160 FORMAT(68H0*** WARNING *** WING NOT STRAIGHT TAPERED. UNIFORM SECTM50O62
     1ION ASSUMED.)                                                     M50O62
 1170 FORMAT(68H0*** WARNING *** H.T. NOT STRAIGHT TAPERED. UNIFORM SECTM50O62
     1ION ASSUMED.)                                                     M50O62
 1180 FORMAT(68H0*** WARNING *** V.T. NOT STRAIGHT TAPERED. UNIFORM SECTM50O62
     1ION ASSUMED.)                                                     M50O62
 1190 FORMAT(68H0*** WARNING *** VFIN NOT STRAIGHT TAPERED. UNIFORM SECTM50O62
     1ION ASSUMED.)                                                     M50O62
      RETURN                                                            M50O62
      END                                                               M50O62
