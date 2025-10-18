      SUBROUTINE PTINT2(XXP,YYP,A,AMU,X,Y,XI,YI,NX,INDXUI,INDXLI,INXUIE,
     1                  INXLIE,ABORT)
C
C***  CALCULATES BOUNDRIES OF MACH LINES ON BODIES
C
      LOGICAL ABORT
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      DIMENSION XXP(1),YYP(1),A(1),X(42),Y(42),XDIF(42),YDIF(42),XI(42),
     1          YI(42),XP(42),YP(42)
C
      ABORT=.FALSE.
      TAMU=SIN(AMU)/COS(AMU)
      TMAMU=-TAMU
      INDXUI=0
      INDXLI=0
      INXUIE=0
      INXLIE=0
      DX = -0.0001*XXP(NX)
      DO 1000 KL=1,NX
         XP(KL)=XXP(KL)
 1000 YP(KL)=YYP(KL)
      XP(NX+1)=XP(NX)
      YP(NX+1)=-YP(NX)
      NX1=NX+1
      NXM1=NX-1
      NX22M1=NX+21
      DO 1010 I=1,NX1
         X(I)=(XP(I)-A(1))*COS(-A(3)/RAD)+(YP(I)-A(2))*SIN(-A(3)/RAD)
         Y(I)=(YP(I)-A(2))*COS(-A(3)/RAD)-(XP(I)-A(1))*SIN(-A(3)/RAD)
         X(I+21)=(XP(I)-A(1))*COS(-A(3)/RAD)
     1          +(-YP(I)-A(2))*SIN(-A(3)/RAD)
         Y(I+21)=(-YP(I)-A(2))*COS(-A(3)/RAD)
     1          -(XP(I)-A(1))*SIN(-A(3)/RAD)
         L=I+21
 1010 CONTINUE
      IF(X(NX) .LT. 0.0 .AND. X(NX) .GT. DX) X(NX) = 0.0
      IF(X(NX22M1) .LT. 0.0 .AND. X(NX22M1) .GT. DX) X(NX22M1) = 0.0
      DO 1020 J=1,NX
         IF(0.0.LE.X(J)) GO TO 1030
 1020 CONTINUE
 1030 INDEXU=J-1
      DO 1040 J=1,NX
         IF(0.0.LE.X(J+21)) GO TO 1050
 1040 CONTINUE
 1050 INDEXL=J+20
      DO 1080 K=INDEXU,NX
         XDIF(K)=X(K+1)-X(K)
         YDIF(K)=Y(K+1)-Y(K)
         IF(XDIF(K) .EQ. 0.0 .AND. YDIF(K) .EQ. 0.0) GO TO 1080
         XI(K)=(Y(K)*XDIF(K)-X(K)*YDIF(K))/(XDIF(K)*TAMU-YDIF(K))
         YI(K)=XI(K)*TAMU
         IF(K.NE.INDEXU) GO TO 1070
         IF(X(NX).LT.0.0) GO TO 1060
         YIN=(YDIF(K)/XDIF(K))*(-X(K))+Y(K)
         IF(YIN.GE.0.0) GO TO 1070
 1060    ABORT=.TRUE.
         RETURN
 1070    CONTINUE
         IF(XI(K).GE.X(K).AND.XI(K).LE.X(K+1)) GO TO 1090
         IF(K.NE.NXM1) GO TO 1080
         IF(XI(K).GT.X(K+1)) INXUIE=NX-1
 1080 CONTINUE
      K=NX
 1090 INDXUI=K
      DO 1120 K=INDEXL,NX22M1
         XDIF(K)=X(K+1)-X(K)
         YDIF(K)=Y(K+1)-Y(K)
         IF(XDIF(K) .EQ. 0.0 .AND. YDIF(K) .EQ. 0.0) GO TO 1120
         XI(K)=(Y(K)*XDIF(K)-X(K)*YDIF(K))/(XDIF(K)*TMAMU-YDIF(K))
         YI(K)=-XI(K)*TAMU
         IF(K.NE.INDEXL) GO TO 1110
         IF(X(NX22M1).LT.0.0) GO TO 1100
         YIN=(YDIF(K)/XDIF(K))*(-X(K))+Y(K)
         IF(YIN.LE.0.0) GO TO 1110
 1100    ABORT=.TRUE.
         RETURN
 1110    CONTINUE
         IF(XI(K).GE.X(K).AND.XI(K).LE.X(K+1)) GO TO 1130
         IF(K.NE.(NX22M1-1)) GO TO 1120
         IF(XI(K).GT.X(K+1)) INXLIE=NX22M1-1
 1120 CONTINUE
      K=NX22M1
 1130 INDXLI=K
      RETURN
      END
