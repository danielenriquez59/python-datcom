      SUBROUTINE EQSPCE(X,R,P,S,NP,XE,RE,PE,SE,NE,DSEDX)
C
C***  TRANSFORMS 4-DIMENSIONAL ARRAY SO THAT THE 3 INDEPENDENT ARRAYS
C     ARE EQUALLY SPACED
C
      DIMENSION ROUTID(2)
      DIMENSION CR(6),CP(6),CS(6)
      DIMENSION X(1),R(1),P(1),S(1),XE(1),RE(1),PE(1),SE(1),
     1          DSEDX(20)
      DIMENSION VAR(4),LGH(4)
      DATA ROUTID/4HEQSP,4HCE  /
C
      IF(NP.EQ.1)GO TO 1020
      IR=0
      IP=0
      IS=0
C
C     ----EQUAL SPACE ON NEW X S
C
      FNE=NE-1
      XIN=(X(NP)-X(1))/FNE
      XE(1)=X(1)
      RE(1)=R(1)
      PE(1)=P(1)
      SE(1)=S(1)
      XE(NE)=X(NP)
      RE(NE) =R(NP)
      PE(NE)=P(NP)
      SE(NE)=S(NP)
      NN=NE-1
      DO 1000 I=2,NN
         XE(I)=XE(I-1)+XIN
         VAR(1)=XE(I)
         LGH(1)=NP
         CALL INTERX(1,X,VAR,LGH,R,RE(I),NP,NP,0,0,0,0,0,0,0,0,
     1               4HRE  ,1,ROUTID)
         CALL INTERX(1,X,VAR,LGH,P,PE(I),NP,NP,0,0,0,0,0,0,0,0,
     1               4HPE  ,1,ROUTID)
         CALL INTERX(1,X,VAR,LGH,S,SE(I),NP,NP,0,0,0,0,0,0,0,0,
     1               4HSE  ,1,ROUTID)
 1000 CONTINUE
      DO 1010 I=1,NE
         CALL TBFUNX(XE(I),DUMM,DSEDX(I),NE,XE,SE,CS,IS,MI,NG,0,0,
     1               4HDSED,1,ROUTID)
 1010 CONTINUE
      RETURN
 1020 DO 1030 I=1,NE
         XE(I)=X(1)
         RE(I)=R(1)
         PE(I)=P(1)
         SE(I)=S(1)
 1030 DSEDX(I)=0.0
      RETURN
      END
