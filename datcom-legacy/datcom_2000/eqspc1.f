      SUBROUTINE EQSPC1(X,S,NP,XE,SE,NE,DSEDX)
C
C     ----COMPUTE S AND DS/DX ON EQUAL SPACED X
C
      DIMENSION ROUTID(2)
      DIMENSION CS(6),X(1),S(1),XE(1),SE(1),DSEDX(1)
      DIMENSION VAR(4),LGH(4)
      DATA ROUTID/4HEQSP,4HC1  /
C
      IF(NP.EQ.1)GO TO 1020
      IS=0
      FNE=NE-1
      XIN=(X(NP)-X(1))/FNE
      XE(1)=X(1)
      SE(1)=S(1)
      XE(NE)=X(NP)
      SE(NE)=S(NP)
      NN=NE-1
      DO 1000 I=2,NN
         XE(I)=XE(I-1)+XIN
         VAR(1)=XE(I)
         LGH(1)=NP
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
         SE(I)=S(1)
 1030 DSEDX(I)=0.0
      RETURN
      END
