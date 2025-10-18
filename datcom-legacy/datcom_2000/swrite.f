      SUBROUTINE SWRITE(LAST,FOR,NW,ICONT,NRPEAT,A,B,C,D,E,F,G,H,O,P,Q
     1                 ,R,S,T,NDMF,NAF)
C
C***  CONTROLS NUMERIC OUTPUT FOR PRINTOUT, WRITES BLANKS,NA OR NDM
C
      DIMENSION A(1),B(1),C(1),D(1),E(1),F(1),G(1),H(1),O(1),P(1),Q(1)
     1         ,R(1),S(1),T(1),ICONT(1)
      DIMENSION Z(14),FOR(NW),ICPAT(14),FORM(60),FORSUB(21),IPPAT(14)
      CHARACTER*4 FOR,FORM,FORSUB
      REAL NDM,NA
      LOGICAL NDMF,NAF
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      DATA BLANK/4H    /,NDM/4HNDM /,NA/4HNA  /
      DATA FORSUB/
     1 ',2X ',',A4,','1X  ',',2X ',',A4,','2X  ',',3X ',',A4,','2X  ',
     2 ',3X ',',A4,','3X  ',',4X ',',A4,','3X  ',',4X ',',A4,','4X  ',
     3 ',5X ',',A4,','4X  '/
C
C                    STORE BASIC FORMAT
C
      DO 1000 I=1,NW
 1000 FORM(I)=FOR(I)
C
C              STORE ZEROS IN PREVIOUS PATTERN
C
      DO 1010 I=1,14
 1010 IPPAT(I)=0
C
C              STORE L-TH ELEMENT OF EACH ARRAY IN A---------------
C                   HORIZONTAL ROW ARRAY
C
      DO 1230 L=1,NRPEAT
          GO TO (1150,1140,1130,1120,1110,1100,1090,1080,1070,1060,1050,
     11040,1030,1020          ),LAST
 1020     Z(14)=T(L)
 1030     Z(13)=S(L)
 1040     Z(12)=R(L)
 1050     Z(11)=Q(L)
 1060     Z(10)=P(L)
 1070     Z(9)=O(L)
 1080     Z(8)=H(L)
 1090     Z(7)=G(L)
 1100     Z(6)=F(L)
 1110     Z(5)=E(L)
 1120     Z(4)=D(L)
 1130     Z(3)=C(L)
 1140     Z(2)=B(L)
 1150     Z(1)=A(L)
C
C               CONSTRUCT FORMAT FOR HORIZONTAL ROW PRINT OUT
C
          DO 1220 I=1,LAST
              IF(Z(I).NE.UNUSED) GO TO 1160
                  ICPAT(I)=1
                  Z(I)=NDM
                  NDMF=.TRUE.
                  GO TO 1190
 1160         IF(Z(I).NE.-UNUSED) GO TO 1170
                  ICPAT(I)=1
                  Z(I)=BLANK
                  GO TO 1190
 1170         IF(Z(I).NE.2*UNUSED) GO TO 1180
                  ICPAT(I)=1
                  Z(I)=NA
                  NAF=.TRUE.
                  GO TO 1190
 1180         CONTINUE
                  ICPAT(I)=0
 1190         CONTINUE
C
C                COMPARE CURRENT PATTERN WITH PREVIOUS PATTERN
C
              IF(ICPAT(I).EQ.IPPAT(I)) GO TO 1220
              K=4*I-1
              IF(IPPAT(I).NE.0) GO TO 1200
                  J=ICONT(I)-6
                  J2=3*J-2
                  FORM(K)=FORSUB(J2)
                  FORM(K+1)=FORSUB(J2+1)
                  FORM(K+2)=FORSUB(J2+2)
                  GO TO 1210
 1200         CONTINUE
                  FORM(K)=FOR(K)
                  FORM(K+1)=FOR(K+1)
                  FORM(K+2)=FOR(K+2)
 1210         IPPAT(I)=ICPAT(I)
 1220     CONTINUE
C
C                         PRINT L-TH HORIZONTAL ROW
C
          WRITE(6,FORM)(Z(J),J=1,LAST)
 1230 CONTINUE
      RETURN
      END
