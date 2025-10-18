      SUBROUTINE XNAM20(IOP)
C
C***  READ OR WRITE NAMELIST VFPLNF
C
      COMMON /VTI/    A5C(316)
      COMMON /VTDATA/ AVT(195),AVF(195)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      EQUIVALENCE (TYPE,A5C(177))
C
      INTEGER PLNF
      LOGICAL EOF
C
      DIMENSION NLNAME(6),XTYPE(4),VFP(75),IEQU(18)
      DIMENSION LENVFP(18),PLNF(92),LOCVFP(18),LDMVFP(18)
C
      DATA XTYPE /4HSTRA,4HDDUB,4HCRAN,4HCURV/
      DATA NLNAME /4HV   ,4HF   ,4HP   ,4HL   ,4HN   ,4HF   /
      DATA LENVFP / 6,6,5,4,6,5,5,6,5,4*6,4,5,4,3,4/
      DATA LDMVFP / 15*1,3*20 /
      DATA LOCVFP / 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,36,56 /
      DATA IEQU / 163,164,165,166,167,168,169,171,172,173,174,175,
     1            176,177,170,257,277,297 /
      DATA PLNF /4HC   ,4HH   ,4HR   ,4HD   ,4HT   ,4HP   ,4HS   ,
     1    4HS   ,4HP   ,4HN   ,4HO   ,4HP   ,4HS   ,4HS   ,4HP   ,
     2    4HN   ,4HE   ,4HS   ,4HS   ,4HP   ,4HN   ,4HC   ,
     3    4HH   ,4HR   ,4HD   ,4HB   ,4HP   ,4HC   ,4HH   ,4HR   ,
     4    4HD   ,4HR   ,4HS   ,4HA   ,4HV   ,4HS   ,4HI   ,4HC   ,
     5    4HH   ,4HS   ,4HT   ,4HA   ,4HT   ,4HS   ,4HW   ,4HA   ,
     6    4HF   ,4HP   ,4HT   ,4HW   ,4HI   ,4HS   ,4HT   ,4HA   ,
     7    4HS   ,4HS   ,4HP   ,4HN   ,4HD   ,4HD   ,4HD   ,4HH   ,
     8    4HD   ,4HA   ,4HD   ,4HI   ,4HD   ,4HH   ,4HD   ,4HA   ,
     9    4HD   ,4HO   ,4HT   ,4HY   ,4HP   ,4HE   ,4HS   ,4HA   ,
     A    4HV   ,4HS   ,4HO   ,4HS   ,4HV   ,4HW   ,4HB   ,4HS   ,
     B    4HV   ,4HB   ,4HS   ,4HV   ,4HH   ,4HB   /
C
      DO 1000 I=1,15
         IE=IEQU(I)
         VFP(I)=A5C(IE)
 1000 CONTINUE
      DO 1020 I=16,18
         IE=IEQU(I)
         IF(I .EQ. 16)K=-1
         IF(I .EQ. 17)K=18
         IF(I .EQ. 18)K=37
         DO 1010 J=1,20
            VFP(I+J+K)=A5C(IE+J-1)
 1010    CONTINUE
 1020 CONTINUE
C
C**   IF IOP EQUAL ZERO READ NAMELIST VTPLNF
C**   IF IOP EQUAL ONE WRITE NAMELIST VTPLNF
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,PLNF,92,LENVFP,18,LDMVFP,VFP,75,
     2            LOCVFP,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,PLNF,92,LENVFP,18,LDMVFP,VFP,
     2            75,LOCVFP)
C
      DO 1030 I=1,15
         IE=IEQU(I)
         IF(IE .EQ. 172)GO TO 1030
         A5C(IE)=VFP(I)
 1030 CONTINUE
      DO 1050 I=16,18
         IE=IEQU(I)
         IF(I .EQ. 16)K=-1
         IF(I .EQ. 17)K=18
         IF(I .EQ. 18)K=37
         DO 1040 J=1,20
            A5C(IE+J-1)=VFP(I+J+K)
 1040    CONTINUE
 1050 CONTINUE
C
      AVF(106)=A5C(7)
      AVF(112)=A5C(8)
      AVF(138)=0.0
C --- CHECK IF TYPE SET IN NAMELIST
      DO 1060 I=1,4
         IF(TYPE .EQ. XTYPE(I))TYPE=I
 1060 CONTINUE
      IF(TYPE .LT. 5.)GO TO 1070
      WRITE(6,1080)TYPE
      TYPE=1.0
 1070 RETURN
 1080 FORMAT(42H ERROR TYPE CANNOT BE GREATER THAN 4 TYPE=,E12.5,
     19H SET TO 1)
      END
