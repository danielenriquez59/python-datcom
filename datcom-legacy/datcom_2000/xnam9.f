      SUBROUTINE XNAM9(IOP)
C
C***  READ OR WRITE NAMELIST VTPLNF
C
      COMMON /VTI/    A9(162)
      COMMON /VTDATA/ AVT(195)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      EQUIVALENCE (TYPE,A9(15))
C
      INTEGER PLNF
      LOGICAL EOF
C
      DIMENSION NLNAME(6),XTYPE(4),VT(75),IEQU(18)
      DIMENSION LENVP(18),PLNF(92),LOCVP(18),LDMVP(18)
C
      DATA XTYPE /4HSTRA,4HDDUB,4HCRAN,4HCURV/
      DATA NLNAME /4HV   ,4HT   ,4HP   ,4HL   ,4HN   ,4HF   /
      DATA LENVP / 6,6,5,4,6,5,5,6,5,4*6,4,5,4,3,4/
      DATA LDMVP / 15*1,3*20 /
      DATA LOCVP / 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,36,56 /
      DATA IEQU / 1,2,3,4,5,6,7,9,10,11,12,13,14,15,8,95,115,135 /
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
         VT(I)=A9(IE)
 1000 CONTINUE
      DO 1020 I=16,18
         IE=IEQU(I)
         IF(I .EQ. 16)K=-1
         IF(I .EQ. 17)K=18
         IF(I .EQ. 18)K=37
         DO 1010 J=1,20
            VT(I+J+K)=A9(IE+J-1)
 1010    CONTINUE
 1020 CONTINUE
C
C**   IF IOP EQUAL ZERO READ NAMELIST VTPLNF
C**   IF IOP EQUAL ONE WRITE NAMELIST VTPLNF
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,PLNF,92,LENVP,18,LDMVP,VT,75,
     2            LOCVP,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,PLNF,92,LENVP,18,LDMVP,VT,75,LOCVP)
C
      DO 1030 I=1,15
         IE=IEQU(I)
C
C  THIS IF CHECK IS BECAUSE VARIABLE 10 IS A DUMMY VARIABLE.
C
         IF(IE .EQ. 10)GO TO 1030
         A9(IE)=VT(I)
 1030 CONTINUE
      DO 1050 I=16,18
         IE=IEQU(I)
         IF(I .EQ. 16)K=-1
         IF(I .EQ. 17)K=18
         IF(I .EQ. 18)K=37
         DO 1040 J=1,20
            A9(IE+J-1)=VT(I+J+K)
 1040    CONTINUE
 1050 CONTINUE
C
      AVT(106)=A9(7)
      AVT(112)=A9(8)
      AVT(138)=0.0
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
