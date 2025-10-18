      SUBROUTINE XNAM7(IOP)
C
C***  READ OR WRITE NAMELIST HTPLNF
C
      COMMON /HTI/    A5B(154)
      COMMON /HTDATA/ AHT(195)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      EQUIVALENCE (TYPE,A5B(15))
C
      INTEGER PLNF
      LOGICAL EOF
C
      DIMENSION PLNF(92),LENHP(18),LDMHP(18),A(75)
      DIMENSION NLNAME(6),XTYPE(4),IEQU(18),LOCHP(18)
C
      DATA NLNAME /4HH   ,4HT   ,4HP   ,4HL   ,4HN   ,4HF   /
      DATA PLNF /4HC   ,4HH   ,4HR   ,4HD   ,4HB   ,4HP   ,
     1    4HC   ,4HH   ,4HR   ,4HD   ,4HR   ,4HC   ,4HH   ,
     2    4HR   ,4HD   ,4HT   ,4HP   ,4HC   ,4HH   ,4HS   ,
     3    4HT   ,4HA   ,4HT   ,4HS   ,4HS   ,4HP   ,4HN   ,
     4    4HS   ,4HS   ,4HP   ,4HN   ,4HE   ,4HS   ,4HS   ,
     5    4HP   ,4HN   ,4HO   ,4HP   ,4HS   ,4HA   ,4HV   ,
     6    4HS   ,4HI   ,4HS   ,4HA   ,4HV   ,4HS   ,4HO   ,
     7    4HS   ,4HW   ,4HA   ,4HF   ,4HP   ,4HT   ,4HW   ,
     8    4HI   ,4HS   ,4HT   ,4HA   ,4HT   ,4HY   ,4HP   ,
     9    4HE   ,4HS   ,4HS   ,4HP   ,4HN   ,4HD   ,4HD   ,
     A    4HD   ,4HH   ,4HD   ,4HA   ,4HD   ,4HI   ,4HD   ,
     B    4HH   ,4HD   ,4HA   ,4HD   ,4HO   ,4HR   ,4HL   ,
     C    4HP   ,4HH   ,4HS   ,4HH   ,4HB   ,4HS   ,4HE   ,
     D    4HX   ,4HT   /
      DATA LENHP /6,5,6,6,4,5,6,3*5,6,4,3*6,4,3,4/
      DATA LDMHP /15*1,3*20/
      DATA LOCHP /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,36,56/
      DATA XTYPE / 4HSTRA,4HDDUB,4HCRAN,4HCURV /
      DATA IEQU /5,6,1,9,4,3,2,7,8,10,11,15,12,13,14,95,115,135/
C
      DO 1000 I=1,15
         IE=IEQU(I)
         A(I)=A5B(IE)
 1000 CONTINUE
      DO 1020 I=16,18
         IE=IEQU(I)
         IF(I .EQ. 16)K=-1
         IF(I .EQ. 17)K=18
         IF(I .EQ. 18)K=37
         DO 1010 J=1,20
            A(I+J+K)=A5B(IE+J-1)
 1010    CONTINUE
 1020 CONTINUE
C
C**   IOP EQUAL ZERO READ NAMELIST HTPLNF
C**   IOP EQUAL ONE WRITE NAMELIST HTPLNF
C
      IF(IOP .EQ. 0)
     1    CALL NAMER(KAND,9,NLNAME,6,PLNF,92,LENHP,18,LDMHP,A,75,
     2              LOCHP,EOF)
      IF(IOP .EQ. 1)
     1    CALL NAMEW(KAND,6,NLNAME,6,PLNF,92,LENHP,18,LDMHP,A,75,LOCHP)
C
      DO 1030 I=1,15
         IE=IEQU(I)
C
C  THIS IF CHECK IS BECAUSE VARIABLE 10 IS A DUMMY VARIBLE.  ELEMENT 10
C   READ IN NAMELIST HTSCHR
         IF(IE .EQ. 10)GO TO 1030
         A5B(IE)=A(I)
 1030 CONTINUE
      DO 1050 I=16,18
         IE=IEQU(I)
         IF(I .EQ. 16)K=-1
         IF(I .EQ. 17)K=18
         IF(I .EQ. 18)K=37
         DO 1040 J=1,20
            A5B(IE+J-1)=A(I+J+K)
 1040    CONTINUE
 1050 CONTINUE
C
      AHT(106)=A5B(7)
      AHT(112)=A5B(8)
      AHT(138)=0.0
C---  CHECK FOR TYPE SET IN NAMELIST
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
