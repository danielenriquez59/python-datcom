      SUBROUTINE XNAM10(IOP)
C
C***  READ OR WRITE NAMELIST VTSCHR
C
      COMMON /VTI/ A9(162)
      COMMON /IBV/ A9B(363)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER SCHR
      LOGICAL EOF
C
      DIMENSION LENVTS(34),LDMVTS(34),SCHR(160),LOCVTS(34),VTS(314)
      DIMENSION NLNAME(6)
      DATA NLNAME / 4HV   ,4HT   ,4HS   ,4HC   ,4HH   ,4HR   /
      DATA LENVTS / 4,6,4,3,6,6,5,6,3,5,4,4,4,5,3,4,5,6,5,5,4,6,3*3,
     1            6,4,5,6,6,4,5,6,6 /
      DATA LDMVTS / 5*1,20,20,-1,16*1,20,1,1,5*50,1,1 /
      DATA LOCVTS / 16,17,18,19,20,21,41,64,61,66,67,62,63,65,61,
     1            67,70,71,69,69,92,68,155,155,72,314,163,
     2            164,214,264,214,264,10,10 /
      DATA SCHR /     4HT   ,4HO   ,4HV   ,4HC   ,4HD   ,4HE   ,
     1 4HL   ,4HT   ,4HA   ,4HY   ,4HX   ,4HO   ,4HV   ,4HC   ,
     2 4HC   ,4HL   ,4HI   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,
     3 4HI   ,4HC   ,4HL   ,4HA   ,4HL   ,4HP   ,4HA   ,4HC   ,
     4 4HL   ,4HM   ,4HA   ,4HX   ,4HC   ,4HA   ,4HM   ,4HB   ,
     5 4HE   ,4HR   ,4HC   ,4HM   ,4H0   ,4HX   ,4HO   ,4HV   ,
     6 4HC   ,4HO   ,4HC   ,4HM   ,4H0   ,4HT   ,4HL   ,4HE   ,
     7 4HR   ,4HI   ,4HL   ,4HE   ,4HR   ,4HO   ,4HT   ,4HO   ,
     8 4HV   ,4HC   ,4HO   ,4HC   ,4HM   ,4HO   ,4HC   ,4HM   ,
     9 4HO   ,4HT   ,4HT   ,4HC   ,4HE   ,4HF   ,4HF   ,4HK   ,
     A 4HS   ,4HH   ,4HA   ,4HR   ,4HP   ,4HC   ,4HL   ,4HA   ,
     B 4HM   ,4H0   ,4HC   ,4HL   ,4HA   ,4HM   ,4HO   ,4HA   ,
     C 4HR   ,4HC   ,4HL   ,4HC   ,4HL   ,4HM   ,4HA   ,4HX   ,
     D 4HL   ,4HY   ,4HC   ,4HM   ,4HC   ,4HL   ,4HD   ,4HX   ,
     E 4HA   ,4HC   ,4HT   ,4HY   ,4HP   ,4HE   ,4HI   ,4HN   ,
     F 4HN   ,4HP   ,4HT   ,4HS   ,4HX   ,4HC   ,4HO   ,4HR   ,
     G 4HD   ,4HY   ,4HU   ,4HP   ,4HP   ,4HE   ,4HR   ,4HY   ,4HL   ,
     H 4HO   ,4HW   ,4HE   ,4HR   ,4HM   ,4HE   ,4HA   ,4HN   ,
     I 4HT   ,4HH   ,4HI   ,4HC   ,4HK   ,4HA   ,4HL   ,4HP   ,
     J 4HH   ,4HA   ,4HO   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,4H0   /
C
      DO 1000 I=1,162
         VTS(I)=A9(I)
 1000 CONTINUE
      VTS(314)=A9B(132)
      DO 1010 I=163,313
         VTS(I)=A9B(I+50)
 1010 CONTINUE
C
C**   IOP EQUAL ZERO READ NAMELIST VTSCHR
C**   IOP EQUAL ONE WRITE NAMELIST VTSCHR
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,SCHR,160,LENVTS,34,LDMVTS,VTS,314,
     2            LOCVTS,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,SCHR,160,LENVTS,34,LDMVTS,VTS,314,
     2            LOCVTS)
C
C     PUT THE VARIABLES IN THE APPROPRIATE COMMON BLOCK
C
      DO 1020 I=1,162
         IF(I .EQ. 155)GO TO 1020
         A9(I)=VTS(I)
 1020 CONTINUE
      A9B(132)=VTS(314)
      DO 1030 I=163,313
         A9B(I+50)=VTS(I)
 1030 CONTINUE
C
      RETURN
      END
