      SUBROUTINE XNAM5(IOP)
C
C***  READ OR WRITE NAMELIST WGSCHR
C
      COMMON /WINGI/  A6A(101)
      COMMON /IBW/    A6D(363)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER SCHR
      LOGICAL EOF
C
      DIMENSION LEN6(36),LDM6(36),SCHR(170),LOC6(36),A6(253)
      DIMENSION NLNAME(6)
      DATA NLNAME / 4HW   ,4HG   ,4HS   ,4HC   ,4HH   ,4HR   /
      DATA LEN6 / 4,6,4,3,6,6,5,6,3,5,4,4,4,5,3,4,5,6,6,5,5,5,4,
     1            3,5,3,3,6,4,5,6,6,4,5,6,6 /
      DATA LDM6 / 5*1,20,20,-1,11*1,6,3*1,20,5*1,5*50,1,1 /
      DATA LOC6 / 16,17,18,19,20,21,41,64,61,66,67,62,63,65,61,
     1            67,70,71,68,95,69,69,92,72,101,93,94,253,102,
     2            103,153,203,153,203,10,10 /
      DATA SCHR /   4HT   ,4HO   ,4HV   ,4HC   ,4HD   ,4HE   ,4HL   ,
     14HT   ,4HA   ,4HY   ,4HX   ,4HO   ,4HV   ,4HC   ,4HC   ,4HL   ,
     24HI   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,4HI   ,4HC   ,4HL   ,
     34HA   ,4HL   ,4HP   ,4HA   ,4HC   ,4HL   ,4HM   ,4HA   ,4HX   ,
     44HC   ,4HA   ,4HM   ,4HB   ,4HE   ,4HR   ,4HC   ,4HM   ,4H0   ,
     54HX   ,4HO   ,4HV   ,4HC   ,4HO   ,4HC   ,4HM   ,4H0   ,4HT   ,
     64HL   ,4HE   ,4HR   ,4HI   ,4HL   ,4HE   ,4HR   ,4HO   ,4HT   ,
     74HO   ,4HV   ,4HC   ,4HO   ,4HC   ,4HM   ,4HO   ,4HC   ,4HM   ,
     84HO   ,4HT   ,4HT   ,4HC   ,4HE   ,4HF   ,4HF   ,4HK   ,4HS   ,
     94HH   ,4HA   ,4HR   ,4HP   ,4HC   ,4HL   ,4HM   ,4HA   ,4HX   ,
     A4HL   ,4HS   ,4HL   ,4HO   ,4HP   ,4HE   ,4HC   ,4HL   ,4HA   ,
     B4HM   ,4HO   ,4HC   ,4HL   ,4HA   ,4HM   ,4H0   ,4HA   ,4HR   ,
     C4HC   ,4HL   ,4HX   ,4HA   ,4HC   ,4HD   ,4HW   ,4HA   ,4HS   ,
     D4HH   ,4HY   ,4HC   ,4HM   ,4HC   ,4HL   ,4HD   ,4HT   ,4HY   ,
     E4HP   ,4HE   ,4HI   ,4HN   ,4HN   ,4HP   ,4HT   ,4HS   ,4HX   ,
     F4HC   ,4HO   ,4HR   ,4HD   ,4HY   ,4HU   ,4HP   ,4HP   ,4HE   ,
     G4HR   ,4HY   ,4HL   ,4HO   ,4HW   ,4HE   ,4HR   ,4HM   ,4HE   ,
     H4HA   ,4HN   ,4HT   ,4HH   ,4HI   ,4HC   ,4HK   ,4HA   ,4HL   ,
     I4HP   ,4HH   ,4HA   ,4HO   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,
     J4H0   /
C
C***  WRITE NAMELIST WGSCHR
C
      DO 1000 I=1,101
         A6(I)=A6A(I)
 1000 CONTINUE
      A6(253)=A6D(132)
      DO 1010 I=102,252
         A6(I)=A6D(I+111)
 1010 CONTINUE
C
C**   IOP EQUAL ZERO READ NAMELIST WGSCHR
C**   IOP EQUAL ONE WRITE NAMELIST WGSCHR
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,SCHR,170,LEN6,36,LDM6,A6,253,LOC6,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,SCHR,170,LEN6,36,LDM6,A6,253,LOC6)
C
C     PUT THE VARIABLES IN THE APPROPRIATE COMMON BLOCK
C
      DO 1020 I=1,101
         A6A(I)=A6(I)
 1020 CONTINUE
      A6D(132)=A6(253)
      DO 1030 I=102,252
         A6D(I+111)=A6(I)
 1030 CONTINUE
C
      RETURN
      END
