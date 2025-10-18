      SUBROUTINE XNAM8(IOP)
C
C***  READ OR WRITE NAMELIST HTSCHR
C
      COMMON /HTI/ A6A(154)
      COMMON /IBH/ A6D(363)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER SCHR
      LOGICAL EOF
C
      DIMENSION LENHS(34),LDMHS(34),SCHR(160),LOCHS(34),HS(306)
      DIMENSION NLNAME(6)
      DATA NLNAME / 4HH   ,4HT   ,4HS   ,4HC   ,4HH   ,4HR   /
      DATA LENHS / 4,6,4,3,6,6,5,6,3,5,4,4,4,5,3,4,5,6,5,5,4,6,3*3,
     1            6,4,5,6,6,4,5,6,6 /
      DATA LDMHS / 5*1,20,20,-1,16*1,20,1,1,5*50,1,1 /
      DATA LOCHS / 16,17,18,19,20,21,41,64,61,66,67,62,63,65,61,
     1            67,70,71,69,69,92,68,93,94,72,306,155,
     2            156,206,256,206,256,10,10 /
      DATA SCHR /   4HT   ,4HO   ,4HV   ,4HC   ,4HD   ,4HE   ,
     14HL   ,4HT   ,4HA   ,4HY   ,4HX   ,4HO   ,4HV   ,4HC   ,
     24HC   ,4HL   ,4HI   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,
     34HI   ,4HC   ,4HL   ,4HA   ,4HL   ,4HP   ,4HA   ,4HC   ,
     44HL   ,4HM   ,4HA   ,4HX   ,4HC   ,4HA   ,4HM   ,4HB   ,
     54HE   ,4HR   ,4HC   ,4HM   ,4H0   ,4HX   ,4HO   ,4HV   ,
     64HC   ,4HO   ,4HC   ,4HM   ,4H0   ,4HT   ,4HL   ,4HE   ,
     74HR   ,4HI   ,4HL   ,4HE   ,4HR   ,4HO   ,4HT   ,4HO   ,
     84HV   ,4HC   ,4HO   ,4HC   ,4HM   ,4HO   ,4HC   ,4HM   ,
     94HO   ,4HT   ,4HT   ,4HC   ,4HE   ,4HF   ,4HF   ,4HK   ,
     A4HS   ,4HH   ,4HA   ,4HR   ,4HP   ,4HC   ,4HL   ,4HA   ,
     B4HM   ,4H0   ,4HC   ,4HL   ,4HA   ,4HM   ,4HO   ,4HA   ,
     C4HR   ,4HC   ,4HL   ,4HC   ,4HL   ,4HM   ,4HA   ,4HX   ,
     D4HL   ,4HY   ,4HC   ,4HM   ,4HC   ,4HL   ,4HD   ,4HX   ,
     E4HA   ,4HC   ,4HT   ,4HY   ,4HP   ,4HE   ,4HI   ,4HN   ,
     F4HN   ,4HP   ,4HT   ,4HS   ,4HX   ,4HC   ,4HO   ,4HR   ,
     G4HD   ,4HY   ,4HU   ,4HP   ,4HP   ,4HE   ,4HR   ,4HY   ,4HL   ,
     H4HO   ,4HW   ,4HE   ,4HR   ,4HM   ,4HE   ,4HA   ,4HN   ,4HT   ,
     I4HH   ,4HI   ,4HC   ,4HK   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,
     J4HO   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,4H0   /
C
      DO 1000 I=1,154
         HS(I)=A6A(I)
 1000 CONTINUE
      HS(306)=A6D(132)
      DO 1010 I=155,305
         HS(I)=A6D(I+58)
 1010 CONTINUE
C
C**   IOP EQUAL ZERO READ NAMELIST HTSCHR
C**   IOP EQUAL ONE WRITE NAMELIST HTSCHR
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,SCHR,160,LENHS,34,LDMHS,HS,306,LOCHS,
     2            EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,SCHR,160,LENHS,34,LDMHS,HS,306,LOCHS)
C
C     PUT THE VARIABLES IN THE APPROPRIATE COMMON BLOCK
C
      DO 1020 I=1,154
         A6A(I)=HS(I)
 1020 CONTINUE
      A6D(132)=HS(306)
      DO 1030 I=155,305
         A6D(I+58)=HS(I)
 1030 CONTINUE
C
      RETURN
      END
