      SUBROUTINE XNAM21(IOP)
C
C***  READ OR WRITE NAMELIST VFSCHR
C
      COMMON /VTI/ A9(316)
      COMMON /IVF/ A9B(363)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER SCHR
      LOGICAL EOF
C
      DIMENSION LENVFS(34),LDMVFS(34),SCHR(160),LOCVFS(34),VFS(476)
      DIMENSION NLNAME(6)
      DATA LENVFS / 4,6,4,3,6,6,5,6,3,5,4,4,4,5,3,4,5,6,5,5,4,6,3*3,
     1            6,4,5,6,6,4,5,6,6 /
      DATA LDMVFS / 5*1,20,20,-1,16*1,20,1,1,5*50,1,1 /
      DATA NLNAME / 4HV   ,4HF   ,4HS   ,4HC   ,4HH   ,4HR   /
C
      DATA LOCVFS / 178,179,180,181,182,183,203,226,223,228,229,224,
     1            225,227,223,229,232,233,231,231,254,230,0,0,234,476,
     2            325,326,376,426,376,426,10,10 /
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
     H4HO   ,4HW   ,4HE   ,4HR   ,4HM   ,4HE   ,4HA   ,4HN   ,
     I4HT   ,4HH   ,4HI   ,4HC   ,4HK   ,4HA   ,4HL   ,4HP   ,
     J4HH   ,4HA   ,4HO   ,4HA   ,4HL   ,4HP   ,4HH   ,4HA   ,4H0   /
C
      DO 1000 I=1,324
         VFS(I)=A9(I)
 1000 CONTINUE
      VFS(476)=A9B(132)
      DO 1010 I=325,475
         VFS(I)=A9B(I-112)
 1010 CONTINUE
C
C**   IOP EQUAL ZERO READ NAMELIST VFSCHR
C**   IOP EQUAL ONE WRITE NAMELIST VFSCHR
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,SCHR,160,LENVFS,34,LDMVFS,VFS,476,
     2            LOCVFS,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,SCHR,160,LENVFS,34,LDMVFS,VFS,476,
     2            LOCVFS)
C
C     PUT THE VARIABLES IN THE APPROPRIATE COMMON BLOCK
C
      DO 1020 I=1,324
         A9(I)=VFS(I)
 1020 CONTINUE
      A9B(132)=VFS(476)
      DO 1030 I=325,475
         A9B(I-112)=VFS(I)
 1030 CONTINUE
C
      RETURN
      END
