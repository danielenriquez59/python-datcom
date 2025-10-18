      SUBROUTINE XNAM16(IOP)
C
      COMMON /FLAPIN/ AI(137)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER SYMFLP
      LOGICAL EOF
C
      DIMENSION LENI(27),LDMI(27),SYMFLP(139),LOCI(27),AIA(117)
      DIMENSION NLNAME(6)
C
      DATA NLNAME / 4HS   ,4HY   ,4HM   ,4HF   ,4HL   ,4HP   /
      DATA LENI / 6*6,3*5,4,2,2,4,5,6,6,4,4,3,8*6 /
      DATA LDMI / 13*1,5*9,1,9,1,4*9,2*1 /
      DATA LOCI / 12,13,14,15,16,61,11,17,62,18,59,60,18,1,39,
     1            49,19,29,63,64,74,75,85,95,105,115,116 /
      DATA SYMFLP /   4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,4HI   ,
     1  4HC   ,4HH   ,4HR   ,4HD   ,4HF   ,4HO   ,4HS   ,4HP   ,
     2  4HA   ,4HN   ,4HF   ,4HI   ,4HS   ,4HP   ,4HA   ,4HN   ,
     3  4HF   ,4HO   ,4HN   ,4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,
     4  4HP   ,4HH   ,4HE   ,4HT   ,4HE   ,4HP   ,4HP   ,4HH   ,
     5  4HE   ,4HT   ,4HE   ,4HF   ,4HT   ,4HY   ,4HP   ,4HE   ,
     6  4HN   ,4HT   ,4HY   ,4HP   ,4HE   ,4HS   ,4HC   ,4HH   ,
     7  4HA   ,4HC   ,4HB   ,4HT   ,4HC   ,4HS   ,4HC   ,4HH   ,
     8  4HD   ,4HD   ,4HE   ,4HL   ,4HT   ,4HA   ,4HC   ,4HP   ,
     9  4HR   ,4HM   ,4HE   ,4HI   ,4HC   ,4HP   ,4HR   ,4HM   ,
     A  4HE   ,4HO   ,4HS   ,4HC   ,4HL   ,4HD   ,4HS   ,4HC   ,
     B  4HM   ,4HD   ,4HC   ,4HM   ,4HU   ,4HD   ,4HE   ,4HL   ,
     C  4HJ   ,4HE   ,4HT   ,4HJ   ,4HE   ,4HT   ,4HF   ,4HL   ,
     D  4HP   ,4HE   ,4HF   ,4HF   ,4HJ   ,4HE   ,4HT   ,4HC   ,
     E  4HA   ,4HP   ,4HI   ,4HN   ,4HB   ,4HC   ,4HA   ,4HP   ,
     F  4HO   ,4HU   ,4HT   ,4HD   ,4HO   ,4HB   ,4HD   ,4HE   ,
     G  4HF   ,4HD   ,4HO   ,4HB   ,4HC   ,4HI   ,4HN   ,4HD   ,
     H  4HO   ,4HB   ,4HC   ,4HO   ,4HT   /
C
      DO 1000 I=1,116
         AIA(I)=AI(I)
 1000 CONTINUE
C
C**   IF IOP EQUAL ZERO READ NAMELIST SYMFLP
C**   IF IOP EQUAL ONE WRITE NAMELIST SYMFLP
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,SYMFLP,139,LENI,27,LDMI,AIA,117,
     2            LOCI,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,SYMFLP,139,LENI,27,LDMI,AIA,117,LOCI)
C
      DO 1010 I=1,116
         AI(I)=AIA(I)
 1010 CONTINUE
C
      RETURN
      END
