      SUBROUTINE XNAM4(IOP)
C
C***  READ OR WRITE NAMELIST WGPLNF
C
      COMMON /WINGI/  B(101)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
      COMMON /WINGD/  WNG(138)
C
      INTEGER PLNF
      LOGICAL EOF
C
      DIMENSION LEN5(15),LDM5(15),PLNF(81),LOC5(15),A(15),
     1          IEQU(15)
      DIMENSION NLNAME(6)
      DATA NLNAME / 4HW   ,4HG   ,4HP   ,4HL   ,4HN   ,4HF   /
      DATA LEN5 / 6,5,6,6,4,5,6,3*5,6,4,3*6/
      DATA LDM5 / 15*1/
      DATA LOC5 / 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
      DATA PLNF /     4HC   ,4HH   ,4HR   ,4HD   ,4HB   ,4HP   ,
     1  4HC   ,4HH   ,4HR   ,4HD   ,4HR   ,4HC   ,4HH   ,4HR   ,
     2  4HD   ,4HT   ,4HP   ,4HC   ,4HH   ,4HS   ,4HT   ,4HA   ,
     3  4HT   ,4HS   ,4HS   ,4HP   ,4HN   ,4HS   ,4HS   ,4HP   ,
     4  4HN   ,4HE   ,4HS   ,4HS   ,4HP   ,4HN   ,4HO   ,4HP   ,
     5  4HS   ,4HA   ,4HV   ,4HS   ,4HI   ,4HS   ,4HA   ,4HV   ,
     6  4HS   ,4HO   ,4HS   ,4HW   ,4HA   ,4HF   ,4HP   ,4HT   ,
     7  4HW   ,4HI   ,4HS   ,4HT   ,4HA   ,4HT   ,4HY   ,4HP   ,
     8  4HE   ,4HS   ,4HS   ,4HP   ,4HN   ,4HD   ,4HD   ,4HD   ,
     9  4HH   ,4HD   ,4HA   ,4HD   ,4HI   ,4HD   ,4HH   ,4HD   ,
     A  4HA   ,4HD   ,4HO   /
      DATA IEQU / 5,6,1,9,4,3,2,7,8,10,11,15,12,13,14 /
C
      DO 1000 I=1,15
         IE=IEQU(I)
         A(I)=B(IE)
 1000 CONTINUE
C
C**   IOP EQUAL ZERO READ NAMELIST WGPLNF
C**   IOP EQUAL ONE WRITE NAMELIST WGPLNF
C
      IF(IOP .EQ. 0)
     1   CALL NAMER(KAND,9,NLNAME,6,PLNF,81,LEN5,15,LDM5,A,15,LOC5,EOF)
      IF(IOP .EQ. 1)
     1   CALL NAMEW(KAND,6,NLNAME,6,PLNF,81,LEN5,15,LDM5,A,15,LOC5)
C
      DO 1010 I=1,15
         IE=IEQU(I)
C
C  THIS IF CHECK IS HERE BECAUSE VARIABLE TEN IS A "DUMB" VARIABLE. ELEM
C  TEN OF COMMON WINGI IS READ IN NAMELIST WGSCHR.
         IF(IE .EQ. 10)GO TO 1010
         B(IE)=A(I)
 1010 CONTINUE
C
      WNG(106) = B(7)
      WNG(112) = B(8)
      WNG(138) = 0.0
      RETURN
      END
