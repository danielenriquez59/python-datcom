      SUBROUTINE XNAM1(IOP)
C
C***  READ OR WRITE NAMELIST FLTCON
C
      COMMON /FLGTCD/ A2(160)
      COMMON /FLOLOG/ FLO(19)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER FLTCON
      LOGICAL EOF
C
      DIMENSION LEN2(18),LDM2(18),FLTCON(81),LOC2(18),A2A(161)
      DIMENSION NLNAME(6)
C
      DATA NLNAME / 4HF   ,4HL   ,4HT   ,4HC   ,4HO   ,4HN   /
      DATA LEN2 / 5,4,6,6,5,3*6,2,3,3*4,2,5,4,4,5 /
      DATA LDM2 / 1,20,1,20,20,-1,3*1,4*20,4*1,20 /
      DATA LOC2 / 1,3,2,23,43,161,94,95,96,97,74,117,137,157,158
     1            ,159,160,23 /
      DATA FLTCON /   4HN   ,4HM   ,4HA   ,4HC   ,4HH   ,4HM   ,
     1  4HA   ,4HC   ,4HH   ,4HN   ,4HA   ,4HL   ,4HP   ,4HH   ,
     2  4HA   ,4HA   ,4HL   ,4HS   ,4HC   ,4HH   ,4HD   ,4HR   ,
     3  4HN   ,4HN   ,4HU   ,4HB   ,4HH   ,4HY   ,4HP   ,4HE   ,
     4  4HR   ,4HS   ,4HS   ,4HT   ,4HM   ,4HA   ,4HC   ,4HH   ,
     5  4HT   ,4HS   ,4HM   ,4HA   ,4HC   ,4HH   ,4HT   ,4HR   ,
     6  4HA   ,4HL   ,4HT   ,4HP   ,4HI   ,4HN   ,4HF   ,4HT   ,
     7  4HI   ,4HN   ,4HF   ,4HV   ,4HI   ,4HN   ,4HF   ,4HW   ,
     8  4HT   ,4HG   ,4HA   ,4HM   ,4HM   ,4HA   ,4HN   ,4HA   ,
     9  4HL   ,4HT   ,4HL   ,4HO   ,4HO   ,4HP   ,4HA   ,4HL   ,
     A  4HP   ,4HH   ,4HA   /
C
      DO 1000 I=1,160
         A2A(I)=A2(I)
 1000 CONTINUE
      A2A(161)=FLO(19)
C
C**   IF IOP EQUAL ZERO READ NAMELIST FLTCON
C**   IF IOP EQUAL ONE WRITE NAMELIST FLTCON
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,9,NLNAME,6,FLTCON,81,LEN2,18,LDM2,A2A,161,
     2            LOC2,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,6,FLTCON,81,LEN2,18,LDM2,A2A,161,LOC2)
C
      DO 1010 I=1,160
         A2(I)=A2A(I)
 1010 CONTINUE
      FLO(19)=A2A(161)
C
      RETURN
      END
