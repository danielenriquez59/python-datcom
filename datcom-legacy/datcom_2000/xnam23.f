      SUBROUTINE XNAM23(IOP)
C
C***  SUBROUTINE TO READ OR WRITE NAMELIST EXPR
C
      COMMON /IBODY/  AHA(141)
      COMMON /IWING/  AHB(141)
      COMMON /IHT/    AHC(141)
      COMMON /IBW/    AHD(141)
      COMMON /IVT/    AHE(2)
      COMMON /IDWASH/ AHF(61)
      COMMON /EXPER/  AHG(116)
      COMMON /SBETA/  STB(135),TRA(108),TRAH(108)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND
C
      INTEGER EXPR
      LOGICAL EOF
C
      DIMENSION LENH(32),LDMH(32),EXPR(131),LOCH(32),AH(469)
      DIMENSION NLNAME(4)
C
C***  NAMELIST EXPR
C
      DATA NLNAME / 4HE   ,4HX   ,4HP   ,4HR   /
      DATA LENH / 3*3,4,4,3*3,4,4,3*3,5*4,5,5,6,6,5,3,5*5,4,5,4 /
      DATA LDMH / 23*20,9*1 /
      DATA LOCH / 1,21,41,61,81,101,121,141,161,181,201,221,241,
     1            261,281,301,321,341,361,381,401,421,441,461,
     2            462,463,464,465,466,467,468,469 /
      DATA EXPR /     4HC   ,4HD   ,4HB   ,4HC   ,4HL   ,4HB   ,
     1  4HC   ,4HM   ,4HB   ,4HC   ,4HL   ,4HA   ,4HB   ,4HC   ,
     2  4HM   ,4HA   ,4HB   ,4HC   ,4HD   ,4HW   ,4HC   ,4HL   ,
     3  4HW   ,4HC   ,4HM   ,4HW   ,4HC   ,4HL   ,4HA   ,4HW   ,
     4  4HC   ,4HM   ,4HA   ,4HW   ,4HC   ,4HD   ,4HH   ,4HC   ,
     5  4HL   ,4HH   ,4HC   ,4HM   ,4HH   ,4HC   ,4HL   ,4HA   ,
     6  4HH   ,4HC   ,4HM   ,4HA   ,4HH   ,4HC   ,4HD   ,4HW   ,
     7  4HB   ,4HC   ,4HL   ,4HW   ,4HB   ,4HC   ,4HM   ,4HW   ,
     8  4HB   ,4HC   ,4HL   ,4HA   ,4HW   ,4HB   ,4HC   ,4HM   ,
     9  4HA   ,4HW   ,4HB   ,4HQ   ,4HO   ,4HQ   ,4HI   ,4HN   ,
     A  4HF   ,4HE   ,4HP   ,4HS   ,4HL   ,4HO   ,4HN   ,4HD   ,
     B  4HE   ,4HO   ,4HD   ,4HA   ,4HC   ,4HD   ,4HV   ,4HA   ,
     C  4HL   ,4HP   ,4HO   ,4HW   ,4HA   ,4HL   ,4HP   ,4HL   ,
     D  4HW   ,4HA   ,4HL   ,4HP   ,4HO   ,4HH   ,4HA   ,4HL   ,
     E  4HP   ,4HL   ,4HH   ,4HA   ,4HC   ,4HL   ,4HM   ,4HW   ,
     F  4HC   ,4HL   ,4HM   ,4HW   ,4HA   ,4HC   ,4HL   ,4HM   ,
     G  4HH   ,4HC   ,4HL   ,4HM   ,4HH   /
C
C**   LOOP TO FILL WORKING ARRAY
C     J PARAMETER IS TO PUT CLA AND CMA FROM COMMON INTO CORRECT POSITIO
C     IN THE WORKING ARRAY (SEE INSERT TO USERS MANUAL MARK UP COPY FOR
C     EXPLANATION OF IOM COMMON BLOCKS WHEN USING EXPERIMENTAL SUBSTIUTI
C
      J=1
      DO 1000 I=1,100
         J=J+1
         AH(I)=AHA(J)
         AH(I+100)=AHB(J)
         AH(I+200)=AHC(J)
         AH(I+300)=AHD(J)
         IF(I .EQ. 60)J=101
 1000 CONTINUE
      DO 1010 I=1,60
         AH(I+400)=AHF(I+1)
 1010 CONTINUE
      AH(461)=AHE(2)
      DO 1020 I=1,4
         AH(I+461)=AHG(I+112)
 1020 CONTINUE
      AH(466)=TRA(37)
      AH(467)=TRA(38)
      AH(468)=TRAH(37)
      AH(469)=TRAH(38)
C
C***  IF IOP EQUAL 0 READ NAMELIST EXPR, IF IOP EQUAL 1 WRITE EXPR
C
      IF(IOP .EQ. 0)
     1 CALL NAMER(KAND,10,NLNAME,4,EXPR,131,LENH,32,LDMH,AH,469,
     2            LOCH,EOF)
      IF(IOP .EQ. 1)
     1 CALL NAMEW(KAND,6,NLNAME,4,EXPR,131,LENH,32,LDMH,AH,469,LOCH)
C
C***  REPLACE DATA READ INTO APPROPRIATE COMMON BLOCKS
C
      J=1
      DO 1030 I=1,100
         J=J+1
         AHA(J)=AH(I)
         AHB(J)=AH(I+100)
         AHC(J)=AH(I+200)
         AHD(J)=AH(I+300)
         IF(I .EQ. 60)J=101
 1030 CONTINUE
      DO 1040 I=1,60
         AHF(I+1)=AH(I+400)
 1040 CONTINUE
      AHE(2)=AH(461)
      DO 1050 I=1,4
         AHG(I+112)=AH(I+461)
 1050 CONTINUE
      TRA(37)=AH(466)
      TRA(38)=AH(467)
      TRAH(37)=AH(468)
      TRAH(38)=AH(469)
C
      RETURN
      END
