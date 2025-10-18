      SUBROUTINE M57O71                                                 M57O71
C                                                                       M57O71
C  PROGRAM TO READ EXTRAPOLATION MESSAGES FROM UNIT 12                  M57O71
C  AND PRINT IN A READABLE FORMAT                                       M57O71
C                                                                       M57O71
      INTEGER ROUT,ROUTL                                                M57O71
      LOGICAL PAIR,PAIRX                                                M57O71
C                                                                       M57O71
C   ROUTINE USES THE IOM AND WORKING DATA ARRAYS TO                     M57O71
C   CONSERVE CORE.  ROUTINE EXECUTES AFTER COMPLETION OF ALL            M57O71
C   CASE COMPUTATIONS.                                                  M57O71
C                                                                       M57O71
      COMMON /IBODY/   PB, AFIG01(121), AFIG02(121)                     M57O71
      COMMON /IWING/   PW, AFIG03(121), AFIG04(121)                     M57O71
      COMMON /IHT/    PHT, AFIG05(121), AFIG06(121)                     M57O71
      COMMON /IVT/    PVT, AFIG07(121), AFIG08(121)                     M57O71
      COMMON /IVF/    PVF, AFIG09(121), AFIG10(121)                     M57O71
      COMMON /IBW/    PBW, AFIG11(121), AFIG12(121)                     M57O71
      COMMON /IBH/    PBH, AFIG13(121), AFIG14(121)                     M57O71
      COMMON /IBV/    PBV, JFIG01(121), JFIG02(121)                     M57O71
      COMMON /IBWH/  PBWH, JFIG03(121), JFIG04(121)                     M57O71
      COMMON /IBWV/  PBWV, JFIG05(121), JFIG06(121)                     M57O71
      COMMON /IBWHV/PBWHV, JFIG07(121), JFIG08(121)                     M57O71
      COMMON /WINGD/       JFIG09(121)                                  M57O71
      COMMON /SBETA/       JFIG10(121), JFIG11(121), JFIG12(121)        M57O71
      COMMON /BDATA/       JFIG13(121), JFIG14(121), JFIG15(121),       M57O71
     1                     JFIG16(121), JFIG17(121)                     M57O71
      COMMON /WHWB/        JFIG18(121)                                  M57O71
      COMMON /WBHCAL/      JFIG19(121)                                  M57O71
      COMMON /HTDATA/      JFIG20(121)                                  M57O71
      COMMON /VTDATA/ LFIGN(121),LFIGO(121),IOVLY,IOVL,NSTQ,NSTP,FINALR,M57O71
     1                NFIG,IFIGN                                        M57O71
      COMMON /WHAERO/ LFIGS(121), IFIGST(20), IEXCD(4), MSSCL(2),       M57O71
     1                ROUT(2), XLL(4), XUL(4), IEXTRL(4,2),             M57O71
     2                IEXTRU(4,2), XVAL(4), LDUM(121), PAIRX            M57O71
      COMMON /ERROR/ DERROR(82),LEXTRP                                  M57O71
      LOGICAL DERROR,LEXTRP                                             M57O71
      DIMENSION IFIG(20,121)                                            M57O71
C                                                                       M57O71
C*** IF PRINT EXTRAP CARD NOT INPUT REURN AND DO NOT DUMP EXTRAPOLATION M57O71
C       MESSAGES                                                        M57O71
      IF(.NOT. LEXTRP)RETURN                                            M57O71
C                                                                       M57O71
C  INITIALIZE ALL ARRAYS AND VARIABLES TO BE USED                       M57O71
C                                                                       M57O71
      CALL CLEARA(IFIG,ROUTL,MSSCLL)                                    M57O71
C                                                                       M57O71
C  WRITE OVERLAY NUMBER OF 999 TO UNIT 12 TO SIGNIFY END OF DATA        M57O71
C  THIS OVERLAY NUMBER HAS BEEN PRE-SET IN SUBROUTINE CLEARA            M57O71
C                                                                       M57O71
      WRITE(12,1090)IOVLY,IOVLY,IOVLY                                   M57O71
C                                                                       M57O71
C   REWIND UNIT 12 AND START READING EXTRAPOLATION MESSAGES             M57O71
C                                                                       M57O71
      REWIND 12                                                         M57O71
      IPRT=60                                                           M57O71
 1000 CALL READXM(IFIG)                                                 M57O71
      IF(IOVLY.NE.IOVL .OR. NFIG.GE.120 .OR. IOVL.EQ.999) GO TO 1020    M57O71
 1010 CONTINUE                                                          M57O71
C                                                                       M57O71
C  TRANSLATE FIGURE NUMBER INTO 9 DIGIT INTEGER WORD                    M57O71
C                                                                       M57O71
      CALL DECFIG(IFIG(1,NFIG),20,IFIGN)                                M57O71
C                                                                       M57O71
C  STORE FIGURE DATA FOR LATER USE                                      M57O71
C                                                                       M57O71
      CALL STORXM                                                       M57O71
      GO TO 1000                                                        M57O71
C                                                                       M57O71
C  SORT FIGURE NUMBERS FOR THIS OVERLAY IN ASCENDING ORDER              M57O71
C                                                                       M57O71
 1020 CONTINUE                                                          M57O71
      IOVM=IOVL                                                         M57O71
      NFIG=NFIG-1                                                       M57O71
      IF(NFIG.LE.0 .AND. IOVLY.EQ.999)GO TO 1060                        M57O71
      CALL SORTER(LFIGN,NFIG,1,0,1,LFIGS,LFIGO,LDUM,IFLAG)              M57O71
C                                                                       M57O71
C  PULL FIGURES IN ASCENDING ORDER                                      M57O71
C                                                                       M57O71
      PAIRX=.TRUE.                                                      M57O71
      DO 1040 I=1,NFIG                                                  M57O71
         IF(IPRT.GE.54)WRITE(6,1070)                                    M57O71
         IF(IPRT.GE.54)WRITE(6,1080)                                    M57O71
         IF(IPRT.GE.54)PAIR=.TRUE.                                      M57O71
         IF(IPRT.GE.54)IPRT=6                                           M57O71
C                                                                       M57O71
C  FIND WHERE THIS FIGURE IS STORED IN THE ARRAYS                       M57O71
C                                                                       M57O71
         JJ=LFIGO(I)                                                    M57O71
         JK=AFIG01(JJ)+0.5                                              M57O71
         CALL WRITXM(IFIG,IPRT,JK,JJ,IOVM,ROUTL,JJ,MSSCLL)              M57O71
         IF(IPRT.GE.54)WRITE(6,1100)                                    M57O71
         IF(IPRT.GE.54)PAIR=.FALSE.                                     M57O71
C                                                                       M57O71
C  SET TEST ARRAYS TO PRESENT FIGURE NUMBER AND ROUTINES                M57O71
C                                                                       M57O71
         DO 1030 K=1,20                                                 M57O71
            IFIGST(K)=IFIG(K,JJ)                                        M57O71
 1030    CONTINUE                                                       M57O71
         MSSCLL=JFIG01(JJ)                                              M57O71
         ROUTL=JFIG03(JJ)                                               M57O71
 1040 CONTINUE                                                          M57O71
      IF(IOVLY.EQ.999 .OR. IOVLY.LE.0)GO TO 1060                        M57O71
      NFIG=NFIG+1                                                       M57O71
      DO 1050 L=1,20                                                    M57O71
         IFIG(L,1)=IFIG(L,NFIG)                                         M57O71
 1050 CONTINUE                                                          M57O71
      NFIG=1                                                            M57O71
C                                                                       M57O71
C  CONTINUE PROCESSING WITH LAST FIGURE READ FROM UNIT 12               M57O71
C                                                                       M57O71
      GO TO 1010                                                        M57O71
 1060 CONTINUE                                                          M57O71
C                                                                       M57O71
C  COMPLETED THE DUMP OF EXTRAPOLATION MESSAGES FOR THIS CASE           M57O71
C  REWIND UNIT 12 TO BEGIN NEXT CASE                                    M57O71
C                                                                       M57O71
      REWIND 12                                                         M57O71
      IF(PAIR)WRITE(6,1100)                                             M57O71
 1070 FORMAT(1H1,50X,29HEXTRAPOLATION MESSAGE SUMMARY/)                 M57O71
 1080 FORMAT(1X,7HOVERLAY,2X,13HFIGURE NUMBER,                          M57O71
     1 26X,35HTYPE OF EXTRAPOLATION (LOWER UPPER)/                      M57O71
     2 11X,11HSUBROUTINES,32X,27HFIGURE LIMITS (LOWER UPPER)/           M57O71
     3 25X,12HFINAL RESULT,20X,21HINDEPENDENT VARIABLES)                M57O71
 1090 FORMAT(3I3)                                                       M57O71
 1100 FORMAT(107H0*NOTE*   MESSAGES HAVE BEEN SORTED WITHIN EACH OVERLAYM57O71
     1 BY FIGURE NUMBER AND DO NOT NECESSARILY CORRESPOND               M57O71
     2 /10X,34HTO THE PROGRAM EXECUTION SEQUENCE.                       M57O71
     3 /9X,53H ** DENOTES THE INDEPENDENT VARIABLE EXCEEDED IN CALL)    M57O71
      RETURN                                                            M57O71
      END                                                               M57O71
