      SUBROUTINE WTGEOM(A,AIN)
C
C     ----SUBROUTINE COMPUTES WING OR TAIL GEOMETRY
C
      DIMENSION A(195) , FN(4) , AIN(77)
      LOGICAL LOGSWT
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      DATA FN / 0.,.25,.50,1./
      EQ1(A1,A2) = 4.*(A1**2)/A2
      EQ2(A1,A2) = 4.*(1.-A1)/(A2*(1.+A1))
      EQ4(A1,A2) = 2.*A1*(1.+A2 + A2**2) /(3.*(1.+A2))
      LOGSWT = .FALSE.
      IF(AIN(2).LT.10.*UNUSED)LOGSWT=.TRUE.
      A(21)=AIN(4)-AIN(2)
      A(23)=AIN(3)-AIN(2)
      A(19)=A(23)/A(21)
      IF(LOGSWT)AIN(5)=AIN(1)
      A(25)=AIN(5)/AIN(6)
      A(10)=AIN(6)*(A(25)+(1.-A(25))*A(19))
      A(26)=AIN(5)/A(10)
      A(28) = 1.
      IF(AIN(5).NE.0.) A(28) = AIN(1)/AIN(5)
      A(27)=A(26)*A(28)
      A(118)=AIN(1)/AIN(6)
      A(1)=(A(10)+AIN(5))*A(23)
      A(2)=(AIN(5)+AIN(1))*AIN(2)
      A(3)=A(1)+A(2)
      A(119)= (AIN(6)+AIN(5))*A(21)
      A(4)=A(119)+A(2)
      A(5)=EQ1(A(23),A( 1))
      A(6)=1.
      IF(.NOT.LOGSWT)A(6)=EQ1(AIN(2),A(2))
      A(7)=EQ1(AIN(3),A(3))
      A(120)=EQ1(AIN(4),A(4))
      CI=EQ2(A(26),A(5))
      CO=0.
      IF(.NOT.LOGSWT) CO=EQ2(A(28),A(6))
C
C     ----SWEEPANGLE SECTION
C
      XOVCO = AIN(66)
      DO 1080 I=1,5
         IF(I.LT.5) GO TO 1000
C
C     ----HERE TO COMPUTE SWEEP ANGLES AT MAXIMUM THICKNESS
C
         INDXE = 175
         INDXI = 187
         INDXO = 181
         TMPI = AIN(9)-A(174)
         TMPO = AIN(9)-XOVCO
         GO TO 1010
 1000    INDXE = 34 + (I-1)*6
         INDXI = INDXE + 24
         INDXO = INDXI + 24
         TMPO = AIN(9)-FN(I)
         TMPI = TMPO
C
C     ----TEST FOR STATION = INPUT STATION
C
 1010    IF(ABS(TMPI).LT.1.E-6) GO TO 1060
         A(INDXI+4)=CI*TMPI+A(110)
         CALL ANGLES(5,A(INDXI))
         IF(LOGSWT) GO TO 1030
         A(INDXO+4)=CO*TMPO+A(116)
         CALL ANGLES(5,A(INDXO))
 1020    A(INDXE+3) = (A(1)*A(INDXI+3) + A(2)*A(INDXO+3))/ A(3)
         CALL ANGLES(4,A(INDXE))
         GO TO 1080
C
C     ----HERE IF NO OUTER PANEL
C
 1030    CALL ZERANG(A(INDXO))
 1040    DO 1050 J=1,6
 1050       A(INDXE+J-1) = A(INDXI+J-1)
         GO TO 1080
C
C     ----HAVE THIS STATION ALREADY
C
 1060    DO 1070 J=1,6
            A(INDXI+J-1) = A(J+105)
 1070    A(INDXO+J-1) = A(J+111)
         IF(LOGSWT)GO TO 1040
         GO TO 1020
 1080 CONTINUE
      A(15) = 2.*A(10)*(1.+A(26)+A(26)**2)/(3.*(1.+A(26)))
      A(121)=2.*AIN(6)*(1.+A(25)+A(25)**2)/(3.*(1.+A(25)))
      A(17) = 2.*AIN(5) *(1.+A(28)+A(28)**2)/(3.*(1.+A(28)))
      A(16) = (A(1)*A(15) + A(2)*A(17))/A(3)
      A(122)=(A(119)*A(121)+A(2)*A(17))/A(4)
      A(32) = A(23)*(1.+2.*A(26))/(3.*(1.+A(26)))
      A(33) = AIN(2)*(1.+2.*A(28))/(3.*(1.+A(28)))+ A(23)
      A(31) = (A(1)*A(32) + A(2)*A(33))/A(3)
      A(30) = A(16)/2. + (A(1)*A(32)*A(62)+ A(2)*(A(23)*A(62)+(A(33)-
     1        A(23))*A(86))) / A(3)
      A(18) = (A(23)*A(62) + AIN(2)*A(86)) /A(10)
      A1=A(23)*A(62)
      A2=AIN(2)*A(86)+A1
      A3=A(23)*A(80)
      A4=AIN(2)*A(104)+A3
      A(29)=A(10)+AMAX1(0.,A4,A3)-AMIN1(0.,A2,A1)
      A(162)=A(23)/AIN(3)
      A(163)=4.*(A(21)**2)/A(119)
      A(164)=A(23)/2.
      A(165)=A(164)+AIN(2)
      IF(LOGSWT) GO TO 1090
C
C     ----HERE ONLY IF OUTBOARD SECTION
C
      A(166)=AIN(1)+A(165)*((AIN(5)-AIN(1))/AIN(2))
      A(167)=(A(166)+AIN(1))*A(165)
      A(168)= 4.*( A(165)**2)/A(167)
      A(169)=AIN(1)/A(166)
 1090 CONTINUE
      A(130)=A(21)*(1.+2.*A(25))/(3.*(1.+A(25)))
      A(133)=0.
      IF(.NOT.LOGSWT) A(133)=AIN(2)*(1.+2.*A(28))/(3.*(1.+A(28)))+A(21)
      A(136)=(A(119)*A(130)+A(2)*A(133))/A(4)
      A(195)=(A(119)*A(130)*A(62)+A(2)*(A(21)*A(62)+(A(133)-A(21))*
     1        A(86)))/A(4)
      A(161)=A(122)/4.+A(195)
      RETURN
      END
