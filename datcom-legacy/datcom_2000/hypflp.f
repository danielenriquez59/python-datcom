      SUBROUTINE HYPFLP
C
C***  COMPUTES HYPERSONIC FLAP CONTROL AERO
C
      REAL LFIOD0,LFI
      REAL MACH,MALP(20),M2
      REAL MADF
      LOGICAL LAMNR
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA
      COMMON /OPTION/ SR,CBAR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLAPIN/ F(16)
      COMMON /FLGTCD/ FLC(93)
      COMMON /IBODY/  PBODY, BODY(200)
      COMMON /IWING/  PWING, WING(200)
      COMMON /IHT/    PHT, HT(200)
      COMMON /IVT/    PVT, VT(200)
      COMMON /BDATA/  HYP(80)
      DIMENSION DCA(200),DCN(200),DCMCA(200),DCMCN(200)
      DIMENSION ALPHA(20),PAOPI(20),TAOTI(20),RAORI(20),A1(8),A2(8),
     1          A3(8)
      EQUIVALENCE (LAMNR,F(15))
      EQUIVALENCE (XHL,F(2)),(TWOTI,F(3))
      EQUIVALENCE (ALPHA(1),FLC(23))
      EQUIVALENCE (DCN(1),BODY(1)),(DCA(1),WING(1)),(DCMCN(1),HT(1))
     1,(DCMCA(1),VT(1)),(PAOPI(1),HYP(1)),(TAOTI(1),HYP(21)),
     2 (MALP(1),HYP(41)),(RAORI(1),HYP(61))
      EQUIVALENCE (ALT,F(1)),(CF,F(4))
C
      CPIP=0.0
      MACH=FLC(I+2)
      RL=FLC(I+42)
      NDELTA=F(16)+.5
      DO 1200 J=1,NALPHA
         IF(ALPHA(J).LT.0.OR.ALPHA(J).GT.20.0)GO TO 1200
         CALL HYPROP(ALT,MACH,ALPHA(J),PAOPI(J),TAOTI(J),MALP(J),
     1               RAORI(J))
         CPIA=(PAOPI(J)-1.)/(.7*MACH**2)
         RLHL=XHL*RAORI(J)*RL
         M2=MALP(J)**2
         DO 1190 N=1,NDELTA
            DELTA=F(N+4)
            SINDF=SIN(DELTA/RAD)
            COSDF=COS(DELTA/RAD)
            IF(.NOT.LAMNR)GO TO 1020
C
C    ----HERE FOR LAMINAR BOUNDARY LAYER
C    ----DETERMINE IF FLOW SEPARATION EXIST
C
            CPINC=2.03*(MALP(J)**2-1.)**(-0.306)/RLHL**0.25
C
C              FIGURE 6.3.1-58 SINGLE SHOCK (EQUATION)
C
            CALL FIG68(MALP(J),DELTA,THETA,IER)
            IF(IER .EQ. 2) THETA=90.0
            P2PASS=(7.*M2*SIN(THETA/RAD)**2-1.)/6.
            IF(DELTA.EQ.0.0)P2PASS=1.0
            IF(MALP(J).GT.6.0)GO TO 1000
C
C              FIGURE 6.3.1-58 ISENTROPIC
C
            ARG=0.2*M2*SINDF/SQRT(M2-1.)
            P2PAIS=(1.+ARG)**7
            P2PA=0.5*(P2PASS+P2PAIS)
            GO TO 1010
 1000       P2PA=P2PASS
 1010       CPA2=(P2PA-1.)/(.7*MALP(J)**2)
            IF(CPA2.GE.CPINC)GO TO 1030
            CPAP=CPIA
            LFI=0.0
            D1=0.0
            D2=0.0
            D3=0.0
            GO TO 1100
 1020       CONTINUE
C
C    ----HERE FOR TURBULENT BOUNDARY LAYER
C    ----DETERMINE IF FLOW SEPARATION EXIST
C              FIGURE 6.3.1-57 (EQUATION)
C
            CPINC=2.2/RLHL**0.1
C
C        FIGURE 6.3.1-58 SINGLE SHOCK (EQUATION)
C
            CALL FIG68(MALP(J),DELTA,THETA,IER)
            P2PASS=(7.*M2*SIN(THETA/RAD)**2-1.)/6.
            IF(IER .EQ. 2) THETA=90.0
            IF(DELTA.EQ.0.0)P2PASS=1.0
            P2PA=P2PASS
            CPA2=(P2PA-1.)/(.7*MALP(J)**2)
            IF(CPA2.GE.CPINC)GO TO 1030
            CPAP=CPIA
            LFI=0.0
            D1=0.0
            D2=0.0
            D3=0.0
            GO TO 1100
 1030       CONTINUE
C
C    ----HERE FOR SEPARATED FLOW
C    ----DETERMINE SEPARATION LOCATION X0
C
            X8=XHL/8.0
            X0=XHL/20.
            ARG=RAORI(J)*RL
            DO 1060 K=1,8
               X0=X0+X8
               RLX0=ARG*X0
               TWOTA=TWOTI/TAOTI(J)
               IF(LAMNR)GO TO 1040
C
C    ----HERE FOR TURBULENT FLOW
C              FIGURE 6.3.1-60 (EQUATION)
C
               CPAP=1.91*(MALP(J)**2-1.)**(-0.309)/RLX0**0.1
C
C**       FIGURE 6.3.1-62 (EQUATION)
C
               ARG4=0.7*CPAP*MALP(J)**2+1.
               D1OD0=1.1E6*(MALP(J)**(-1.67)*(ARG4-1.))**8.55
C
C--       FIGURE 6.3.1-64 (EQUATION)
C
               IF(K.EQ.1)RLSTRL= (0.28+0.5*TWOTA+0.22*(1.+0.2*MALP(J)
     1                   **2*0.72**0.3333))
               RLSTR=RLSTRL*ARG
C
C              FIGURE 6.3.1-66 (EQUATION)
C
               D0=0.154/RLSTR**0.142857*X0**0.85714
               GO TO 1050
 1040          CONTINUE
C
C    ----HERE FOR LAMINAR FLOW
C           FIGURE 6.3.1-59 (EQUATION)
C
               ARG4=1.56*(MALP(J)**2-1.)**(-0.262)
               CPAP=ARG4/RLX0**0.25
C
C           FIGURE 6.3.1-61 (EQUATION)
C
               ARG4=0.7*CPAP*MALP(J)**2+1.
               D1OD0=5.69E5*MALP(J)**(-4.1)*(ARG4-1.)**3.5
C
C           FIGURE 6.3.1-63 (EQUATION)
C
               IF(K.EQ.1)RLSTRL=(0.28+0.5*TWOTA+0.22*(1.+0.2*MALP(J)
     1                         **2*0.72**0.3333))**(-1.67)
               RLSTR=RLSTRL*ARG
C
C           FIGURE 6.3.1-65 (EQUATION)
C
               D0=5.2*SQRT(X0)/SQRT(RLSTR)
 1050          D1=D1OD0*D0
               A1(K)=X0+D1
               A2(K)=XHL
               A3(K)=X0
 1060       CONTINUE
            CALL SIMUL2(A3,A2,A1,8,X0,X0PD1)
            IF(X0.EQ.(-1000.))X0=0.0
            IF(X0.EQ.0.0)GO TO 1110
            D1=X0PD1-X0
C
C    ----CALCULATE PLATEAU-PRESSURE LEVEL AND INTERACTION LENGTHS
C
            RLX0=ARG*X0
            ARG1=(MALP(J)/MACH)**2*PAOPI(J)
            RLSTR=RLSTRL*ARG
            IF(LAMNR)GO TO 1070
C
C    ----HERE FOR TURBULENT
C              FIGURE 6.3.1-66 (EQUATION)
C
            D0=0.154/RLSTR**0.142857*X0**0.85714
C
C              FIGURE 6.3.1-60 (EQUATION)
C
            CPAP=1.91*(MALP(J)**2-1.)**(-0.309)/RLX0**0.1
            CPIP=CPAP*ARG1+CPIA
C
C              FIGURE 6.3.1-68 (EQUATION)
C
            ARG4=0.7*CPAP*MALP(J)**2+1.
            LFIOD0=1.84E4*((ARG4-1.)/MALP(J)**1.325)**8.4
            IF(MALP(J).GE.3.9)LFIOD0=D1/D0
C
C              FIGURE 6.3.1-62 (EQUATION)
C
            D1OD0=1.1E6*(MALP(J)**(-1.67)*(ARG4-1.))**8.55
            GO TO 1080
 1070       CONTINUE
C
C    ----HERE FOR LAMINAR
C              FIGURE 6.3.1-65 (EQUATION)
C
            D0=5.2*SQRT(X0)/SQRT(RLSTR)
C
C              FIGURE 6.3.1-59 (EQUATION)
C
            ARG4=1.56*(M2-1.)**(-0.262)
            CPAP=ARG4/RLX0**0.25
            CPIP=CPAP*ARG1+CPIA
C
C              FIGURE 6.3.1-67 (EQUATION)
C
            PPPA=0.7*CPAP*M2+1.
            LFIOD0=2.47E5*MALP(J)**(-4.2)*(PPPA-1.)**3.45
C
C              FIGURE 6.3.1-61 (EQUATION)
C
            D1OD0=5.69E5*MALP(J)**(-4.1)*(PPPA-1.)**3.5
 1080       CONTINUE
            LFI=LFIOD0*D0
            D1=D1OD0*D0
            IF(LAMNR)GO TO 1090
C
C    ----SET D3=0 FOR TURBULENT B.L.
C
            D3=0.0
            GO TO 1100
 1090       PPPA=0.7*CPAP*MALP(J)**2+1.
            ARG=6.0*PPPA+1.0
            ARG=SQRT(ARG/(7.*MALP(J)**2))
            THETA=ATAN(ARG/SQRT(1.-ARG**2))
            ARG=M2*SIN(2.*THETA)-2./TAN(THETA)
            ARG1=10.+M2*(7.+5.*COS(2.*THETA))
            TANPHE=5.*ARG/ARG1
            PHE=ATAN(TANPHE)*RAD
 1100       CPI2=CPA2*(MALP(J)/MACH)**2*PAOPI(J)+CPIA
 1110       CONTINUE
            D3=0.0
            IF(CPA2.LT.CPINC)GO TO 1180
            CFOD1=CF/D1
            MADF=MALP(J)*DELTA/RAD
            IF(.NOT.LAMNR)GO TO 1130
C
C              FIGURE 6.3.1-69 (EQUATIONS)
C
            IF(MADF.GE.5.)GO TO 1120
            D2=D1*(0.545-0.0403*MADF)*SQRT(CFOD1)
            IF(CFOD1.GE.1.)D2=D1*(0.545-0.04*MADF)
            IF(CFOD1.LE.0.25)D2=D1*(.273-0.02*MADF)
            GO TO 1150
 1120       D2=.344*D1*SQRT(CFOD1)
            IF(CFOD1.GE.1.)D2=D1*0.344
            IF(CFOD1.LE.0.25)D2=D1*0.172
            GO TO 1150
C
C              FIGURE 6.3.1-70 (EQUATIONS)
C
 1130       IF(MADF.GE.2.4)GO TO 1140
            D2=(1.16-0.33*MADF)*SQRT(CFOD1)*D1
            IF(CFOD1.GE.1.)D2=D1*(1.16-0.33*MADF)
            IF(CFOD1.LE.0.25)D2=D1*(0.58-0.165*MADF)
            GO TO 1150
 1140       D2=D1*0.37*SQRT(CFOD1)
            IF(CFOD1 .GE. 1.0)D2=D1*0.37
            IF(CFOD1 .LE. 0.25)D2=D1*0.37*0.5
 1150       IF(X0 .NE. 0.0)GO TO 1160
C
C    ----HERE FOR SEPARATION AT L.E. OF FLAT PLATE
C
            D1=XHL
            LFI=0.0
            PHE=ATAN(D2*SINDF)/(D1+D2*COSDF)*RAD
            CALL FIG68(MALP(J),PHE,THETA,IER)
            IF(IER .EQ. 2) THETA=90.0
            PPPA=(7.*M2*SIN(THETA/RAD)**2-1.)/6.
            IF(PHE.EQ.0.0)PPPA=1.0
            CPAP=(PPPA-1.)/(0.7*M2)
            CPIP=CPAP*(MALP(J)/MACH)**2*PAOPI(J)+CPIA
 1160       IF(LAMNR)GO TO 1170
C
C    ----SET D3=0 FOR TURB. B.L.
C
            D3=0.0
            GO TO 1180
 1170       ARG1=TAN(DELTA/RAD)-TAN(PHE/RAD)
            D3OD1=(TAN(PHE/RAD)/ARG1)/COS(DELTA/RAD)
            D3=D1*D3OD1
 1180       CONTINUE
C
C    ----NORMAL-FORCE INCREMENT
C
            ARG1=(CPIP-CPIA)*(D1-LFI/2.+CF*COSDF)
            ARG2=(CPI2-CPIP)*COSDF*(CF-0.5*(D2+D3))
            K=(J-1)*10+N
            DCN(K)=(ARG1+ARG2)/SR
C
C    ----AXIAL-FORCE INCREMENT
C
            ARG1=SIN(DELTA/RAD)/SR
            ARG2=CF*(CPIP-CPIA)
            ARG3=(CPI2-CPIP)*(CF-0.5*(D2+D3))
            DCA(K)=ARG1*(ARG2+ARG3)
C
C    ----INCREMENT OF PITCHING MOMENT ABOUT H.L. DUE TO NORMAL FORCE
C
            ARG1=(CPIP-CPIA)*(LFI**2/6.+.5*(-LFI*D1+D1**2-
     1           (CF*COSDF)**2))
            ARG2=(CPI2-CPIP)*COSDF**2*(CF**2/2.-( D2**2+D2*D3+D3**2)/6.)
            DCMCN(K)=(ARG1-ARG2)/(SR*CBAR)
C
C    ----INCREMENT OF PITCHING MOMENT ABOUT H.L. DUE TO AXIAL FORCE
C
            ARG1=-SIN(DELTA/RAD)**2/(SR*CBAR)
            ARG2=(CPIP-CPIA)*CF**2/2.0
            ARG3=(CPI2-CPIP)*(CF**2/2.0-(D2**2+D2*D3+D3**2)/6.)
            DCMCA(K)=ARG1*(ARG2+ARG3)
 1190    CONTINUE
 1200 CONTINUE
      RETURN
      END
