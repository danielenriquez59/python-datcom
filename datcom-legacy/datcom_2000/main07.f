      SUBROUTINE MAIN07                                                 MAIN07
C                                                                       MAIN07
C***  SUPERSONIC HIGH LIFT AND CONTROL DEVICES EXECUTIVE                MAIN07
C                                                                       MAIN07
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,  MAIN07
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,   MAIN07
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,        MAIN07
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART              MAIN07
C                                                                       MAIN07
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,         MAIN07
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,          MAIN07
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,               MAIN07
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART                     MAIN07
C                                                                       MAIN07
      IF(SYMFP) CALL M41O51                                             MAIN07
      IF(ASYFP) CALL M53O65                                             MAIN07
      RETURN                                                            MAIN07
      END                                                               MAIN07
