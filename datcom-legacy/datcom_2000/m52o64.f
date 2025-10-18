      SUBROUTINE M52O64                                                 M52O64
C                                                                       M52O64
C***  EXEC FOR OVERLAY 52, SUBSONIC LATERAL CONTROL/FLAP AERO           M52O64
C                                                                       M52O64
      COMMON /OVERLY/  IJKDUM(8),NOVLY                                  M52O64
      NOVLY=52                                                          M52O64
      CALL LATFLP                                                       M52O64
      RETURN                                                            M52O64
      END                                                               M52O64
