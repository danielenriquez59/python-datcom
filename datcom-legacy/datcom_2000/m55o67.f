      SUBROUTINE M55O67                                                 M55O67
C                                                                       M55O67
C***  EXEC FOR OVERLAY 55, CALCULATE INCREMENTAL AFFECTS DUE TO JET FLAPM55O67
C                                                                       M55O67
      COMMON /OVERLY/  IJKDUM(8),NOVLY                                  M55O67
      NOVLY=55                                                          M55O67
      CALL JETFP                                                        M55O67
      RETURN                                                            M55O67
      END                                                               M55O67
