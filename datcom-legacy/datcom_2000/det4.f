      SUBROUTINE DET4(A,D)
C
C     ----COMPUTES 4X4 DETERMINANT
C
      DIMENSION A(16),A3(9)
      P=0.
      DO 1020 M=1,4
         ISKP=M+4
         KK=1
         K=1
         DO 1010 I=5,16
            IF(I.NE.ISKP)GO TO 1000
            ISKP=ISKP+4
            GO TO 1010
 1000       A3(K)=A(I)
            K=K+1
 1010    KK=KK+1
         PP=A3(1)*(A3(5)*A3(9)-A3(6)*A3(8))-A3(2)*(A3(4)*A3(9)-A3(6)
     1      *A3(7))+ A3(3)*(A3(4)*A3(8)-A3(5)*A3(7))
         IF(M.EQ.2.OR.M.EQ.4) PP=-PP
 1020 P=P+A(M)*PP
      D=P
      RETURN
      END
