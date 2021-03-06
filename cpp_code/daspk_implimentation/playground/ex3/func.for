      REAL FUNCTION HORNER(A,N,X)
C     This function returns the value of the polynomial 
C     y = a_0 + a_1 x + a_2 x^2 + … + a_n x^n
C     using Horner's method.
      INTEGER I,N
      REAL*4    A(0:N),X
      
      HORNER = A(N)
      DO 10 I = N-1,0,-1
      HORNER = A(I) + HORNER*X
 10   CONTINUE
      END
