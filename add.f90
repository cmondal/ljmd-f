Program add
   Implicit none
   Integer :: i
   Integer, Parameter :: n=10
   Real :: sum1

   sum1 = 0.0
   Do i = 1, n
    sum1 = sum1 + real(i)
   Enddo

   Print*, "sum of numbers from 1 to",n, "is=", sum1

End Program add
