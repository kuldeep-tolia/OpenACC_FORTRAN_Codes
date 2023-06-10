! OpenACC tutorial for kernel construct
! Check how it auto-parallelizes
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER :: a(10), i
  
  !$acc kernels
  PRINT*, "One"

  DO i = 1, 10
     a(i) = i
  END DO
  
  PRINT*, "Two: a(10) = ", a(10)

  DO i = 1, 10
     a(i) = a(i) * 2
  END DO

  PRINT*, "Three: a(10) = ", a(10)
  !$acc end kernels

END PROGRAM MAIN

