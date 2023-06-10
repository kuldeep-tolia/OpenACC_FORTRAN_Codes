! Hello world program with OpenACC

PROGRAM MAIN
  IMPLICIT NONE

  INTEGER :: i
  
  !$acc parallel loop
  DO i = 1, 1000
     PRINT*, "Hello World ", i
  END DO
    
END PROGRAM MAIN
