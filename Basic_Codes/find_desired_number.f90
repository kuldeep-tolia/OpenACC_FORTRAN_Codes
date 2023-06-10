! OpenACC tutorial: Check and find the desired number in an array
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER, PARAMETER    :: N = 1000 ! set the array size
  INTEGER               :: myNumber, j
  INTEGER, DIMENSION(N) :: numbers

  CALL INITIALIZE_NUM()

  myNumber = numbers(N/2)

  !$acc parallel loop
  DO j = 1, N

     IF(numbers(j) .EQ. myNumber) PRINT*, "Number found at ", j

  END DO
  !$acc end parallel loop

CONTAINS

  SUBROUTINE INITIALIZE_NUM()
    IMPLICIT NONE

    INTEGER :: i

    DO i = 1, N

       numbers(i) = (i + 5) / 5

    END DO
    
  END SUBROUTINE INITIALIZE_NUM
  
END PROGRAM MAIN

