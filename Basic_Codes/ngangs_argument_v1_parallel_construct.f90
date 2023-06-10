! OpenACC tutorial for parallel construct
! num_gangs taken as input argument from user

PROGRAM MAIN
  IMPLICIT NONE

  INTEGER :: ngangs = 1
  CHARACTER(LEN = 100) :: numchar, nameq

  CALL GET_COMMAND_ARGUMENT(0, nameq)
  CALL GET_COMMAND_ARGUMENT(1, numchar)
  READ(numchar, *) ngangs

  !$acc parallel num_gangs(ngangs)
  PRINT*, "Hello World"
  PRINT*, "Bye World"
  !$acc end parallel

  PRINT*, "Host"

  !$acc parallel num_gangs(ngangs/2)
  CALL PRINT_SUB()
  !$acc end parallel

CONTAINS

  SUBROUTINE PRINT_SUB()
    IMPLICIT NONE

    PRINT*, "This is print from function call"
    
  END SUBROUTINE PRINT_SUB
  
END PROGRAM MAIN
