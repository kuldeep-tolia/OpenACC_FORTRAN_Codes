! OpenACC tutorial for parallel construct
! Check async host-device printing
  
PROGRAM MAIN
  IMPLICIT NONE

  !$acc parallel num_gangs(10) async
  PRINT*, "Hello"
  PRINT*, "Bye"
  !$acc end parallel

  PRINT*, "Host"

  !$acc parallel num_gangs(3)
  PRINT*, "One"
  PRINT*, "Two"
  !$acc end parallel

END PROGRAM MAIN

