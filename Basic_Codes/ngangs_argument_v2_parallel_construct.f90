! OpenACC tutorial for parallel construct
! Check how many times printf is executed depending on #gangs used

PROGRAM MAIN
  IMPLICIT NONE

  !$acc parallel num_gangs(2)
  PRINT*, "How many times am I printed"
  !$acc end parallel
  
END PROGRAM MAIN
  
