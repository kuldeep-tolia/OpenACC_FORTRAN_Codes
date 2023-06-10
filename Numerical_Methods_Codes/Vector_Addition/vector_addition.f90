! OpenACC tutorial for vector addition
! Use -Minfo=accel to see the compiler behaviour
! Use export PGI_ACC_TIME=1 to profile the time behaviour of the code
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 1000
  INTEGER, DIMENSION(:) :: A(N), B(N), C(N)
  INTEGER :: i

  !$acc data copyin(A, B) copyout(C)
  CALL INIT()
  CALL ADD()
  !$acc end data

  DO i = 1, N
     WRITE(*, *) C(i)
  END DO
  
CONTAINS

  SUBROUTINE INIT()
    IMPLICIT NONE

    INTEGER :: i

    !$acc parallel loop present(A, B)
    DO i = 1, N
       A(i) = i
       B(i) = i
    END DO
    !$acc end parallel    

  END SUBROUTINE INIT

  SUBROUTINE ADD()
    IMPLICIT NONE

    INTEGER :: i

    !$acc parallel loop present(A, B, C)
    DO i = 1, N
       C(i) = A(i) + B(i)
    END DO
    !$acc end parallel
    
  END SUBROUTINE ADD
  
END PROGRAM MAIN
