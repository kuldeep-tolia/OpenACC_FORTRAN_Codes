! OpenACC stencil computation code using in-situ update
! Code profiling done using PGI_ACC_TIME environment variable

PROGRAM MAIN
  IMPLICIT NONE

  INTEGER, PARAMETER    :: N = 256             ! set number of rows
  INTEGER, PARAMETER    :: M = 256             ! set number of columns
  INTEGER, PARAMETER    :: num_iterations = 10 ! set number of iterations

  INTEGER, DIMENSION(:, :), ALLOCATABLE :: A, A_old

  INTEGER :: iter
  DOUBLE PRECISION :: start_t, end_t

  ALLOCATE (A(1:N, 1:M), A_old(1:N, 1:M))

  CALL CPU_TIME(start_t)
  !$acc data copy(A(1:N, 1:M), A_old(1:N, 1:M))
  CALL INITIALIZE(A)
  CALL INITIALIZE(A_old)
  IF(N .LE. 10 .AND. M .LE. 10) THEN
     WRITE(*, '(A)') "Initial matrix elements"
     CALL DISPLAYMAT(A)
  END IF

  DO iter = 1, num_iterations

     CALL STENCIL_COMPUTE(A, A_old)
     CALL UPDATE_PREVIOUS(A, A_old)

  END DO
  
  IF(N .LE. 10 .AND. M .LE. 10) THEN
     WRITE(*, '(A, I2, A)') "After stencil computations for ", num_iterations, " iterations"
     CALL DISPLAYMAT(A)
  END IF
  !$acc end data
  CALL CPU_TIME(end_t)
  
  DEALLOCATE (A, A_old)

  WRITE(*, '(A, F8.5, A)') "Time measured = ", end_t - start_t, " seconds"
  
CONTAINS

  SUBROUTINE INITIALIZE(matrix)
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER, DIMENSION(:, :) :: matrix(1:N, 1:M)

    !$acc parallel loop collapse(2) present(matrix(1:N, 1:M))
    DO j = 1, M
       DO i = 1, N

          matrix(i, j) = MOD(i+j, 10)

       END DO
    END DO
    !$acc end parallel loop
    
  END SUBROUTINE INITIALIZE

  SUBROUTINE DISPLAYMAT(matrix)
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER, DIMENSION(:, :) :: matrix(1:N, 1:M)

    !$acc serial loop present(matrix(1:N, 1:M))
    DO i = 1, N
       
       WRITE(*, *) (matrix(i, j), j = 1, M)
       
    END DO
    !$acc end serial loop

  END SUBROUTINE DISPLAYMAT

  SUBROUTINE STENCIL_COMPUTE(matA, matA_old)
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER, DIMENSION(:, :) :: matA(1:N, 1:M), matA_old(1:N, 1:M)

    !$acc parallel loop collapse(2) present(matA(1:N, 1:M), matA_old(1:N, 1:M))
    DO j = 1, M-1
       DO i = 2, N

          matA(i, j) = matA_old(i, j) + matA_old(i-1, j) + matA_old(i-1, j+1)

       END DO
    END DO
    !$acc end parallel loop

  END SUBROUTINE STENCIL_COMPUTE

  SUBROUTINE UPDATE_PREVIOUS(matA, matA_old)
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER, DIMENSION(:, :) :: matA(1:N, 1:M), matA_old(1:N, 1:M)

    !$acc parallel loop collapse(2) present(matA(1:N, 1:M), matA_old(1:N, 1:M))
    DO j = 1, M
       DO i = 1, N

          matA_old(i, j) = matA(i, j)

       END DO
    END DO
    !$acc end parallel loop

  END SUBROUTINE UPDATE_PREVIOUS
  
END PROGRAM MAIN

