! Serial code for matrix addition
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 300   ! set number of rows
  INTEGER, PARAMETER :: M = 300   ! set number of columns
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: A, B, C

  DOUBLE PRECISION :: start_t, end_t

  ALLOCATE(A(1:N, 1:M), B(1:N, 1:M), C(1:N, 1:M))

  CALL CPU_TIME(start_t)

  CALL INITIALIZE()

  IF(N .LE. 10 .AND. M .LE. 10) THEN

     WRITE(*, '(A)') "Matrix-A"
     CALL DISPLAYMAT(A)
     WRITE(*, '(A)') "Matrix-B"
     CALL DISPLAYMAT(B)

  END IF

  CALL MATADD()
  
  IF(N .LE. 10 .AND. M .LE. 10) THEN

     WRITE(*, '(A)') "Matrix-C"
     CALL DISPLAYMAT(C)

  END IF
  
  CALL CPU_TIME(end_t)

  WRITE(*, '(A, I8)') "Sanity check: C(N/2)(M/2) = ", C(N/2, M/2)
  
  DEALLOCATE(A, B, C)

  WRITE(*, '(A, F8.5, A)') "Time measured = ", end_t - start_t, " seconds"
  
CONTAINS

  SUBROUTINE INITIALIZE()
    IMPLICIT NONE

    INTEGER :: i, j

    DO j = 1, M
       DO i = 1, N

          A(i, j) = i * N + j
          B(i, j) = i * N + j

       END DO
    END DO

  END SUBROUTINE INITIALIZE

  SUBROUTINE DISPLAYMAT(matrix)
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER, DIMENSION(1:N, 1:M) :: matrix
    
    DO i = 1, N

       WRITE(*, *) (matrix(i, j), j = 1, M)

    END DO

  END SUBROUTINE DISPLAYMAT

  SUBROUTINE MATADD()
    IMPLICIT NONE

    INTEGER :: i, j

    DO j = 1, M
       DO i = 1, N

          C(i, j) = A(i, j) + B(i, j)

       END DO
    END DO
        
  END SUBROUTINE MATADD
    
END PROGRAM MAIN
