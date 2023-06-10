! Paralle code using OpenACC for red-black mapping
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 100   ! define the size of matrix
  INTEGER            :: iter
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: A
  DOUBLE PRECISION                               :: error

  DOUBLE PRECISION :: start_t, end_t

  ALLOCATE(A(1:N, 1:N))

  CALL CPU_TIME(start_t)
  
  !$ACC DATA CREATE(A(1:N, 1:N))

  CALL INITIALIZE()

  IF(N .LE. 10) THEN

     WRITE(*, '(A)') "Matrix-A"
     CALL DISPLAYMAT(A)

  END IF

  DO iter = 1, 100

     error = 0.0d0
     CALL REDBLACK()

  END DO
  
  !$ACC END DATA
  
  CALL CPU_TIME(end_t)

  WRITE(*, '(A, I8, A, F10.7)') "Sanity check: iter = ", iter, " error = ", error
  
  DEALLOCATE(A)

  WRITE(*, '(A, F8.5, A)') "Time measured = ", end_t - start_t, " seconds"
  
CONTAINS

  SUBROUTINE INITIALIZE()
    IMPLICIT NONE

    INTEGER :: i, j

    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(A(1:N, 1:N))
    DO j = 1, N
       DO i = 1, N

          A(i, j) = i * 0.50d0 + j * 0.30d0

       END DO
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE INITIALIZE

  SUBROUTINE DISPLAYMAT(matrix)
    IMPLICIT NONE

    INTEGER :: i, j
    DOUBLE PRECISION, DIMENSION(1:N, 1:N) :: matrix

    !$ACC SERIAL LOOP PRESENT(matrix(1:N, 1:N))
    DO i = 1, N

       WRITE(*, *) (matrix(i, j), j = 1, N)

    END DO
    !$ACC END SERIAL LOOP

  END SUBROUTINE DISPLAYMAT

  SUBROUTINE REDBLACK()
    IMPLICIT NONE

    INTEGER :: i, j
    DOUBLE PRECISION :: temp

    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(A(1:N, 1:N)) COPYIN(error, temp)

    DO j = 2, N-1
       DO i = 2, N-1

          IF(MOD(i+j, 2) .EQ. 0) THEN

             temp = A(i, j)
             A(i, j) = 0.250d0 * (A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1))
             error = error + (A(i, j) - temp)
             
          END IF
          
       END DO
    END DO

    !$ACC END PARALLEL LOOP

    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(A(1:N, 1:N), error) COPYOUT(error)

    DO j = 2, N-1
       DO i = 2, N-1

          IF(MOD(i+j, 2) .EQ. 1) THEN

             temp = A(i, j)
             A(i, j) = 0.250d0 * (A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1))
             error = error + (A(i, j) - temp)
             
          END IF
          
       END DO
    END DO

    !$ACC END PARALLEL LOOP
        
  END SUBROUTINE REDBLACK
    
END PROGRAM MAIN
