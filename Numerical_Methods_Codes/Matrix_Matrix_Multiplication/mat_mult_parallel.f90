! Paralle code using OpenACC for square matrix multiplication
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 50   ! define the size of matrix
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: A, B, C

  DOUBLE PRECISION :: start_t, end_t

  ALLOCATE(A(1:N, 1:N), B(1:N, 1:N), C(1:N, 1:N))

  CALL CPU_TIME(start_t)
  
  !$ACC DATA CREATE(A(1:N, 1:N), B(1:N, 1:N), C(1:N, 1:N)) COPYOUT(C(1:N, 1:N))

  CALL INITIALIZE()

  !$ACC KERNELS
  C(:, :) = 0
  !$ACC END KERNELS

  IF(N .LE. 10) THEN

     WRITE(*, '(A)') "Matrix-A"
     CALL DISPLAYMAT(A)
     WRITE(*, '(A)') "Matrix-B"
     CALL DISPLAYMAT(B)

  END IF

  CALL MATMULT()
  
  IF(N .LE. 10) THEN

     WRITE(*, '(A)') "Matrix-C"
     CALL DISPLAYMAT(C)

  END IF
  
  !$ACC END DATA
  
  CALL CPU_TIME(end_t)

  WRITE(*, '(A, I8)') "Sanity check: C(1)(1) = ", C(1, 1)
  
  DEALLOCATE(A, B, C)

  WRITE(*, '(A, F8.5, A)') "Time measured = ", end_t - start_t, " seconds"
  
CONTAINS

  SUBROUTINE INITIALIZE()
    IMPLICIT NONE

    INTEGER :: i, j

    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(A(1:N, 1:N), B(1:N, 1:N))
    DO j = 1, N
       DO i = 1, N

          A(i, j) = i + j
          B(i, j) = i + j

       END DO
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE INITIALIZE

  SUBROUTINE DISPLAYMAT(matrix)
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER, DIMENSION(1:N, 1:N) :: matrix

    !$ACC SERIAL LOOP PRESENT(matrix(1:N, 1:N))
    DO i = 1, N

       WRITE(*, *) (matrix(i, j), j = 1, N)

    END DO
    !$ACC END SERIAL LOOP

  END SUBROUTINE DISPLAYMAT

  SUBROUTINE MATMULT()
    IMPLICIT NONE

    INTEGER :: i, j, k

    !$ACC PARALLEL LOOP COLLAPSE(3) PRESENT(A(1:N, 1:N), B(1:N, 1:N), C(1:N, 1:N)) REDUCTION(+:C(1:N, 1:N))

    DO j = 1, N
       DO k = 1, N
          DO i = 1, N

             C(i, j) = C(i, j) + A(i, k) * B(k, j)
             
          END DO
       END DO
    END DO
    !$ACC END PARALLEL LOOP
        
  END SUBROUTINE MATMULT
    
END PROGRAM MAIN
