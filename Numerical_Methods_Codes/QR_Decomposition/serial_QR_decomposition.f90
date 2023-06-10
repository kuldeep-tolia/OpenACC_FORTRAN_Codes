! Serial code for QR decomposition using the Gram-Schmidt method: A = QR

PROGRAM QR_DECOMPOSITION
  IMPLICIT NONE

  INTEGER, PARAMETER :: M = 1500   ! Number of rows of matrix-A
  INTEGER, PARAMETER :: N = 1000   ! Number of columns of matrix A; N < M
  INTEGER            :: i, j, k

  DOUBLE PRECISION, PARAMETER :: eps = 1.0-6 ! A small value to prevent division by zero error
  DOUBLE PRECISION   :: sum_temp
  DOUBLE PRECISION   :: start_time, end_time
  
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: A, Q, R

  ! Allocate memory
  ALLOCATE(A(1:M, 1:N), Q(1:M, 1:N), R(1:N, 1:N))

  ! Initialize input matrix-A
  CALL RANDOM_SEED()
  DO j = 1, N
     DO i = 1, M

        CALL RANDOM_NUMBER(A(i, j))

     END DO
  END DO

  ! A = RESHAPE([1.d0, 1.d0, 1.d0, 1.d0, -1.d0, 4.d0, 4.d0, -1.d0, 4.d0, -2.d0, 2.d0, 0.d0], SHAPE(A))     ! for smaller matrix dimensions and checking the program

  ! Specify the format specifier, update the format specifier replace with the value of N
151 FORMAT(3F10.6)
  

  ! Print the matrix-A for smaller dimensions only
  IF(M .LE. 10 .AND. N .LE. 10) THEN

     WRITE(*, '(A)') "Printing Matrix-A:"
     DO i = 1, M

        WRITE(*, 151) A(i, :)

     END DO
     
  END IF

  CALL CPU_TIME(start_time)

  ! Perform QR decomposition using Gram-Schmidt process
  DO k = 1, N

     DO i = 1, M

        Q(i, k) = A(i, k)
        
     END DO
     
     DO j = 1, k-1
        
        sum_temp = 0.0d0
        DO i = 1, M

           sum_temp = sum_temp + Q(i, j) * A(i, k)
           
        END DO

        R(j, k) = sum_temp

        DO i = 1, M

           Q(i, k) = Q(i, k) - R(j, k) * Q(i, j)
           
        END DO

     END DO

     sum_temp = 0.0d0
     DO i = 1, M

        sum_temp = sum_temp + Q(i, k)**2.0
        
     END DO

     R(k, k) = sum_temp**0.50

     DO i = 1, M

        IF(R(k, k) .GT. eps) Q(i, k) = Q(i, k) / R(k, k)
        
     END DO
     
  END DO
  
  CALL CPU_TIME(end_time)

  ! Print the matrix-Q and matrix-R for smaller dimensions only
  IF(M .LE. 10 .AND. N .LE. 10) THEN

     WRITE(*, '(A)') "Printing Matrix-Q:"
     DO i = 1, M

        WRITE(*, 151) Q(i, :)

     END DO

     WRITE(*, '(A)') "Printing Matrix-R:"
     DO i = 1, N

        WRITE(*, 151) R(i, :)

     END DO
     
  END IF

  WRITE(*, '(A)') "QR decomposition is computed..."
  WRITE(*, '(A, 2I6)') "Matrix dimensions = ", M, N
  WRITE(*, '(A, F15.8, A)') "Compute time for QR decomposition = ", end_time - start_time,  " seconds"
  
  ! Deallocate memory
  DEALLOCATE(A, Q, R)
  
END PROGRAM QR_DECOMPOSITION
