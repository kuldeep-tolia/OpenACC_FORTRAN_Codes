! OpenACC parallelized program to calculate first derivative for Pade scheme using LU-decomposition

PROGRAM OPENACC_LUD
  USE OPENACC
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 5001  ! Number of grid points
  INTEGER            :: i, j

  DOUBLE PRECISION, PARAMETER :: Lx  = 3.0d0 ! Length of domain
  DOUBLE PRECISION, PARAMETER :: eps = 1.0d-6 ! define the tolerance level to avoid division by zero error
  DOUBLE PRECISION            :: h, start_time, end_time

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE    :: xgrid, df, df_exact, rhs, f, y
  DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: A

  h = Lx / (N - 1.0d0)			! grid spacing

  ! Allocate memory
  ALLOCATE(xgrid(1:N), f(1:N), df(1:N), y(1:N), df_exact(1:N), rhs(1:N))
  ALLOCATE(A(1:N, 1:N))

  ! Initialize vectors and matrix
  DO i = 1, N

     xgrid(i)    = (i - 1.0d0) * h
     f(i)        = SIN(5.0d0 * xgrid(i))
     df_exact(i) = 5.0d0 * COS(5.0d0 * xgrid(i))
     df(i)       = 0.0d0
     y(i)        = 0.0d0
     
  END DO

  A = 0.0d0

  CALL ACC_INIT(acc_device_nvidia)
  CALL CPU_TIME(start_time)

  !$ACC DATA COPYIN(A, rhs, f, h, df, y) COPYOUT(df)
  CALL INIT_MAT(A)
  CALL INIT_VEC(rhs, f, h)
  CALL LU_FACTOR(A)
  CALL DECOMPOSE(A, rhs, df, y)
  !$ACC END DATA
  
  CALL CPU_TIME(end_time)

  ! Post-processing
  OPEN(11, file = "openacc_derivative.txt")
  DO i = 0, N

     WRITE(11, *) xgrid(i), df_exact(i), df(i)

  END DO
  
  WRITE(*, '(A, I6)') "Problem size = ", N
  WRITE(*, '(A, F15.8, A)') "Compute time taken = ", end_time - start_time, " seconds"
  WRITE(*, *) "---------- Program ran successfully. Exiting!!! ----------"
  
  ! Deallocate memory
  DEALLOCATE(xgrid, f, df, y, df_exact, rhs)
  DEALLOCATE(A)

CONTAINS

  SUBROUTINE INIT_MAT(mat)
    IMPLICIT NONE

    INTEGER :: i, j

    DOUBLE PRECISION, DIMENSION(1:N, 1:N) :: mat

    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(mat) NUM_GANGS(10)
    DO j = 1, N
       DO i = 1, N

          IF(i .EQ. j .AND. i .EQ. 1) THEN

             mat(i, j)   = 1.0d0
             mat(i, j+1) = 2.0d0

          ELSE IF(i .EQ. j .AND. i .GT. 1 .AND. i .LT. N) THEN

             mat(i, j-1) = 1.0d0
             mat(i, j)   = 4.0d0
             mat(i, j+1) = 1.0d0

          ELSE IF(i .EQ. j .AND. i .EQ. N) THEN

             mat(i, j)   = 1.0d0
             mat(i, j-1) = 2.0d0

          END IF

       END DO
    END DO
    !$ACC END PARALLEL LOOP
    
  END SUBROUTINE INIT_MAT

  SUBROUTINE INIT_VEC(b, fvar, dx)
    IMPLICIT NONE

    INTEGER :: i

    DOUBLE PRECISION                 :: dx
    DOUBLE PRECISION, DIMENSION(1:N) :: b, fvar

    !$ACC PARALLEL LOOP PRESENT(b, fvar, dx) NUM_GANGS(10)
    DO i = 1, N

       IF(i .EQ. 1) THEN

          b(i) = (1.0d0 / dx) * (-2.50d0 * fvar(i) + 2.0d0 * fvar(i+1) + 0.50d0 * fvar(i+2))

       ELSE IF (i .EQ. N) THEN

          b(i) = (1.0d0 / dx) * ( 2.50d0 * fvar(i) - 2.0d0 * fvar(i-1) - 0.50d0 * fvar(i-2))

       ELSE

          b(i) = (3.0d0 / dx) * (fvar(i+1) - fvar(i-1))

       END IF

    END DO
    !$ACC END PARALLEL LOOP
    
  END SUBROUTINE INIT_VEC
  
  SUBROUTINE LU_FACTOR(mat)
    IMPLICIT NONE

    INTEGER :: i, j, k

    DOUBLE PRECISION, DIMENSION(1:N, 1:N) :: mat

    !$ACC KERNELS PRESENT(mat)
    DO j = 1, N-1

       !$ACC LOOP
       DO i = j+1, N
          
          mat(i, j) = mat(i, j) / mat(j, j)

          !$ACC LOOP
          DO k = j+1, N
             
             mat(i, k) = mat(i, k) - mat(i, j) * mat(j, k)

          END DO
          !$ACC END LOOP

       END DO
       !$ACC END LOOP
       
    END DO
    !$ACC END KERNELS
    
  END SUBROUTINE LU_FACTOR

  SUBROUTINE DECOMPOSE(mat, vec, v1, v2)
    IMPLICIT NONE

    ! mat --> input LU coefficient matrix
    ! vec --> rhs vector
    ! v1  --> solution vector: x
    ! v2  --> intermediate vector: y
    
    INTEGER          :: i, j

    DOUBLE PRECISION :: w_temp
    
    DOUBLE PRECISION, DIMENSION(1:N, 1:N) :: mat
    DOUBLE PRECISION, DIMENSION(1:N)      :: vec, v1, v2

    !$ACC KERNELS PRESENT(vec, mat, v2)
    DO i = 1, N                 ! Computing L*v2=vec

       w_temp = vec(i)

       !$ACC LOOP REDUCTION(+:w_temp)
       DO j = 1, i

          w_temp = w_temp + -mat(i, j) * v2(j)
          
       END DO
       !$ACC END LOOP

       v2(i) = w_temp
       
    END DO
    !$ACC END KERNELS

    !$ACC KERNELS PRESENT(vec, mat, v2, v1)
    DO i = 1, N                 ! Computing U*v1=v2

       w_temp = v2(N+1-i)

       !$ACC LOOP REDUCTION(+:w_temp)
       DO j = N+1-i, N

          w_temp = w_temp + -mat(N+1-i, j) * v1(j)
          
       END DO
       !$ACC END LOOP

       IF(A(N+1-i, N+1-i) .GE. eps) THEN

          v1(N+1-i) = w_temp / mat(N+1-i, N+1-i)

       ELSE

          v1(N+1-i) = w_temp

       END IF
       
    END DO
    !$ACC END KERNELS
    
  END SUBROUTINE DECOMPOSE
  
END PROGRAM OPENACC_LUD
