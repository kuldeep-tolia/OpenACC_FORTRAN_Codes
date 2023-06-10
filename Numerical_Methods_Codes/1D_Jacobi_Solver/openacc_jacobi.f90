! OpenACC parallelized code for a 1D Jacobi iterative solver

PROGRAM OPENACC_JACOBI
  USE OPENACC
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 131072
  INTEGER, PARAMETER :: ITER_MAX = 100000
  INTEGER            :: i, iter

  DOUBLE PRECISION :: dx, comp_time, t1, t2, L, error, diff
  DOUBLE PRECISION, DIMENSION(0:N) :: xgrid, phi, phi_new, phi_exact, b

  L  = 1.0d0
  dx = L / N
  
  ! Initialize vectors
  DO i = 0, N

     xgrid(i)     = i * dx
     phi_exact(i) = ((-4.0d0 * xgrid(i)**3) + 1204.0d0 * xgrid(i) + 300.0d0) / 3.0d0
     b(i)         = dx * dx * 8.0d0 * xgrid(i)
     phi(i)       = 100.0d0
     phi_new(i)   = phi(i)
     
  END DO
  
  ! Setting boundary conditions
  phi(0)     = 100.0d0
  phi(N)     = 500.0d0
  phi_new(0) = phi(0)
  phi_new(N) = phi(N)

  CALL ACC_INIT(acc_device_nvidia)

  CALL CPU_TIME(t1)

  !$ACC DATA COPY(phi, b, phi_new) CREATE(diff)
  DO iter = 1, ITER_MAX

     error = 0.0d0
     
     ! Jacobi iterative solver
     !$ACC PARALLEL LOOP PRESENT(phi, b, phi_new)
     DO i = 1, N-1
        
        phi_new(i) = 0.50d0 * (b(i) + phi(i-1) + phi(i+1))

     END DO
     !$ACC END PARALLEL LOOP

     ! Monitors
     IF (MOD(iter, 10000) .EQ. 0) THEN

        !$ACC PARALLEL LOOP PRESENT(phi, phi_new, diff) REDUCTION(MAX:error) PRIVATE(diff)
        DO i = 1, N-1

           diff = phi_new(i) - phi(i)
           IF(diff .LT. 0.0d0) diff = -diff
           error = MAX(error, diff)
           
        END DO
        !$ACC END PARALLEL LOOP

        WRITE(*, '(A, I10, A, F15.9)') "Iteration =", iter, " Error =", error
        
     END IF
     
     ! Updating the vector   
     !$ACC PARALLEL LOOP PRESENT(phi, phi_new)
     DO i = 1, N-1
        
        phi(i) = phi_new(i)
        
     END DO
     !$ACC END PARALLEL LOOP
     
  END DO
  !$ACC END DATA
  
  CALL CPU_TIME(t2)
  comp_time = t2 - t1
  
  ! Post-processing
  OPEN(UNIT = 11, FILE = "openacc_solution.txt")
  DO i = 0, N
     
     WRITE(11, *) xgrid(i), phi_exact(i), phi(i)

  END DO
  
  WRITE(*, *) "----------------------------"
  WRITE(*, *) "N =", N
  WRITE(*, 11) "Total time taken =", comp_time
  WRITE(*, 11) "Final error =", error
11 format (A, F15.9)
  WRITE(*, *) "Total iterations =", ITER_MAX
  
END PROGRAM OPENACC_JACOBI
