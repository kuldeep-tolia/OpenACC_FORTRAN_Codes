! OpenMP parallelized code for a 1D Jacobi iterative solver

PROGRAM OPENMP_JACOBI
  USE OMP_LIB
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 131072
  INTEGER, PARAMETER :: ITER_MAX = 100000
  INTEGER            :: i, iter
  INTEGER            :: num_thread

  CHARACTER(100)     :: numchar, nameq
  
  DOUBLE PRECISION :: dx, comp_time, t1, t2, L, error
  DOUBLE PRECISION, DIMENSION(0:N) :: xgrid, phi, phi_new, phi_exact, b
  
  L  = 1.0d0
  dx = L / N
  
  IF(COMMAND_ARGUMENT_COUNT() .NE. 1) THEN
     
     WRITE(*, *) "Command line argument of number of threads is required."
     STOP
     
  END IF
  
  CALL GET_COMMAND_ARGUMENT(0, nameq)	
  CALL GET_COMMAND_ARGUMENT(1, numchar)
  READ(numchar, *) num_thread
  
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

  t1 = OMP_GET_WTIME()
  
  !$OMP PARALLEL NUM_THREADS(num_thread) DEFAULT(shared) PRIVATE(i)
  DO iter = 1, ITER_MAX

     error = 0.0d0
     
     ! Jacobi iterative solver
     !$OMP DO REDUCTION(MAX:error)
     DO i = 1, N-1
        
        phi_new(i) = 0.50d0 * (b(i) + phi(i-1) + phi(i+1))
        error      = MAX(error, ABS(phi_new(i) - phi(i)))

     END DO
     !$OMP END DO

     ! Monitors
     !$OMP SINGLE
     IF (MOD(iter, 1000) .EQ. 0) WRITE(*, '(A, I10, A, F15.9)') "Iteration =", iter, " Error =", error
     !$OMP END SINGLE
     
     ! Updating the vector
     !$OMP DO
     DO i = 1, N-1

        phi(i) = phi_new(i)

     END DO
     !$OMP END DO
     
  END DO
  !$OMP END PARALLEL
  
  t2 = OMP_GET_WTIME()
  comp_time = t2 - t1
  
  ! Post-processing
  OPEN(UNIT = 11, FILE = "openmp_solution.txt")
  DO i = 0, N
     
     WRITE(11, *) xgrid(i), phi_exact(i), phi(i)

  END DO
  
  WRITE(*, *) "----------------------------"
  WRITE(*, *) "N =", N
  WRITE(*, 11) "Total time taken =", comp_time
  WRITE(*, 11) "Final error =", error
11 format (A, F15.9)
  WRITE(*, *) "Total iterations =", ITER_MAX
  
END PROGRAM OPENMP_JACOBI
