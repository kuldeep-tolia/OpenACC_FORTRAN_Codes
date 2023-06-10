! Serial code for a 1D Jacobi iterative solver

PROGRAM SERIAL_JACOBI
  IMPLICIT NONE

  INTEGER, PARAMETER :: N = 131072
  INTEGER, PARAMETER :: ITER_MAX = 100000
  INTEGER            :: i, iter

  DOUBLE PRECISION :: dx, comp_time, t1, t2, L, error
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

  CALL CPU_TIME(t1)
  
  DO iter = 1, ITER_MAX

     error = 0.0d0
     
     ! Jacobi iterative solver
     DO i = 1, N-1
        
        phi_new(i) = 0.50d0 * (b(i) + phi(i-1) + phi(i+1))
        error      = MAX(error, ABS(phi_new(i) - phi(i)))

     END DO

     ! Monitors
     IF (MOD(iter, 1000) .EQ. 0) WRITE(*, '(A, I10, A, F15.9)') "Iteration =", iter, " Error =", error
     
     ! Updating the vector
     DO i = 1, N-1

        phi(i) = phi_new(i)

     END DO
     
  END DO
  
  CALL CPU_TIME(t2)
  comp_time = t2 - t1
  
  ! Post-processing
  OPEN(UNIT = 11, FILE = "serial_solution.txt")
  DO i = 0, N
     
     WRITE(11, *) xgrid(i), phi_exact(i), phi(i)

  END DO
  
  WRITE(*, *) "----------------------------"
  WRITE(*, *) "N =", N
  WRITE(*, 11) "Total time taken =", comp_time
  WRITE(*, 11) "Final error =", error
11 format (A, F15.9)
  WRITE(*, *) "Total iterations =", ITER_MAX
  
END PROGRAM SERIAL_JACOBI
