! OpenACC tutorial for kernel construct
! Check how it auto-parallelizes using Minfo
! comments added based on information from -Minfo=accel
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER :: a(10)
  INTEGER, VALUE :: i
  INTEGER :: sum

  a = 0
  
  !$acc kernels async          ! generate a copy of a(:) on GPU
  sum = 0                      ! shared variable amongst thread is created
  PRINT*, "One"                ! run on CPU

  DO i = 1, 10                 ! parallelizable loop, launches 1 gang with 32 vectors, also performing reduction (+:sum) operation on sum variable
     sum = sum + a(i) + i    
  END DO
  
  PRINT*, "Two: ", sum         ! avoids data transfer between CPU-GPU and thus launches a serial kernel(1 gang 32 vectors but only 1 vector lane is operating) on GPU

  DO i = 1, 10                 ! parallelizable loop, launches 1 gang with 32 vectors
     a(i) = a(i) + sum
  END DO

  PRINT*, "Three"              ! serial kernel launched on GPU 

  DO i = 1, 10                 ! serial kernel launched on GPU
     WRITE(*, *) i, a(i)
  END DO

  PRINT*, "Four"               ! serial kernel launched on GPU
  !$acc end kernels

  PRINT*, "Five"               ! run on CPU

END PROGRAM MAIN

