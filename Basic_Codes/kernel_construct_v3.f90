! OpenACC tutorial for kernel construct
! Check how it auto-parallelizes using Minfo
! comments added based on information from -Minfo=accel
  
PROGRAM MAIN
  IMPLICIT NONE

  INTEGER :: a(10)
  INTEGER, VALUE :: i
  
  !$acc kernels                 ! generate a copy of a(:) on GPU
  PRINT*, "One"                 ! run on CPU

  DO i = 1, 10                  ! parallelizable loop, launches 1 gang with 32 vectors
     a(i) = i  
  END DO
  
  PRINT*, "Two: a(10) =", a(10) ! avoids data transfer between CPU-GPU and thus launches a serial kernel(1 gang 32 vectors but only 1 vector lane is operating) on GPU

  DO i = 1, 10                  ! parallelizable loop, launches 1 gang with 32 vectors
     a(i) = a(i) * 2
  END DO

  PRINT*, "Three"               ! Launched a serial kernel on GPU
  !$acc end kernels

END PROGRAM MAIN

