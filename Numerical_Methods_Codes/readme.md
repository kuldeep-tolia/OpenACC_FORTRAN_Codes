
-> This section has computer programs that implements solution procedures for solving linear and non-linear algebraic equations, ordinary differential equations (ODE) and partial differential equations (PDE).  
-> Finite Difference Method (FDM) forms the basis of the discretization of the governing equations.  
-> The programs are parallelized using OpenACC and the speedup performance of the parallelized codes is measured.    
-> To check the syntax/operation of a particular OpenACC clause/construct, you may refer to the following website:  
- https://www.openacc.org/  

-> To compile and run the codes, I have used PGI compiler.  
-> The program is profiled using the environment variable PGI_ACC_TIME.  
-> Compiling and running a FORTRAN program:
- $ pgf90 -acc -ta=nvidia -Minfo=accel file_name.f90 -o ./output_name.out  
- export PGI_ACC_TIME=1                          (Note: To turn on the environment variable for profiling purposes)  
- $ ./output_name.out
