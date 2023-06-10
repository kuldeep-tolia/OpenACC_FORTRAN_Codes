-> This section discusses about basic OpenACC clauses, constructs like parallel construct, kernels construct, async clause, num_gangs clause.    
-> The computer programs only demonstrates the behaviour of a particular OpenACC clause/construct. To check the syntax/operation of a particular OpenACC clause/construct, you may refer to the following OpenACC website:  
- https://www.openacc.org/  

-> To compile and run the codes, I have used PGI compiler.  
-> The program is profiled using the environment variable PGI_ACC_TIME.  
-> Compiling and running a FORTRAN program:
- $ pgf90 -acc -ta=nvidia -Minfo=accel file_name.f90 -o ./output_name.out  
- export PGI_ACC_TIME=1                          (Note: To turn on the environment variable for profiling purposes)  
- $ ./output_name.out
