Problem Description:  

-> In this problem, the objective is to develope the serial, OpenMP, and OpenACC versions of the Jacobi-iterative method to solve a 1D Poisson equation.  
-> To simplify the problem, I have considered a system of linear equations $Ax=b$. The matrix $A$ is a tri-diagonal matrix and the coefficients on the sub, main, and super diagonals are constant and equal to (-1, 2, -1).  
-> The right hand side known vector $b$ is filled with suitable assumed source term and boundary conditions. The user can play with this settings.  
-> The serial,OpenMP-parallelized, and OpenACC-parallelized programs are written and the speedup performance is measured for maximum iterations = $10^5$ with $N = 131072$, where $N$ is the size of the system.  
-> The GPU-accelerated program gave a speedup of about 9x compared to the serial solver and about 4.53x compared to the OpenMP-parallelized solver with 4 threads.
