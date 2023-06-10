Problem Description:  

-> The tri-diagonal matrix system obtained from, say, using implicit finite-difference derivative schemes (like Pade's scheme) can be solved using LU-decomposition method.  
-> This is an OpenACC program to compute the LU decomposition, and subsequently solve the system of linear equations using forward/backward substitutions.  
-> Given a matrix $A[N][N]$ is decomposed into a lower and upper triangular matrix such that $A = L * U$.  
-> Consider the calculation of the derivative of the following function,  
$$f(x) = sin(5x), \hspace{4mm} 0 \leq x \leq 5$$  
using $4^{th}$ order accurate Pade scheme for the interior and $3^{rd}$ order accurate one-sided Pade scheme near the boundaries.  
-> Pade's scheme is given as:  
$$g_0 + 2g_1 = \frac{1}{h} \left( -\frac{5}{2}f_0 + 2f_1 + \frac{1}{2}f_2 \right)$$  

$$g_{j+1} + 4g_j + g_{j-1} = \frac{3}{h} \left( f_{j+1} - f_{j-1} \right)$$  

$$g_n + 2g_{n-1} = \frac{1}{h} \left( \frac{5}{2}f_n - 2f_{n-1} - \frac{1}{2}f_{n-2} \right)$$  

where $g(x) = f'(x)$, $h$ is the grid spacing, $n$ is the number of grid points in $x$ direction and $j = 1,2,...,n-1$.  
-> The analytical solution is plotted with the obtained numerical solution.  
-> Speedup plots are given for $n = 1000$ and $n = 5000$ to observe the parallel performance of the OpenACC-parallelized LU-decomposition method.  
-> A speedup of about 2.5x is obtained with $n = 1000$ and about 4.23x is obtained with $n = 5000$ when compared with the serial solver.  
