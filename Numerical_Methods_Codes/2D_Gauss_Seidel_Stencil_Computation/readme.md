Problem Description:  

-> In this problem, the objective is to develope the OpenACC version of the Gauss-Seidel iterative method to perform some 2D stencil computation.  
-> The Jacobi-iterative method uses the older values at a given instant. Thus, there is no loop-carried dependency and it is straight forward to parallelize the Jacobi method solver.  
-> The Gauss-Seidel iterative method uses the updated/latest values at a given instant which creates a loop-carried dependency. So it cannot be parallelized directly.  
-> Consider the following discrete stencil computation:  
$$A_{i,j} = f \left( A_{i+1,j}, A_{i-1,j}, A_{i,j+1}, A_{i,j-1} \right)$$    
-> The matrix-A is initialized using random floating point values.  
-> For simplicity, the boundary points are avoided in the stencil computation.  
-> As the Gauss-Seidel method cannot be parallelized directly, it is parallelized using the Red-Black coloring approach.  
-> In this method, each alternate points are colored as red and black. So in the $1^{st}$ pass, we only traverse through red points and in the next pass, we traverse through only black points. With this type of traversing method, the loop-carried dependency will be gone.  
