Problem Description:  
  
-> A stencil computation is implemented wherein each iteration, cells compute their values using the neighboring ones, according to the following formula:  
$$A[i][j] += A[i-1][j] + A[i-1][j+1]$$  

-> The computation boundaries remain intact.  
-> The stencil computation is performed for $10$ outer iterations.  
-> The matrix $A$ is initialized using random integer values.  
-> Another version of the same stencil computation is developed, where the in-situ update is relaxed. Both the program versions are compared with the serial solver.
