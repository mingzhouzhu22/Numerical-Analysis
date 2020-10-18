#input:  m -- number of cols, rows
#        d -- column vector of entries on diagonal
#        s -- column vector of entries on subdiagonal and superdiagonal
#output: *R -- Cholesky factor, lower triangular, not return this one but calculated
#        diag -- vector of diagonal entries of R
#        sup  -- vector of superdiagonal entries of R (same as subdiagonal)

def choltri(m,d,s):
    R = np.zeros((m,m))
    A = np.zeros((m,m)) #tridiagonal matrix we want to solve for
    
    for i in range(m):
        A[i][i] = d[i]   #fill out diagonal with values in d
    for i in range(m-1):
        A[i][i+1] = s[i] #fill out superdiagonal
        A[i+1][i] = s[i] #fill out subdiagonal
        
    for i in range(m):
        for k in range(i+1):
            s = sum(R[i][j] * R[k][j] for j in range(k)) #sum of previous terms
            if (i == k): # Diagonal elements, formula 1
                if (A[i][i]-s) > 0: #positive definite:
                    R[i][k] = sqrt(A[i][i]-s)
                else:
                    raise ValueError("given matrix is not positive definite")
            else: #subdiagonal/superdiagonal, same to each other, formula 2
                R[i][k] = (1.0/R[k][k]*(A[i][k]-s))
              
    print(R)
    
    diag = []
    sup  = []
    for i in range(m):
        diag.append(R[i][i])
    for i in range(1,m-1):
        sup.append(R[i][i-1])
        
    return diag,sup
