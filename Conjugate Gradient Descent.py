def conj_grad(n):
    nodes = [] #list of nodes
    f = [] #list of f
    N = n
    h = 1/N
    iter = 0

    #nodes.append(0)
    for j in range(1,N):
        nodes.append(0+j*h)
    #nodes.append(1)

    for i in nodes:
        f.append(2*(pi**2)*sin(pi*i))
    #f[0] = 0
    #f[N] = 0

    A = np.zeros((N-1,N-1)) #create a n x n matrix with 0 as all entries
    for i in range(N-1): #fill out diagonal entries
        A[i,i] = 2 + (pi**2)*(h**2)
    for j in range(N-2): #fill out entries below diagonal
        A[j,j+1] = -1
    for k in range(1,N-1): #fill out entries above diagonal
        A[k,k-1] = -1
    A = A/(h**2)

    b = np.array(f)
    b.shape = (N-1,1)

    xi = np.ones((N-1,1)) #initial guess, x^0

    xlist = []
    tlist = []
    rlist = []
    slist = []
    vlist = []

    xlist.append(xi)
    rlist.append(np.subtract(b,np.dot(A,xlist[0]))) #r0
    vlist.append(b-np.dot(A,xlist[0])) #v0 = r0
    tlist.append((np.dot(np.transpose(rlist[0]),rlist[0]))/(np.dot(np.transpose(vlist[0]),np.dot(A,vlist[0]))))

    while (sqrt(np.sum(rlist[-1]**2))) > tol:
        xlist.append(xlist[-1]+tlist[-1]*vlist[-1])
        rlist.append(rlist[-1]-tlist[-1]*(np.dot(A,vlist[-1])))
        slist.append((np.dot(np.transpose(rlist[-1]),rlist[-1]))/(np.dot(np.transpose(rlist[-2]),rlist[-2])))
        vlist.append(rlist[-1]+slist[-1]*vlist[-1])
        iter+=1

    approx = np.transpose(xlist[-1])

    exact = []

    for x in nodes:
        exact.append(sin(pi*x))

    exact = np.array(exact)

    error = exact - approx
    error = error**2
    error = np.sum(error)
    error = sqrt(error)
    print(error)
    print("Required iteration is",iter)
    
conj_grad(100)