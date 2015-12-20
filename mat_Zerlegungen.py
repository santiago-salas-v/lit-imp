import numpy as np
from sympy import pprint, Matrix, nsimplify

def lrpd(A):
    """L,R,P,D,DA aus der LR = PDA Zerlegung
    Methode aus Dahmen W., Reusken A.; Numerik fuer Ingenieure und Naturwissenschaeftler; Springer [DE] S. 79
    :param A numpy.matrix NxN
    """
    n = A.shape[0]
    p = np.matrix(np.eye(n).astype(dtype=float))
    d = np.matrix(np.eye(n).astype(dtype=float))
    da = np.matrix(np.zeros_like(p))
    l = np.matrix(np.eye(n))
    r = np.matrix(np.zeros_like(p))
    indexes_r = range(n)
    max_aij = 0.0
    for i in range(n):
        d[i,i] = 1/abs(A[i]).sum()
        # Skalierung
        da[i] = d[i,i]*A[i]
    for j in range(n):
        max_aij = max(abs(da[j:,j]))
        indexes_r[j] = j + abs(da[j:,j]).argmax()
        # Spaltenpivotisierung
        da[[j,indexes_r[j]]] = da[[indexes_r[j],j]]
        p[[j,indexes_r[j]]] = p[[indexes_r[j],j]]
        for i in range(j+1,n):
            # neue Eintraege in L
            da[i,j] = da[i,j]/da[j,j]
            for k in range(j+1,n):
                # neue Eintraege in R
                da[i,k] = da[i,k]-da[i,j]*da[j,k]
    for i in range(n):
        l[i+1:n,i] = da[i+1:n,i]
        r[i,i:n] = da[i,i:n]
    return l,r,p,d,da

A = np.matrix([[1,5,0],[2,2,2],[-2,0,2]]).astype(dtype = float)
A = np.matrix([[2,-1,-3,3],[4,0,-3,1],[6,1,-1,6],[-2,-5,4,1]]).astype(dtype = float)
l,r,p,d,da = lrpd(A)
for mat in [l,r,p,d,da,l*r]:
    mat = Matrix(mat)
    for i in range(mat.shape[0]):
        mat[i,:] = Matrix(map(nsimplify,mat[i,:])).T
    pprint(mat)
    print "\n"