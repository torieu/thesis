from sage.all import *

def Iji_matrix(n, j, i):
    return random_matrix(ZZ,n,j, algorithm='echelonizable', rank = min(n, j), upper_bound=i)

def hard_instance_lattice(n, r, q):
    '''
    implemented according to page 2. of https://people.csail.mit.edu/vinodv/CS294/ajtai99.pdf
        part: 
        Assume that  
        ...
        The distribution of L 
    '''
    assert  r < n
    q = Integer(q)
    assert q.is_prime() == True
    A = Iji_matrix(n=n, j=r,  i=q)
    b = matrix(GF(q), [0]*n).transpose()
    A_q = A.change_ring(GF(q))
    b_q = b.change_ring(GF(q))
    K_space = kernel(A_q)
    K_q = matrix(K_space.basis())
    K = K_q.change_ring(ZZ) 
    assert K*A_q == matrix(n-r,r)
    
    return K

hard_instance_lattice(n=4, r=1, q=113)


