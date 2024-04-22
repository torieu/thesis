from sage.all import *

def Iji_matrix(rows, cols, component_bound):
    return random_matrix(ZZ,rows,cols, algorithm='echelonizable', rank = min(rows, cols), upper_bound=component_bound)
    # n rows, j cols, nejvyšši čisla jsou i

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
    A = Iji_matrix(rows=n, cols=r,  component_bound=q)
    b = matrix(GF(q), [0]*n).transpose()
    A_q = A.change_ring(GF(q))
    b_q = b.change_ring(GF(q))
    K_space = kernel(A_q)
    K_q = matrix(K_space.basis())
    K = K_q.change_ring(ZZ) 
    assert K*A_q == matrix(n-r,r)
    
    return K

print(hard_instance_lattice(n=4, r=1, q=113))

