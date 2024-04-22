
from sage.all import *
from generate import random_special_matrix
from sage.modules.misc import gram_schmidt
from sage.modules.free_module_element import vector
from sage.rings.rational_field import QQ

def mu(b_i, b_j):
    return b_i.dot_product(b_j) / b_j.dot_product(b_j) if b_j.dot_product(b_j) != 0 else 0

def full_LLL(basis, delta=0.75):
    n = basis.nrows()
    basis = [vector(basis[i]) for i in range(n)]
    while_it = 0

    # Initially orthogonalizing the basis
    ortho, _ = gram_schmidt(basis)  # Ensure the orthogonalized vectors are returned
    ortho = [vector(QQ, v) for v in ortho]  # Convert tuples to vectors in QQ
    
    k = 1
    while k < n:
        while_it += 1
        for j in range(k - 1, -1, -1):
            proj = mu(basis[k], ortho[j])
            if abs(proj) > 1/2:
                basis[k] = basis[k] - basis[j] * round(proj)
                ortho, _ = gram_schmidt(basis)  # Re-orthogonalize after modification
                ortho = [vector(QQ, v) for v in ortho]  # Ensure conversion to QQ

        if ortho[k].dot_product(ortho[k]) >= (delta - mu(basis[k], ortho[k-1])**2) * ortho[k-1].dot_product(ortho[k-1]):
            k += 1
        else:
            basis[k], basis[k-1] = basis[k-1], basis[k]
            ortho, _ = gram_schmidt(basis)  # Orthogonalize after swapping
            ortho = [vector(QQ, v) for v in ortho]  # Convert to rational vectors

            k = max(k - 1, 1)
    basis = matrix(basis)
    return basis, while_it
