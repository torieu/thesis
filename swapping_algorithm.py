#!/usr/bin/env python
# coding: utf-8

# # vector swap
# 
# 

# In[224]:


from generate import *
from closest_vs_minimal import *

# def gram_matrix(matrix) -> matrix:
#     return matrix * matrix.transpose()

# def vector_to_list(vector) -> list:
#     return [float(num) for num in vector.list()]

def find_real_minimum_genfile(G, current_row, inserted_component) -> vector:
    """
    returns linear combination of the minimum over real numbers
    the >current_row< parameter is the coordinate of fixed col/row idk
    the >K< parameter is the number inserted into this column
    """
    G = matrix_to_list(G)
    dimension=len(G)
    matrixA = matrix(dimension - 1, dimension - 1, 0) # square matrix of size (dimension - 1) x (dimension - 1), filled with zeros.
    matrixB = matrix(dimension - 1, 1, 0) # column matrix of size (dimension - 1) x 1, filled with zeros.
    # matrixA[0,0] = 1
    # print("G",matrix(G))
    a, b = 0, 0
    for row in range(dimension):
        if row != current_row:
            matrixA[a] = [G[row][j] for j in range(len(G[row])) if j != current_row]
            matrixB[b] = sum([inserted_component * G[row][j] for j in range(len(G[row])) if j == current_row])
            a += 1
            b += 1
    # insert indices
    # print("A", matrixA,"\n","b", -(1)*matrixB, )
    result = vector_to_list(matrixA.solve_right((-1)*(inserted_component) * matrixB))
    result.insert(current_row-1, inserted_component)
    print(result)
    return vector(result).n(digits=5)


# In[225]:


# %run closest_versus_minimal_tools.ipynb
# def swapping_alg_ver1(dimension, perimeter):
#     """ Verse s for cyklem """
#     basis = randomMatrix(dimension, perimeter)
#     basis.sort by length

#     for i in range(dimension):
#         perform a cube with fixed index i
#         realmin = vysledek cube tedy něco jkao [1, 0.54, 7.2345]
#         sv = zaokrouhleni realmin * B
#         basis[i] = sv




def swapping(basis_as_matrix):
    # basis_as_matrix= matrix([[5, 3, 4], [-2, 4, 4], [-10, -3, -7]]
    # basis_as_matrix = matrix([[3, 10, 4], [-4, 2, 7], [-4, -8, -3]]) # asi funguje??
    # basis_as_matrix = matrix([[-8, -9, 7], [10, 9, -6], [-8, 4, 10]]) # asi funguje??
    # print(basis_as_matrix)
    basis_as_list= matrix_to_list(basis_as_matrix)
    dim = len(basis_as_list)
    print("B = ",basis_as_list)
    used_vectors =[]
    used_vectors.extend(basis_as_list)

    G = gram_matrix(basis_as_matrix)

    basis_as_list.sort(key=lambda x: euclidean_norm(x))
    print("Basis na začatku: ")
    for i in range(dim):
            print("  b_",i, " ", euclidean_norm(basis_as_list[i]), basis_as_list[i])
    print("--------------")
    basis_as_matrix=matrix(basis_as_list)
    # print(basis_as_matrix)
    # print(basis_as_list)

    i = 0
    while i < dim:
#     for j in range(0,dimension):
        # print("\n")
        # print(i)
        real_min = find_real_minimum_genfile(G, i, 1)

        # print(evaluate_norms_at_cube(cube_points(real_min), basis_as_list))
        
        shortestLC = shortest_lc_in_cube(real_min, basis_as_list)
        # print(real_min)
        # print(shortestLC)
        shortest_vector = vector(shortestLC)*basis_as_matrix
        print("sv", shortest_vector,(shortest_vector).norm().n())
        # print("used vectors",used_vectors)
        if vector_to_list(shortest_vector) in used_vectors:
              print("ERROR: I already used the vector ", shortest_vector)
              i+=1
              continue
        if shortest_vector.norm() > basis_as_matrix[-1].norm():
              print("ERROR, LONGER VECTOR")
              i+=1
              continue
        print(type(shortest_vector))
        shortest_vector = vector_to_list(shortest_vector)
        if shortest_vector == [0 * dim]:
              print("ERROR: Zero vector")
              i+=1
              continue
        used_vectors.append(shortest_vector)
        # print("basis to be thrown away", basis_as_list[dimension-1])
        basis_as_list[dim-1] = shortest_vector
        # print(basis_as_list)
        # print("basis new, should be SV", basis_as_list[dimension-1])
        basis_as_list.sort(key=lambda x: euclidean_norm(x))
        # print(basis_as_list)
        for i in range(dim):
                print("  b_",i, " ", euclidean_norm(basis_as_list[i]), basis_as_list[i])
        basis_as_matrix = matrix(basis_as_list)
        G = gram_matrix(basis_as_matrix)
    print("reduced lll",  (matrix_to_list(matrix(ZZ, basis_as_list).LLL())), euclidean_norm(first_nonzero_norm_vector(matrix_to_list(matrix(ZZ, basis_as_list).LLL()))))
    return euclidean_norm(basis_as_list[0])



        

# swapping(matrix([[-8, -9, 7], [10, 9, -6], [-8, 4, 10]]))

# def swapping_alg_ver2(dimension, perimeter):
#     """ Verse s for cyklem """
#     basis = randomMatrix(dimension, perimeter)
#     i=0
#     while ještě pořád dělám rozumne změny:
#         perform a cube with fixed index i
#         realmin = vysledek cube tedy něco jkao [1, 0.54, 7.2345]
#         sv = zaokrouhleni realmin * B
#         basis[i] = sv
#         i= i + 1 mod dimension
#         basis.sort by length



# CO KDYŽ VE VĚTŠÍCH DIMENZÍCH SE STREFUJU DO VEKTORU KRATŠÍHO NEŽ LLL???? JE TO POSSIBLE !!!! ZKONTROLOVAT
    
# možná se ještě vyplatí opravit cube takto: co když sii jako real-min uložím něco, co ale je moc daleko a  
# fixací JINÉ souřadnice ale dostanu už správnej výsledek? ještě na to mrkni viko


# 

# In[230]:


def counter(dimension : int, perimeter):
    res = {"True":0, "False":0}
    for i in range(100):
        mojematice = randomMatrix(dimension, perimeter)
        print(mojematice)
        swapresult = swapping(mojematice)
        if swapresult == euclidean_norm(first_nonzero_norm_vector(matrix_to_list(mojematice.LLL()))):
            res["True"]+=1
            continue
        res["False"]+=1
    return res

counter(2, 100)


# In[227]:


def first_nonzero_norm_vector(vectors):
    for v in vectors:
        if euclidean_norm(v) != 0:
            return v
        
def testik(matice):
    matice = matrix(matice)
    lllnorm = euclidean_norm(first_nonzero_norm_vector(matrix_to_list(matrix(ZZ, matice).LLL())))
    return swapping(matice) == lllnorm

# B_1 = [[-8, -9, 7], [10, 9, -6], [-8, 4, 10]]
# assert testik(B_1)

# B_2 = [[9, -1, -6], [7, -9, -5], [9, 6, -2]]
# assert testik(B_2) # NEJEDE

B_3 = [[3, 10, 4], [-4, 2, 7], [-4, -8, -3]]
assert testik(B_3) # FUNGUJE

B_4 =  [[6, -4, -5], [9, 6, -10], [2, 9, 2]]
assert testik(B_4) # FUNGUJE

