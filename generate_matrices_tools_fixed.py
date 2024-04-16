#!/usr/bin/env python
# coding: utf-8

# # This one fixes one component to the same as LLL
# 

# In[2]:


# UNINTERESTING TOOLS
from random import randint, seed

def shortestVector(matrix): 
    """
    :returns: 
        - n - norm of the shortest vector
        - v - the shortest vector (vector)
        - i - index of the SV
    """
    return sorted([(matrix[idx].norm().n(), matrix[idx], idx) for idx in range(matrix.nrows())])[0][1]

def randomMatrix(dimension, per) -> matrix:
    '''
    Returns a random square matrix with full rank.
    INPUT:
    dimension: dimension of the random matrix
    per: perimeter of the components
    '''
    list = [randint(-per, per) for _ in range(dimension**2)]
    M = matrix(ZZ, dimension, dimension, list)
    while M.rank() != dimension:
        list = [randint(-per, per) for _ in range(dimension**2)]
        M = matrix(ZZ, dimension, dimension, list)
    return M

def gram_matrix(matrix) -> matrix:
    return matrix * matrix.transpose()

def matrix_to_list(A) -> list:
    return [[int(A.list()[row * A.ncols() + col]) for col in range(A.ncols())] for row in range(A.nrows())]

def vector_to_list(vector) -> list:
    return [float(num) for num in vector.list()]


# ### Interesting math tools. 
# 
# Following functions actually have some thoughts behind them. 
# 
# 
# `abnormality_data`
# Given a matrix $B$ and it's shortest vector given by the LLL algorithm, returns __a dictionary__ with data concerting the SV given by the cube algorithm, if it's far enough from the LLL solution. 
# 
# `find_real_minimum` Given a matrix $B$ and it's shortest vector given by the LLL algorithm, returns the __linear combination of the SV given by the cube algorithm__.
# 
# 
# 

# ### Formatting tools
# 
# These functions provide transfering the data to the output jsonfile as described above.

# In[3]:


### output formatting
import json
import numpy as np

def format_data(output_data):
    for dic in output_data:
        for key, value in dic.items():
            if isinstance(value, sage.matrix.matrix_integer_dense.Matrix_integer_dense):
                dic[key] = matrix_to_list(value)
            elif isinstance(value, (int, float, sage.rings.real_mpfr.RealNumber)):
                dic[key] = float(value) 
            else:
                dic[key] = vector_to_list(value)
    return output_data


def into_json(data, jsonfilename):
    with open(jsonfilename, 'w') as f:
        json.dump(data, f)

def from_json(filename):
    out_file = open(filename)
    return json.load(out_file)

def print_listdict(list) -> None:
    """
    :param list: list of dictionaries
    """
    for dictionary in list:
        for pair in dictionary.items():
            print(pair[1], ": ", pair[0])
        print()

def norm_G(lin_comb, Gram_matrix):
    tmp = [sum(x * y for x, y in zip(lin_comb, col)) for col in zip(*Gram_matrix)]
    return sum(x * y for x, y in zip(tmp, lin_comb))
    
def into_dict(B, lcLLL: vector, lcCube: int) -> dict:
    result = {}
    G = gram_matrix(B)
    result["B"] = B
    result["G"] = G
    result["lincomb_LLL"] = lcLLL
    result["lincomb_cube"] = lcCube
    result["LLL.norm"] = vector(lcLLL*B).norm().n(digits=5)
    result["cube.norm"] = (vector(lcCube)*B).norm().n(digits=5)
    result["sv_LLL"] = lcLLL * B
    result["sv_cube"] = lcCube * B
    return result


# ### Main function
# 
# This cell actually generates the new examples and checks wheter a matrix is *dysfunctional*. 
# 
# `generate_new_examples`
# The main function. Based on input parameters above, generates *some* number of *dysfunctional* matrices, computes their invariants such as Gram matrix, LLL/cube linear combinations etc. and creates a json file with this information (specified above).
# 
# `is_dysfunctional`
# Given a matrix $B$, checks whether its dysfunctional and if so, returns the data describing the case.
# 

# In[4]:


# MAIN 
def generate_new_examples(iterations, dimension, perimeter, rng, jsonfilename, printing = False, functioning = True) -> None:
    '''
    generates >iterations< of >dis/functioning< matrices and saves them in a json file.
    returns number of suitable cases
    '''
    output_data = []
    it = 0
    for i in range(iterations):
    # while it < iterations:
        B = randomMatrix(dimension, perimeter)
        G = gram_matrix(B)

        # Solve the SVP exactly
        # FIXME
        v_min = shortestVector(B.LLL())
        lcLLL = B.solve_left(v_min)

        # Compare the result with solution of the cube algorithm
        this_case_works, lcCube = compare_with_cube(lcLLL, G, rng)
        if this_case_works != functioning:
            case_info = into_dict(B, lcLLL, lcCube)
            output_data.append(case_info)
            it +=1
    
    # Format output data
    if printing: print_listdict(output_data)
    dict = format_data(output_data)
    into_json(dict, jsonfilename)
    return len(output_data)


def compare_with_cube(lcLLL, G, rng) -> bool:
    working_case = None
    for current_row in range(dimension):
        # for K in range(2, rng):
        LLL_component = lcLLL[current_row]
        lcCube = find_real_minimum(G, current_row, LLL_component)
        working_case = lcCube                

        # Check for invalid cases
        if norm_G(lcCube, G) > norm_G(lcLLL, G) or lcCube == zero_vector(SR, dimension):
            continue

        # Check if the two results are the same
        for i in range(dimension - 1):  
            difference = abs(lcLLL[i]) - abs(lcCube[i])
            if abs(difference) >= sensitivity:
                return True, lcCube
    return False, working_case
        
def find_real_minimum(G, current_row, LLL_component) -> vector:
    """
    returns linear combination of the minimum over real numbers
    the >current_row< parameter is the coordinate of fixed col/row idk
    the >K< parameter is the number inserted into this column
    """
    matrixA = matrix(dimension - 1, dimension - 1, 0) # square matrix of size (dimension - 1) x (dimension - 1), filled with zeros.
    matrixB = matrix(dimension - 1, 1, 0) # column matrix of size (dimension - 1) x 1, filled with zeros.
    matrixA[0,0] = 1
    a, b = 0, 0
    for row in range(dimension):
        if row != current_row:
            matrixA[a] = [G[row, j] for j in range(len(G[row])) if j != current_row]
            matrixB[b] = sum([LLL_component * G[row,j] for j in range(len(G[row])) if j == current_row])
            a += 1
            b += 1
    # insert indices
    result = (matrixA.solve_right((-1) * matrixB)).list()
    result.insert(current_row, LLL_component)
    return vector(result).n(digits=5)


# In[18]:


dimension = 6
perimeter = 1000
sensitivity = 1
tries = 10000
jsonfilename = "basic_cube_nonf_matrices_%sx%s_%s_tries_%s_perimeter.json"% (dimension, dimension, tries ,perimeter)
range_for_K = 4
generate_new_examples(tries, dimension, perimeter, range_for_K, jsonfilename, printing=False, functioning=False)


# In[19]:


dimension = 6
perimeter = 10000
sensitivity = 1
tries = 10000
jsonfilename = "basic_cube_nonf_matrices_%sx%s_%s_tries_%s_perimeter.json"% (dimension, dimension, tries ,perimeter)
range_for_K = 4
generate_new_examples(tries, dimension, perimeter, range_for_K, jsonfilename, printing=False, functioning=False)

