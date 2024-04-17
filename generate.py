# UNINTERESTING TOOLS
from random import randint, seed
from sage.all import *


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


def vector_to_int_list(vector) -> list:
    return [int(num) for num in vector.list()]

def euclidean_norm(vector):
    return math.sqrt(sum(x**2 for x in vector))




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


def into_dict(B, lcLLL: vector, lcCube: int, lc_exact) -> dict:
    result = {}
    G = gram_matrix(B)
    lc_shortest_in_cube=  vector(shortest_lc_in_cube(lcCube, B))
    result["B"] = B
    result["G"] = G
    result["lincomb_exact"] = lc_exact
    result["lincomb_LLL"] = lcLLL
    result["lincomb_cube"] = lcCube
    result["lincomb_cube_shortest"] = lc_shortest_in_cube
    result["exact.norm"] = vector(lc_exact*B).norm().n(digits=5)
    result["LLL.norm"] = vector(lcLLL*B).norm().n(digits=5)
    result["shortest_cube.norm"] = (vector(lc_shortest_in_cube)*B).norm().n(digits=5)
    result["sv_LLL"] = lcLLL * B
    result["sv_exact"] = lc_exact * B
    result["sv_shortst_cube"] = lc_shortest_in_cube* B
    return result


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
        
def find_real_minimum(G, current_row, inserted_component) -> vector:
    dimension = len(G)
    print(type(dimension))
    matrixA = matrix(dimension - 1, dimension - 1, 0) # square matrix of size (dimension - 1) x (dimension - 1), filled with zeros.
    matrixB = matrix(dimension - 1, 1, 0) # column matrix of size (dimension - 1) x 1, filled with zeros.
    matrixA[0,0] = 1
    a, b = 0, 0
    for row in range(dimension):
        if row != current_row:
            matrixA[a] = [G[row, j] for j in range(len(G[row])) if j != current_row]
            matrixB[b] = sum([inserted_component * G[row,j] for j in range(len(G[row])) if j == current_row])
            a += 1
            b += 1
    # insert indices
    print(A,"A")
    result = (matrixA.solve_right((dimension) * matrixB)).list()
    result.insert(current_row, inserted_component)
    return vector(result).n(digits=5)

# generate_new_examples(1, 3, 10, 3, "dummyfile.json", printing = False, functioning = False)
