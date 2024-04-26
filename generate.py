# UNINTERESTING TOOLS
from random import randint, normalvariate
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




### formatting tools
import json

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



def random_special_matrix(dim, rng):
    """
    generates a random matrix of dimension >dim< in the form
    (d_1    1   0   ...)
    (d_1    0   1   ...)
    (...    ... ... ...)
    (d_dim    0   ... 0  )

    where ds are relatively large integers in the range >rng<
    """
    M = matrix(ZZ, dim) 
    while M.rank() != dim:
        for i in range(dim):
            M[i, 0] = randint(-rng, rng)  
            if i > 0:
                M[i-1, i] = 1
    return M
