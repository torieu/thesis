import json
from random import randint, seed
from sage.all import *

# Utilitary math tools

def shortestVector(matrix):
    """Returns the shortest vector of the given matrix."""
    return sorted([(matrix[idx].norm().n(), matrix[idx], idx) for idx in range(matrix.nrows())])[0][1]

def randomMatrix(dimension, per):
    """Returns a random square matrix with full rank."""
    list = [randint(-per, per) for _ in range(dimension**2)]
    M = matrix(ZZ, dimension, dimension, list)
    while M.rank() != dimension:
        list = [randint(-per, per) for _ in range(dimension**2)]
        M = matrix(ZZ, dimension, dimension, list)
    return M

def gram_matrix(matrix):
    """Returns the Gram matrix of the given matrix."""
    return matrix * matrix.transpose()

def matrix_to_list(A):
    """Converts a matrix into a list of lists."""
    return [[int(A.list()[row * A.ncols() + col]) for col in range(A.ncols())] for row in range(A.nrows())]

def vector_to_list(vector):
    """Converts a vector into a list."""
    return [float(num) for num in vector.list()]

# Formatting tools

def format_data(output_data):
    """Formats the data for JSON output."""
    for dic in output_data:
        for key, value in dic.items():
            if isinstance(value, sage.matrix.matrix_integer_dense.Matrix_integer_dense):
                dic[key] = matrix_to_list(value)
            elif isinstance(value, float):
                dic[key] = float(value)
            else:
                dic[key] = vector_to_list(value)
    return output_data

def into_json(output_data, jsonfilename):
    """Saves data to a JSON file."""
    with open(jsonfilename, "w+") as out_file:
        json.dump(output_data, out_file)

def from_json(filename):
    """Loads data from a JSON file."""
    with open(filename) as out_file:
        return json.load(out_file)

def print_listdict(list):
    """Prints a list of dictionaries."""
    for dictionary in list:
        for key, value in dictionary.items():
            print(f"{key}: {value}")
        print()

def into_dict(B, lcLLL, rowindex: int, dimension) -> dict:
    """Converts data into a dictionary format."""
    result = {}
    G = gram_matrix(B)
    lcCube = find_real_minimum(G, rowindex, lcLLL[rowindex], dimension)
    result["B"] = B
    result["G"] = G
    result["lincomb_LLL"] = lcLLL
    result["lincomb_cube"] = lcCube
    result["sv_LLL"] = lcLLL * B
    result["sv_cube"] = lcCube * B
    result["lincomb_diff"] = lcLLL - vector(lcCube)
    result["LLL.norm"] = float(vector(lcLLL*B).norm().n(digits=5))
    result["cube.norm"] = float((vector(lcCube)*B).norm().n(digits=5))
    return result

# Main functions

def generate_new_examples(iterations, dimension, perimeter, sensitivity, jsonfilename, printing=False, functioning=True):
    """
    Generates new examples and checks if a matrix is dysfunctional.
    Output: number of dys matrices.
    """
    output_data = []
    for _ in range(iterations):
        B = randomMatrix(dimension, perimeter)
        v_min = shortestVector(B.LLL())
        lcLLL = B.solve_left(v_min)
        G = gram_matrix(B)
        booltmp, rowindex = is_dysfunctional(B, v_min, lcLLL, G, dimension, sensitivity)
        if booltmp != functioning:
            case_info = into_dict(B, lcLLL, rowindex, dimension)
            output_data.append(case_info)
    if printing:
        print_listdict(output_data)
    into_json(format_data(output_data), jsonfilename)
    return len(output_data)

def is_dysfunctional(B, v_min, lcLLL, G, dimension, sensitivity):
    """Checks if the matrix is dysfunctional."""
    nonzero_ind = 0
    for current_row in range(dimension):
        lcCube = find_real_minimum(G, current_row, lcLLL[current_row], dimension)
        if lcCube == zero_vector(SR, dimension):
            continue
        nonzero_ind = current_row
        for i in range(dimension - 1):
            difference = lcLLL[i] - lcCube[i]
            if abs(difference) >= sensitivity:
                return True, current_row
    return False, nonzero_ind

def find_real_minimum(G, current_row, lcLLL, dimension):
    """Finds the real minimum of the given matrix."""
    matrixA = matrix(dimension - 1, dimension - 1, 0)
    matrixB = matrix(dimension - 1, 1, 0)
    matrixA[0, 0] = 1
    a, b = 0, 0
    for row in range(dimension):
        if row != current_row:
            matrixA[a] = [G[row, j] for j in range(len(G[row])) if j != current_row]
            matrixB[b] = sum([lcLLL * G[row, j] for j in range(len(G[row])) if j == current_row])
            a += 1
            b += 1
    result = (matrixA.solve_right((-1) * matrixB)).list()
    result.insert(current_row, lcLLL)
    return vector(result).n(digits=5)
