# UNINTERESTING TOOLS
from random import randint 
from sage.all import *
from copy import deepcopy
import json
import math


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





############ CLOSEST VS MINIMAL




### closest_versus_minimal

def wholify_floatInListAtGivenIndex(inlist, i) -> list:
    """
    Given a list, with a float on the i-th position, 
    """
    outlist = deepcopy(inlist)
    for _ in range(len(outlist)):
        vector = outlist.pop(0)
        vector[i] = math.floor(vector[i])
        outlist.append(deepcopy(vector))
        vector[i] = vector[i] + 1
        outlist.append(vector)
    return outlist




def cube_points(lc_cube) -> list:
    """
    Returns a list of linear combinations/vertices of the cube.
    """
    points = [lc_cube]
    for i in range((len(lc_cube))):
        if points[0][i] % 1 != 0: # if current number is not whole 
            points = wholify_floatInListAtGivenIndex(points, i) 
    points = [list(map(lambda x : int(x), point)) for point in points]
    return points



    
def evaluate_norms_at_cube(vertices, B):
    """
    Given a list of vertices and an integer cube and returns the norms of its vertices and corresponding vertices. 
    """
    vertices = [vertice for vertice in vertices if not all(v == 0 for v in vertice)]
    norms_and_vertices = []
    for vertice in vertices:
        norms_and_vertices.append([(matrix(vertice) * matrix(B)).norm().n(), vertice])
    return sorted(norms_and_vertices)


def closest_point_in_cube(lc_cube, lc_exact) -> list:
    """
    Returns the vertice of the cube closest to the lincomb given by LLL.
    """
    vertices = cube_points(lc_cube)
    lc_closest = None
    min_distance = 10000
    for vertice in vertices:
        this_distance = 0
        for i in range(len(vertice)):
            this_distance += abs(lc_exact[i]-vertice[i])
        if this_distance < min_distance:
            lc_closest = vertice
            min_distance= this_distance
    return lc_closest




def closest_versus_minimal(jsonfilename) -> list:
    "Returns a list of a form [int, int]. The first number indicates  whether the closest vector was the one with shortest norm."
    cases = from_json(jsonfilename)
    stats = [0, 0]
    for case in cases:
        B, lc_cube, lc_exact = case["B"], case["lincomb_cube"], case["lincomb_exact"]
        minimal = shortest_lc_in_cube(lc_cube, B)
        closest = closest_point_in_cube(lc_cube, lc_exact)
        minus_minimal = [-1* i for i in minimal]
        if minimal == closest or minus_minimal==closest:
            stats[0] += 1
        else:
            stats[1] += 1
    return stats




def closest_versus_minimal_withComments(jsonfilename) -> list:
    "Returns a list of a form [int, int]. The first number indicates  whether the closest vector was the one with shortest norm."
    cases = from_json(jsonfilename)
    stats = [0, 0]
    for case in cases:
        B, lc_cube, lc_LLL = case["B"], case["lincomb_cube"], list(map(lambda x : int(x), case["lincomb_LLL"] ))
        print("lc_LLL: ", list(map(lambda x : int(x), lc_LLL )))
        minimal = evaluate_norms_at_cube(cube_points(lc_cube), B)[0][1] # just a linear combination!
        print("minimal:", minimal)
        closest = closest_point_in_cube(lc_cube, lc_LLL)
        print("closest:", closest)
        print("LLL-min:", [lc_LLL[i] - minimal[i] for i in range(len(minimal))])
        print("LLL-clo:", [lc_LLL[i] - closest[i] for i in range(len(minimal))])
        print("Does minimal equal the closest?", minimal == closest)
        # print(cube_points(lc_cube))
        print()
        if minimal == closest:
            stats[0] += 1
        else:
            stats[1] += 1
    return stats


def shortest_lc_in_cube(lc_cube, B) -> list:
    return evaluate_norms_at_cube(cube_points(lc_cube), B)[0][1]
