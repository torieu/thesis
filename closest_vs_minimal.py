### closest_versus_minimal
from sage.all import *
from copy import deepcopy
import math
from generate import *

def closest_versus_grammatrix(jsonfilename):
    cases = from_json(jsonfilename)
    for case in cases:
        closest = closest_point_in_cube(case["lincomb_cube"], case["lincomb_LLL"])
        point_versus_grammatrix(case["G"], case["lincomb_LLL"], closest)
    


    
def point_versus_grammatrix(G: list, lc_LLL: list, point: list) -> None:
    detailed_product, sum = matrix_multiplication_detailed(matrix(G), vector(point)-vector(lc_LLL))
    print("lc in question: {}, lc LLL: {}". format(point, lc_LLL))
    print("difference between the two lcs/points:", vector(point)-vector(lc_LLL))
    print(detailed_product)
    print(sum)
    print()
        

def cube_versus_grammatrix(jsonfilename):
    cases = from_json(jsonfilename)
    for i,case in enumerate(cases):
        print("\n \033[1m ------ NEXT CASE --------------\033[0m\n")
        print("the gram matrix:\n", matrix(case["G"]), "\n")
        for point in cube_points(case["lincomb_cube"]):
            if point == closest_point_in_cube(case["lincomb_cube"], case["lincomb_LLL"]):
                print('\033[1m' + "THIS IS THE CLOSEST POINT IN THIS CASE" + '\033[0m')
            if point == evaluate_norms_at_cube(cube_points(case["lincomb_cube"]), case["B"])[0][1]:
                print('\033[1m' + "THIS IS THE MINIMAL POINT IN THIS CASE" + '\033[0m')    
            point_versus_grammatrix(case["G"], case["lincomb_LLL"], point)
        if i > 10:
            break


### closest_versus_minimal

def wholify_floatInListAtGivenIndex(inlist, i):
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

def evaluate_norms_at_cube_efficient(floored_lc, B, fixed_index):
    """
    B: matrix as a list
    floor: lc_cube but every float is rounded down
    """
    dim = len(B)
    B = matrix(B)
    # floored_times_B = rozepsane_nasobeni(floored_lc, B)
    floored_lc = vector(floored_lc)
    floored_vector = floored_lc * B
    norms_and_vertices = [floored_vector]

    for i in range(1, 2**(dim-1)):
        bin_str = bin(i)[2:].zfill(dim-1)
        tuple_result = [int(char) for char in bin_str]
        tuple_result.insert(fixed_index, 0)
        print(tuple_result,"tuple")
        indices = [idx for idx, char in enumerate(bin_str) if char == '1']
        this_vector=deepcopy(floored_vector)
        print(this_vector, "this vector")
        for indice in indices:
            print(B[indice])
            this_vector += B[indice]
        print(this_vector, "this vector")
        norms_and_vertices.append(vector(this_vector))
        print()
    return norms_and_vertices




def rozepsane_nasobeni(v, A):
    """
    tvaru row vector * matice, ale vysledkem neni row, ale matice
    kde se musi sečist čisilka ve sloupci a dostanu tak rows
    """
    A = matrix(A)
    rows = A.nrows()
    cols = A.ncols()
    newA = matrix(rows, cols, 0)
    for i in range(rows):
        for j in range(cols):
            newA[i,j] = A[i,j] * v[i]
    return newA

def secist_rozepsane_nasobeni_do_vektoru(A):
    """
    input: sage matrix A
    output: sage vector 
    absolutlely terrible implementation, is it actually helpful?
    """
    cols = [0] * A.ncols()
    for i in range(A.ncols()):
        for j in range(A.nrows()):
            cols[i] += A[i,j]
    return vector(cols)



B = [
    [8,-7, -5],
    [-1, 1, 2],
    [6, -4, -7]
]

floored_lc = [-1, 2,1]



cube = [-0.3, 2.9, 1]
print("efficient returned: ",evaluate_norms_at_cube_efficient(floored_lc, B,2))

for i in range(4):
    print("  ", vector(evaluate_norms_at_cube(cube_points(cube), B)[i][1])*matrix(B))

def shortest_lc_in_cube(lc_cube, B) -> list:
    return evaluate_norms_at_cube(cube_points(lc_cube), B)[0][1]


def closest_point_in_cube(lc_cube, lc_LLL) -> list:
    """
    Returns the vertice of the cube closest to the lincomb given by LLL.
    """
    lc_closest = []
    for i in range(len(lc_cube)):
        if lc_LLL[i] > lc_cube[i]:
            lc_closest.append(math.ceil(lc_cube[i]))
        else:
            lc_closest.append(math.floor(lc_cube[i]))
    return lc_closest




def closest_versus_minimal(jsonfilename) -> list:
    "Returns a list of a form [int, int]. The first number indicates  whether the closest vector was the one with shortest norm."
    cases = from_json(jsonfilename)
    stats = [0, 0]
    for case in cases:
        B, lc_cube, lc_LLL = case["B"], case["lincomb_cube"], case["lincomb_LLL"]
        minimal = shortest_vector_in_cube(lc_cube, B)
        closest = closest_point_in_cube(lc_cube, lc_LLL)
        if minimal == closest:
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