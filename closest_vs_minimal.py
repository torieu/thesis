### closest_versus_minimal
from sage.all import *
from copy import deepcopy
import math

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
    norms = []
    for vertice in vertices:
        norms.append([(matrix(vertice) * matrix(B)).norm().n(), vertice])
    # print("     (evaluate norms at cube funct)", sorted(norms))
    return sorted(norms)
    
def shortest_lc_in_cube(lc_cube, B):
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