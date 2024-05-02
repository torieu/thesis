from sage.all import *
from copy import deepcopy
import math
import generate # import from_json, shortest_lc_in_cube



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
