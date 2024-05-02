import itertools
from  closest_vs_minimal import *


def descent(jsonfilename):
    """
    Returns: reached: list of number of succesful and unsuccesful trials
    pathways: list of the centers during the alg
    """
    cases = from_json(jsonfilename)
    reached = []
    pathways = []
    i=1
    for case in cases:
        i+=1
        if i==42 or True:
            B, lc_cube, lc_exact = case["B"], case["lincomb_cube"], [int(i) for i in case["lincomb_exact"]]
            r, pathway = descent_single(B, lc_cube, lc_exact)
            reached.append(r)
            pathways.append(pathway)
            if not r:
                plot_pathway_and_nearest_norms(pathway, case["lincomb_cube"], case["lincomb_exact"], case["B"], 4)
        
    return reached, pathways



def descent_single(B, lc_cube, lc_exact):
    pathway = []

    # Append the first, real-valued point of the path. This is used later for plotting.
    pathway.append(lc_cube)
    exact_norm = (matrix(lc_exact) * matrix(B)).norm().n()

    # Find the point of the small cube with smallest norm and set it as a center
    min_norm, center = evaluate_norms_at_cube(cube_points(lc_cube), B)[0]

     # If I reached the lcExact within the first cube, end
    if min_norm == exact_norm:
        pathway.append(center)
        reached = True
        pathway.append(lc_exact)
        return reached, pathway

        
    while True :
        pathway.append(center)
        big_cube = big_cube_points(center)
        bigCubeMinimum_norm, bigCubeMinimum_lc = evaluate_norms_at_cube(big_cube, B)[0]
        # If I reached the lcExact, end with success status
        if bigCubeMinimum_norm == exact_norm:
            pathway.append(bigCubeMinimum_lc)
            reached = True
            break


        # If the new minimum is bigger than the old
        elif bigCubeMinimum_norm > min_norm:

            reached = False
            break

        # Else, continue with the descent
        min_norm, center = bigCubeMinimum_norm, bigCubeMinimum_lc
    return reached, pathway



def different_components(pathway):
    diff_indices = set(range(len(pathway[0])))
    
    for i in range(len(pathway[0])):
        unique_values = {point[i] for point in pathway}
        if len(unique_values) == 1:
            diff_indices.remove(i)
            break

    return list(diff_indices)



def big_cube_points(int_combination):
    points = []
    for i in range((len(int_combination))):
            dupe = int_combination[:]
            dupe[i] += 1
            points.append(copy(dupe))
            dupe[i] -= 2
            points.append(copy(dupe))
    points = [list(map(lambda x : int(x), point)) for point in points]
#     print("big cube points", points)
    return points


def get_index_to_remove(lcCube, lcLLL):
    '''
    Find the index of the identical coordinate shared by lcCube and lcLLL.

    INPUTS:
        lcCube (list): The lcCube point represented as a list of coordinates.
        lcLLL (list): The lcLLL point represented as a list of coordinates.
    '''
    if len(lcCube) != len(lcLLL):
        raise ValueError("Both input lists must have the same length")

    for index, (elem_lcCube, elem_lcLLL) in enumerate(zip(lcCube, lcLLL)):
        if elem_lcCube == elem_lcLLL:
            return index
    return ValueError("The points do not share identical coordinate.")