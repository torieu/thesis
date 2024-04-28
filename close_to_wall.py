from generate import *

def detect_rounding(jsonfilename, limits):
    '''
    OUTPUT: dict in the form \{ limit: [succesful attempts, unsuc att.],...\}
    [# of succ rounding, # of unsucc rounding, #unsuc in correct direction, # unsucc in wrong direction]
    '''
    cases = from_json(jsonfilename)
    stats = {}
    for limit in limits:
        stats[limit] = [0] * 4
    for case in cases:
        for i in range(len(case["lincomb_cube"])):

            num_cube = case["lincomb_cube"][i]
            num_exact = case["lincomb_exact"][i]
            diff = abs(num_cube) - abs(num_exact)

            #skip equal components
            if diff == 0:
                continue
        
            for limit in limits:
                is_near_int, direction =  is_nearest_integer(num_cube, limit)

                if is_near_int: # if its almost an int
                    if round(num_cube) == num_exact: 
                        stats[limit] = increment(stats, limit, 0)                        # rounding lead directly to lcLLL
                    else:
                        stats[limit] = increment(stats, limit, 1)                        # rounding didnt lead directly to lcLLL
                        if lc_exact_is_in_same_direction(direction, num_cube, num_exact):
                            stats[limit] = increment(stats, limit, 2)                    # but it was in the same direction
                        else: stats[limit] = increment(stats, limit, 3)                  # not even the same direction
    return stats

def is_nearest_integer(num, tolerance):
    nearest_int = round(num)
    is_within_tolerance = abs(num - nearest_int) <= tolerance
    if not is_within_tolerance:
        return False, None
    direction = "up" if num - nearest_int < 0 else "down"
    return True, direction


def increment(dic, limit, ind):
    lst = [dic.get(limit,0)[0], dic.get(limit,0)[1],  dic.get(limit,0)[2],  dic.get(limit,0)[3]]
    lst[ind] += 1
    return lst

def lc_exact_is_in_same_direction(direction, num_cube, num_exact):
    return abs(num_cube - num_exact) - abs(round(num_cube) - num_exact) > 0
