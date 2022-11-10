import numpy as np
import sympy as sp
import math
import copy

minimal_norm = math.inf
shortest_vector = None
best_coefficients = None


def find_shortest_vector(matrix):
    """
    Find the shortest vector in matrix, its norm and coefficients of linear combination (1 for one vector)
    :param matrix: matrix of vectors
    :return: the norm of the shortest vector, its norm and relevant coefficients
    """
    rows, cols = matrix.shape
    min_vector = None
    min_norm = math.inf
    min_index = math.inf
    for i in range(rows):
        norm = np.linalg.norm(matrix[i])
        if norm < min_norm:
            min_norm = norm
            min_vector = matrix[i]
            min_index = i

    coefficients = []
    for i in range(rows):
        if i == min_index:
            coefficients.append(1)
        else:
            coefficients.append(0)
    return min_vector, min_norm, coefficients


def enumeration(matrix):
    """
    Find the shortest vector in matrix of vectors by recursive help function.
    :param matrix: matrix of vectors in lower triangular shape (done by cholesky decomposition).
    :return: norm of the shortest vector, the shortest vector and relevant coefficients
    """
    rows, cols = matrix.shape

    global minimal_norm
    global shortest_vector
    global best_coefficients

    best_coefficients = []
    shortest_vector = []

    vector, minimal_norm, coefficients = find_shortest_vector(matrix)

    best_coefficients.append(coefficients)
    shortest_vector.append(vector)

    # upper and lower boundary, then round to integer
    from_to = abs(minimal_norm) / abs(matrix[rows - 1][rows - 1])
    do_from = round(-1 * from_to)
    do_to = round(from_to)

    matrix_copy = copy.deepcopy(matrix)
    # call recursion function on the last vector of triangular matrix
    enum_help(matrix_copy, do_from, do_to, rows - 1, None, [None for _ in range(rows)])

    # check correctness of the result
    for vector in shortest_vector:
        assert (minimal_norm == np.linalg.norm(vector))
    for each in best_coefficients:
        check_result = np.zeros(rows)
        for i in range(rows):
            for j in range(rows):
                check_result[i] += each[j] * matrix[j][i]
        is_ok = False
        for vector in shortest_vector:
            if np.allclose(check_result, vector):
                is_ok = True
            if np.allclose((-1) * check_result, vector):
                is_ok = True
        if not is_ok:
            print("ERROR!!!")
            exit(1)

    return minimal_norm, shortest_vector, best_coefficients


def enum_help(matrix, do_from, do_to, row, vector_before, coefficients):
    global minimal_norm
    global shortest_vector
    global best_coefficients

    # the first vector, end of recursion
    if row == 0:
        for i in range(do_from, do_to + 1):
            # current vector is current vector multiplied by i
            current_vector = i * matrix[row]

            # add i as a coefficient of linear combination
            coefficients[row] = i

            # vector before plus current vector, vector before is a linear combination of all vectors before so far
            vector = copy.deepcopy(vector_before + current_vector)

            # norm of vector part we can not change
            temp_sum = np.linalg.norm(vector[row:])

            # vector of the same minimal norm found
            if temp_sum == minimal_norm and temp_sum != 0:
                best_coefficients.append(copy.deepcopy(coefficients))
                shortest_vector.append(copy.deepcopy(vector))
            # shorter vector found
            if temp_sum < minimal_norm and temp_sum != 0:
                minimal_norm = temp_sum
                shortest_vector = [copy.deepcopy(vector)]
                best_coefficients = [copy.deepcopy(coefficients)]
    # other vectors
    else:
        for i in range(do_from, do_to + 1):
            # current vector is current vector multiplied by i
            current_vector = i * matrix[row]

            # add i as a coefficient of linear combination
            coefficients[row] = i

            # bigger norm than minimal norm
            if abs(current_vector[-1]) >= minimal_norm:
                continue

            # vector is the current vector plus vector before
            if row + 1 == matrix.shape[0]:
                vector = copy.deepcopy(current_vector)
            else:
                vector = copy.deepcopy(vector_before + current_vector)

            # norm of vector part we can not change
            temp_sum = np.linalg.norm(vector[row:])

            # if temp sum is less than minimal_norm
            if temp_sum < minimal_norm:
                # norm we can add to norm of vector
                remain_norm = math.sqrt(minimal_norm ** 2 - temp_sum ** 2)

                # solve equation to gain lower and upper boundary
                x = sp.Symbol('x', real=True)
                result = sp.solve(
                    sp.Abs(vector[row - 1] + matrix[row - 1][row - 1] * x) - remain_norm, x)
                if not result:
                    continue

                # lower bound and upper bound to integer from real number
                lower_bound, upper_bound = result[0], result[1]
                if lower_bound > 0:
                    lower_bound = math.ceil(lower_bound)
                else:
                    lower_bound = math.ceil(lower_bound)
                if upper_bound > 0:
                    upper_bound = math.floor(upper_bound)
                else:
                    upper_bound = math.floor(upper_bound)

                # recursion call on upper vector
                enum_help(matrix, lower_bound, upper_bound, row - 1, copy.deepcopy(vector), coefficients)
