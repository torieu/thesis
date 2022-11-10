import copy
import numpy as np
from copy import deepcopy
from math import floor, ceil

from cholesky_decomposition import cholesky_decomposition
from zero_one_combination import zero_one_combination

from integer_files import lattices_10_1


def minimum_norm(matrix):
    """
    Return norm of the shortest vector in the matrix.
    :param matrix: matrix where to find the norm
    :return: the minimal norm
    """
    rows, cols = matrix.shape
    magnitudes = []

    for i in range(rows):
        magnitudes.append(np.linalg.norm(matrix[i]))

    return min(magnitudes)


def calculate_exact_coefficients(indices, given_coefficients, matrix):
    """
    Calculates exact coefficients for matrix and given indexes.
    :param indices: indexes already defined for vectors
    :param given_coefficients: array of exact given coefficients
    :param matrix: matrix of vectors
    :return: array of exact coefficients
    """
    dot_matrix = gram_matrix(matrix)
    rows, cols = dot_matrix.shape
    matrix_a = np.empty([rows - len(indices), rows - len(indices)])
    matrix_b = np.empty([rows - len(indices), 1])
    a_index, b_index = 0, 0
    for row in range(rows):
        if row not in indices:
            matrix_a[a_index] = [dot_matrix[row][x] for x in range(len(dot_matrix[row])) if x not in indices]
            matrix_b[b_index] = sum(
                [given_coefficients[x] * dot_matrix[row][x] for x in range(len(dot_matrix[row])) if x in indices])
            a_index += 1
            b_index += 1
    result = np.linalg.solve(matrix_a, -1 * matrix_b)
    return result


def gram_matrix(matrix):
    """
    Calculate Gram matrix of given matrix.
    :param matrix: given matrix
    :return: Gram matrix of matrix
    """
    return np.matmul(matrix, matrix.transpose())


def reduce_vector_binary(matrix: np.ndarray, do_from: int, do_to: int) -> object:
    """
    Find better linear combination of vectors (except the first one) to the first vector in matrix.
    :param matrix: matrix of vectors
    :param do_from: lower boundary of interval
    :param do_to: upper boundary of interval
    :return: array of vector with the minimal norm, coefficients and the minimal norm
    """
    rows, cols = matrix.shape

    best_coefficients = []
    best_vectors = [copy.deepcopy(matrix[0])]

    first_vector = copy.deepcopy(matrix[0])
    min_magnitude = np.linalg.norm(matrix[0])

    # generate numbers to create good linear combinations
    zero_one_array = zero_one_combination(rows - 1)

    for k in range(do_from - 1, do_to + 2):
        # multiply first vector with all possible coefficients (THE IMPROVEMENT)
        # temp_first_vector = k * first_vector
        # matrix[0] = copy.deepcopy(temp_first_vector)

        # calculate gram matrix with multiplied first vector
        dot_matrix = gram_matrix(matrix)

        # coefficients, vector of real numbers
        array = [0 for _ in range(rows)]
        array[0] = k
        coefficients = calculate_exact_coefficients([0], array, matrix)

        # coefficients, vector of integers, from coefficients items - floored and ceiled
        rounded_coefficients = [[] for _ in range(len(coefficients))]
        for i in range(len(coefficients)):
            rounded_coefficients[i] = [floor(coefficients[i]), floor(coefficients[i] + 1)]

        starting_coefficients = [x[0] for x in rounded_coefficients]
        starting_vector = copy.deepcopy(matrix[0] * k)

        for i in range(rows - 1):
            starting_vector += starting_coefficients[i] * matrix[i + 1]

        temp_vector = copy.deepcopy(starting_vector)
        vector = [0 for _ in range(rows)]
        # all coefficients generated with the help of good sorted zero_one array
        for each in zero_one_array:
            coeff = copy.deepcopy(starting_coefficients)
            if each < 0:
                temp_vector -= matrix[abs(each) + 1]
                vector[abs(each) + 1] = 0
            else:
                temp_vector += matrix[abs(each) + 1]
                vector[each + 1] = 1

            if np.linalg.norm(temp_vector) < min_magnitude and not abs(np.linalg.norm(temp_vector - 0)) <= 1e-5:
                best_vectors = [copy.deepcopy(temp_vector)]
                min_magnitude = np.linalg.norm(temp_vector)
                coeff.insert(0, k)
                best_coefficients = [np.add(vector, coeff)]
            elif np.linalg.norm(temp_vector) == min_magnitude:
                if not any((temp_vector == x).all() for x in best_vectors):
                    best_vectors.append(copy.deepcopy(temp_vector))
                coeff.insert(0, k)
                if not any((np.add(vector, coeff) == x).all() for x in best_coefficients):
                    best_coefficients.append(np.add(vector, coeff))

    matrix[0] = copy.deepcopy(first_vector)

    # assert to check the norm
    for i in range(len(best_vectors)):
        assert (min_magnitude == np.linalg.norm(best_vectors[i]))

    return best_vectors, best_coefficients, min_magnitude


def main_reduce_function(matrix):
    """
    Iterates over all vectors and compares norms of vectors from reduce_function.
    :param matrix: matrix of vector
    :return: array of shortest found vectors, relevant coefficients and the minimal norm
    """

    rows, cols = matrix.shape

    # variables for resulting shortest vector, its norm and its coefficients
    result_vector = []
    best_coefficients = []
    min_magnitude = minimum_norm(matrix)

    matrix_copy = copy.deepcopy(matrix)

    for i in range(rows):
        if i != 0:
            temp = deepcopy(matrix[0])
            matrix[0] = matrix[i]
            matrix[i] = temp

        temp = deepcopy(matrix_copy[rows - 1])
        matrix_copy[rows - 1] = matrix_copy[i]
        matrix_copy[i] = temp

        # decompose matrix into lower_triangular matrix
        lower_triangular = cholesky_decomposition(gram_matrix(matrix_copy))

        # compute lower and upper boundaries
        from_to = min_magnitude / lower_triangular[rows - 1][rows - 1]
        do_from = round(-1 * from_to)
        do_to = round(from_to)

        # find the best linear combination for i-th vector with interval
        best_vector, coefficients, magnitude = reduce_vector_binary(matrix, do_from, do_to)

        for j in range(len(coefficients)):
            coefficients[j][0:i + 1] = list(np.roll(coefficients[j][0:i + 1], i))

        if magnitude < min_magnitude:
            result_vector = copy.deepcopy(best_vector)
            best_coefficients = copy.deepcopy(coefficients)
            min_magnitude = magnitude
        elif magnitude == min_magnitude:
            for each in best_vector:
                result_vector.append(deepcopy(each))
            for each in coefficients:
                best_coefficients.append(copy.deepcopy(each))

    matrix = np.roll(matrix, -1, axis=0)

    for each in best_coefficients:
        check_result = np.zeros(cols)
        for i in range(rows):
            check_result += each[i] * matrix[i]
        is_ok = False
        for vector in result_vector:
            if np.allclose(check_result, vector):
                is_ok = True
            if np.allclose((-1) * check_result, vector):
                is_ok = True
        if not is_ok:
            print("ERROR!!!")
            exit(1)

    return min_magnitude, result_vector, best_coefficients
