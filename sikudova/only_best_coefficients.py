import numpy as np


def gram_matrix(matrix):
    """
    Calculate Gram matrix of given matrix.
    :param matrix: given matrix
    :return: Gram matrix of matrix
    """
    return np.matmul(matrix, matrix.transpose())


def calculate_exact_coefficients(indices, given_coefficients, matrix):
    """
    Calculates exact coefficients for matrix and given indexes.
    :param indices: indexes already defined for vectors
    :param given_coefficients: array of exact given coefficients
    :param matrix: matrix of vectors, vectors on indexes already multiplied
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


def change_matrix(indices, matrix):
    """
    Multiply given indices (vectors) with given numbers
    :param indices: array of 2D-arrays
                                - first index: index (row) of matrix
                                - second index: number to multiplication
    :param matrix: matrix of vectors to change
    :return: vector (array) of given coefficients
    """
    all_given_coefficients = [0 for _ in range(len(matrix[0]))]
    for each in indices:
        all_given_coefficients[each[0]] = each[1]
    return all_given_coefficients


def paste_defined_into_exact_coefficients(exact_coefficients, indices, multiply_by, matrix):
    """
    Paste given coefficients for specific vectors into exact coefficients.
    :param exact_coefficients: exact coefficients for vectors except indices
    :param indices: row of vectors with given coefficients
    :param multiply_by: given coefficients for vectors
    :param matrix: matrix of vectors
    :return: array of best coefficients
    """
    best_coefficients = []
    index_a = 0
    index_b = 0
    for i in range(matrix.shape[0]):
        if i in indices:
            best_coefficients.append(multiply_by[index_a])
            index_a += 1
        else:
            best_coefficients.append(exact_coefficients[index_b][0])
            index_b += 1
    return best_coefficients

