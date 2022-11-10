import math
from copy import deepcopy
import numpy as np
import cholesky_decomposition


def generate_all_possible_options(original_vector, temp_vector, index, count_positives, temp_sum, all_options):
    """
    Generates all possible combinations of vectors with given first vector.
    Combination depends only on the first vector sum.
    :param original_vector: sorted (decreasingly) first row of lattice
    :param temp_vector: coefficients for linear combination
    :param index: position (where are we right now?)
    :param count_positives: number of positives elements in the first = original vector
    :param temp_sum: temporary sum, moving to the right side
    :param all_options: consecutively appending possible = acceptable temp vectors
    :return: nothing, change all_options array
    """

    if index == len(original_vector) - 1:
        if temp_sum + 2 * original_vector[index] < 0:
            temp_vector[index] = 1
            all_options.append(temp_vector)
        if temp_sum < 0:
            temp_vector[index] = None
            all_options.append(temp_vector)
        return

    if index >= count_positives:
        negatives_sum_from = index + 1
    else:
        negatives_sum_from = count_positives

    if temp_sum + 2 * original_vector[index] + 2 * sum(original_vector[negatives_sum_from:]) < 0:
        temp_vector[index] = 1
        generate_all_possible_options(original_vector, deepcopy(temp_vector), index + 1, count_positives,
                                      temp_sum + 2 * original_vector[index], all_options)
        temp_vector[index] = None
        generate_all_possible_options(original_vector, deepcopy(temp_vector), index + 1, count_positives,
                                      temp_sum, all_options)
    else:
        temp_vector[index] = None
        generate_all_possible_options(original_vector, deepcopy(temp_vector), index + 1, count_positives,
                                      temp_sum, all_options)


def check_negative_sum(indices, lattice):
    """
    Checks if the sum of elements on given indices is negative or not.
    :param indices: given indices
    :param lattice: lattice (matrix of vectors)
    :return:
    """
    limit = sum(filter(None, indices)) - 2
    count = 0
    for row in range(1, len(indices)):
        temp_sum = 0
        index = row
        if count >= limit:
            break
        if row == index and indices[row] == 1:
            count += 1
            temp_sum += lattice[row][row]
            for i in range(index + 1, len(indices)):
                if indices[i] == 1:
                    temp_sum += 2 * lattice[row][i]
        if temp_sum > 0:
            return False
    return True


def count_positive_numbers(vector):
    """
    Count positive elements in given vector (row of matrix).
    :param vector: given vector
    :return: number of positive elements
    """
    count = 0
    for each in vector:
        if each > 0:
            count += 1
    return count


def check_temp_sum(all_options_unsorted, first_row, lattice):
    """
    Checks if we created ok coefficients.
    :param all_options_unsorted: coefficients
    :param first_row: first vector of lattice
    :param lattice: lattice (matrix of vectors)
    """
    for each in all_options_unsorted:
        temp_sum = 0
        for i in range(len(each)):
            if i == 0 and each[i] == 1:
                temp_sum += first_row[0]
            if each[i] == 1 and i != 0:
                temp_sum += 2 * lattice[0][i]
        if temp_sum >= 0:
            print(temp_sum)
            print("error")
            exit(1)


def part_gram_matrix_sum(gram_matrix, indices):
    """
    Calculates norm^2 with the relevant gram matrix and found coefficients indices.
    :param gram_matrix: Gram matrix
    :param indices: given row indices = coefficients
    :return: norm^2 as a sum of relevant g(i,i) elements
    """
    result_sum = 0
    for i in range(len(indices)):
        for j in range(len(indices)):
            if indices[i] == 1 and indices[j] == 1:
                result_sum += gram_matrix[i][j]
    return result_sum


def add_single_vectors(acceptable_vectors, gram_matrix):
    for i in range(len(gram_matrix[0])):
        acceptable_vectors.append(gram_matrix[i][i])


def get_minimal_norm_single_vectors(gram_matrix):
    """
    Finds minimal norm^2 in the lattice with the help of gram matrix.
    Iterates over diagonal elements.
    :param gram_matrix: relevant gram matrix
    :return: minimal norm^2
    """

    min_magnitude = math.inf
    for i in range(len(gram_matrix[0])):
        if gram_matrix[i][i] < min_magnitude:
            min_magnitude = gram_matrix[i][i]
    return min_magnitude


def main_function_help(lattice):
    """
    Finds minimal norm in given lattice, call recursive function.
    Generate all possible linear combinations for the first row given with alpha_i coefficients in {0,1}
    :param lattice: (sub)lattice where to find minimal norm
    :return: return minimal norm found and relevant coefficients
    """

    # sort first vector, decreasing
    first_row_tuple = sorted(enumerate(lattice[0]), key=lambda i: i[1], reverse=True)
    first_row = []
    sort_index = []

    for each in first_row_tuple:
        first_row.append(each[1])
        sort_index.append(each[0])

    # count positives elements
    count_positives = count_positive_numbers(first_row)
    temp_sum = first_row[0]

    temp_vector = [None for _ in range(len(first_row))]
    all_options = []

    # first vector given
    temp_vector[0] = 1

    # generate all possible options that complies conditions
    generate_all_possible_options(deepcopy(first_row), temp_vector, 1, count_positives, temp_sum, all_options)

    all_options_unsorted = []

    for each in all_options:
        index_array = [None for _ in range(len(each))]
        for i in range(len(each)):
            if each[i] == 1:
                index_array[sort_index[i]] = 1
        all_options_unsorted.append(index_array)

    check_temp_sum(all_options_unsorted, first_row, lattice)

    acceptable_indexes = []

    for each in all_options_unsorted:
        if check_negative_sum(each, lattice):
            acceptable_indexes.append(each)

    min_magnitude = math.inf
    min_coefficients = None

    for each in acceptable_indexes:
        magnitude = part_gram_matrix_sum(lattice, each)
        if magnitude < min_magnitude:
            min_magnitude = magnitude
            min_coefficients = each

    return min_magnitude, min_coefficients


def main_function(lattice):
    """
    Call recursive function *main_function_help* on each sub-lattice lattice[i:], i = {1,2,...,n}
    and return  minimal found norm in LLL reduced lattice.
    :param lattice: LLL reduced basis (e.g. output from fpylll wrapper)
    :return: return minimal found norm^2
    """
    min_coeff = None
    min_magnitude = get_minimal_norm_single_vectors(cholesky_decomposition.gram_matrix(lattice))
    for i in range(len(lattice) - 1):
        magnitude, coeff = main_function_help(cholesky_decomposition.gram_matrix(lattice[i:]))
        if magnitude < min_magnitude:
            min_magnitude = magnitude
            min_coeff = coeff
    return min_magnitude, min_coeff
