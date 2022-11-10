import numpy as np
import math


def gram_matrix(matrix):
    """
    Calculate Gram matrix of given matrix.
    :param matrix: given matrix
    :return: Gram matrix of matrix
    """
    return np.matmul(matrix, matrix.transpose())


def cholesky_decomposition(matrix):
    """
    Decompose a squared matrix (Gram matrix) into lower and upper triangular matrix.
    :param matrix: squared matrix -- Gram matrix of matrix of vectors
    :return: lower triangular matrix of matrix
    """
    rows, cols = matrix.shape
    lower_tri = np.array([[0.0] * cols] * rows, ndmin=2)

    for row in range(0, rows):
        for col in range(0, rows):
            temp_sum: float = 0
            for k in range(0, col):
                if abs(lower_tri[col, k] - 0) >= 0.01:
                    temp_sum += lower_tri[row, k] * lower_tri[col, k]
            if row == col:
                lower_tri[row, col] = (math.sqrt(abs(matrix[row, row] - temp_sum)))
            else:
                if abs(lower_tri[col, col] - 0) >= 0.01:
                    lower_tri[row, col] = (matrix[row, col] - temp_sum) / lower_tri[col, col]
    return lower_tri
