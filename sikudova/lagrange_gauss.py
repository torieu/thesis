import numpy as np


def lagrange_gauss_reduction(basis):
    """
    Perform Lagrange-Gauss algorithm on 2D-basis
    :param basis: basis = (b1, b2)
    :return: reduced basis = (b1, b2)
    """

    magnitude1 = np.linalg.norm(basis[0])
    magnitude2 = np.linalg.norm(basis[1])

    # vectors swapping, vector1 is the smaller one, vector2 is the bigger one
    if magnitude1 < magnitude2:
        vector1, vector2 = basis[0], basis[1]
    else:
        vector2, vector1 = basis[0], basis[1]

    magnitude1 = np.linalg.norm(vector1)
    magnitude2 = np.linalg.norm(vector2)

    # basis is lagrange reduced iff: |b1| <= |b2| <= |b2 +- b1|
    while not (magnitude1 <= magnitude2 <= np.linalg.norm(vector2 + vector1)) or not (
            magnitude1 <= magnitude2 <= np.linalg.norm(vector2 - vector1)):
        coefficient = round(magnitude2 / magnitude1)
        vector2 = vector2 - coefficient * vector1
        magnitude1 = np.linalg.norm(vector1)
        magnitude2 = np.linalg.norm(vector2)

        if magnitude1 < magnitude2:
            vector1, vector2 = vector1, vector2
        else:
            vector2, vector1 = vector1, vector2

    return np.array([vector1, vector2])


def lagrange_gauss_reduction_2(basis):
    """
    Perform Lagrange-Gauss algorithm on 2D-basis.
    :param basis: basis = (b1, b2)
    :return: reduced basis = (b1, b2)
    """
    vector1 = basis[0]
    vector2 = basis[1]

    magnitude1 = np.linalg.norm(vector1) ** 2
    coefficient = np.dot(vector1, vector2) / magnitude1
    vector2 = vector2 - round(coefficient) * vector1
    magnitude2 = np.linalg.norm(vector2) ** 2

    while magnitude2 < magnitude1:
        vector1, vector2 = vector2, vector1  # swapping vectors
        magnitude1 = magnitude2
        coefficient = np.dot(vector1, vector2) / magnitude1
        vector2 = vector2 - round(coefficient) * vector1  # reduction of bigger vector
        magnitude2 = np.linalg.norm(vector2) ** 2

    return np.array([vector1, vector2])
