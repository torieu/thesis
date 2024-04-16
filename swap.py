from generate import *
from closest_vs_minimal import *
from sage.all import *


def find_real_minimum_genfile(G, current_row, inserted_component) -> vector:
    """
    returns linear combination of the minimum over real numbers
    the >current_row< parameter is the coordinate of fixed col/row idk
    the >K< parameter is the number inserted into this column
    """
    G = matrix_to_list(G)
    matrixA = matrix(dimension - 1, dimension - 1, 0) # square matrix of size (dimension - 1) x (dimension - 1), filled with zeros.
    matrixB = matrix(dimension - 1, 1, 0) # column matrix of size (dimension - 1) x 1, filled with zeros.
    # matrixA[0,0] = 1
    print("G",matrix(G))
    a, b = 0, 0
    for row in range(dimension):
        if row != current_row:
            matrixA[a] = [G[row][j] for j in range(len(G[row])) if j != current_row]
            matrixB[b] = sum([inserted_component * G[row][j] for j in range(len(G[row])) if j == current_row])
            a += 1
            b += 1
    # insert indices
    print("A", matrixA,"\n","b", -(1)*matrixB, )
    result = vector_to_list(matrixA.solve_right((-1)*(inserted_component) * matrixB))
    result.insert(current_row-1, inserted_component)
    print(result)
    return vector(result).n(digits=5)


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# %run closest_versus_minimal_tools.ipynb




def swapping(dimension, perimeter):
    basis_as_matrix = (randomMatrix(dimension, perimeter)) #  matrix
#     basis_as_matrix= matrix([[5, 3, 4], [-2, 4, 4], [-10, -3, -7]])
    # print(basis_as_matrix)
    basis_as_list= matrix_to_list(basis_as_matrix)
    print(basis_as_list)
    used_vectors =[]
    used_vectors.extend(basis_as_list)

    G = gram_matrix(basis_as_matrix)

    basis_as_list.sort(key=lambda x: euclidean_norm(x))
    print("Basis na začatku: ")
    for i in range(dimension):
            print("  b_",i, " ", euclidean_norm(basis_as_list[i]), basis_as_list[i])
    print("--------------")
    basis_as_matrix=matrix(basis_as_list)
    # print(basis_as_matrix)
    # print(basis_as_list)



    for i in range(0,dimension):
        print("\n\n")
        print(i)
        real_min = find_real_minimum_genfile(G, i, 1)

        # print(evaluate_norms_at_cube(cube_points(real_min), basis_as_list))
        
        shortestLC = shortest_lc_in_cube(real_min, basis_as_list)
        # print(real_min)
        # print(shortestLC)
        shortest_vector = vector(shortestLC)*basis_as_matrix
        # print("sv", shortest_vector,(shortest_vector).norm().n())
        print("used vectors",used_vectors)
        if vector_to_list(shortest_vector) in used_vectors:
              print("ERROR: I already used the vector ", shortest_vector)
              continue
        if shortest_vector == zero_vector(SR, dimension):
              print("ERROR: Zero vector")
              continue
        shortest_vector = vector_to_list(shortest_vector)
        used_vectors.append(shortest_vector)
        print("basis to be thrown away", basis_as_list[dimension-1])
        basis_as_list[dimension-1] = shortest_vector
        print(basis_as_list)
        print("basis new, should be SV", basis_as_list[dimension-1])
        basis_as_list.sort(key=lambda x: euclidean_norm(x))
        print(basis_as_list)
        for i in range(dimension):
                print("  b_",i, " ", euclidean_norm(basis_as_list[i]), basis_as_list[i])
        basis_as_matrix = matrix(basis_as_list)
        G = gram_matrix(basis_as_matrix)

dimension = 3 
swapping(dimension, 10)


# def swapping_alg_ver2(dimension, perimeter):
#     """ Verse s for cyklem """
#     basis = randomMatrix(dimension, perimeter)
#     i=0
#     while ještě pořád dělám rozumne změny:
#         perform a cube with fixed index i
#         realmin = vysledek cube tedy něco jkao [1, 0.54, 7.2345]
#         sv = zaokrouhleni realmin * B
#         basis[i] = sv
#         i= i + 1 mod dimension
#         basis.sort by length



# CO KDYŽ VE VĚTŠÍCH DIMENZÍCH SE STREFUJU DO VEKTORU KRATŠÍHO NEŽ LLL???? JE TO POSSIBLE !!!! ZKONTROLOVAT
    
# možná se ještě vyplatí opravit cube takto: co když sii jako real-min uložím něco, co ale je moc daleko a  
# fixací JINÉ souřadnice ale dostanu už správnej výsledek? ještě na to mrkni viko