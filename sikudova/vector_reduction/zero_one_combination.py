def zero_one_combination(power: int):
    """
    Generates array of numbers, -x means subtract x-th vector, x means add x-th vector
    :param power: power of 2
    :return: array of numbers
    """
    # start from 0
    starting_vector = [0 for _ in range(power)]
    starting_vector[len(starting_vector) - 1] = 1

    coeff_array = [power - 1]

    for i in range(1, power):
        starting_vector[power - i - 1] = 1
        starting_vector[power - i] = 1
        coeff_array.append(power - i - 1)
        coeff_array = zero_one_combination_help(starting_vector, 2 ** i, power - 1, coeff_array)

    return coeff_array


def zero_one_combination_help(vector, count: int, index: int, coeff_array):
    """
    Help recursive function.
    :param vector: number, expressed binary
    :param count: length
    :param index: position in array
    :param coeff_array: result array of numbers
    :return:
    """
    index_array = []
    moving_index = index
    sign = 1
    for i in range(count):
        if vector[moving_index] == 1:
            vector[moving_index] = 0
            if i != count - 1:
                sign = -1
        else:
            vector[moving_index] = 1
            if i != count - 1:
                sign = 1
        if i != count - 1:
            coeff_array.append(sign * moving_index)
        if moving_index != index:
            index_array.append(moving_index)
            for idx in range(len(index_array) - 1):
                if index_array[idx] > moving_index:
                    index_array[idx] = None
            index_array = [_ for _ in index_array if _]
            moving_index = index
        else:
            temp_value = index - 1
            if len(index_array) == 0:
                moving_index = index - 1
            else:
                expected_value = temp_value
                for j in range(len(index_array) - 1, -1, -1):
                    if index_array[j] != expected_value:
                        moving_index = expected_value
                        break
                    else:
                        expected_value -= 1
                        moving_index = expected_value

    return coeff_array
