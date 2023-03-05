import numpy as np
import math

def create_matrix(file_h, file_J):
    """
    Convert text file into a numpy array.
    :param file_h: Output of bmDCA called parameters_h_final.txt
    :param file_J: Output of bmDCA called parameters_J_final.txt
    :return: numpy array of the same matrix
    """
    J_parameters_file = open(file_J, "r")
    h_parameters_file = open(file_h, "r")
    J_parameters = J_parameters_file.read().splitlines()
    h_parameters = h_parameters_file.read().splitlines()


    for index, row in enumerate(h_parameters):
        h_parameters[index] = list(map(float, row.split(" ")))

    h_matrix = np.matrix(h_parameters)

    for index, row in enumerate(J_parameters):
        J_parameters[index] = list(map(float, row.split(" ")))

    J_matrix = np.matrix(J_parameters)

    J_parameters_file.close()
    h_parameters_file.close()
    return J_matrix, h_matrix


def save_matrix(J_matrix, h_matrix, path_J, path_h):
    """
    Saving of converted matrix into a .npy file, loading is faster than conversion.
    :param J_matrix: Numpy array created by create_matrix function
    :param h_matrix: Numpy array created by create_matrix function
    :param path_J: Path to .npy file
    :param path_h: Path to .npy file
    :return:
    """
    try:
        np.save(path_J, J_matrix)
        np.save(path_h, h_matrix)
        return True
    except:
        return False


def load_matrix(path_h, path_J):
    """
    Loading of previously saved .npy array for speed and convenience.
    Watch out when debugging to regenerate the matrix, rather than load. The
    toggle create_new_matrix_from_file is made for that purpose.
    :param path_h: Path to .npy file
    :param path_J: Path to .npy file
    :return: numpy arrays
    """
    return np.load(path_h), np.load(path_J)


def create_h_matrix(seq_length, aa_count):
    """
    Mock up of the h_matrix created to ensure the correct cell is accessed.
    For debugging purposes only.
    :param seq_length: sequence length of the msa
    :param aa_count: number of amino acids taken into account, default it 21
    :return: numpy array
    """
    # initialize matrix with given dimensions
    h_matrix = np.zeros((seq_length, aa_count + 1))
    # iterate over rows
    for row in range(0, h_matrix.shape[0]):
        # iterate over column
        h_matrix[row, 0] = row
        for index in range(1, len(h_matrix[row, :])):
            # set value of column to its index
            h_matrix[row, index] = index
    return h_matrix

def create_J_matrix(seq_length, aa_count):
    """
    Mock up of the J_matrix created to ensure the correct cell is accessed.
    For debugging purposes only.
    :param seq_length: sequence length of the msa
    :param aa_count: number of amino acids taken into account, default it 21
    :return: numpy array
    """
    for_matrix = []
    for position_1 in range(0, seq_length - 1):
        for position_2 in range(position_1 + 1, seq_length):
            per_row = []
            per_row.extend([position_1, position_2])
            for amino_acid_a in range(0, aa_count):
                for amino_acid_b in range(0, aa_count):
                    per_row.append(f"{amino_acid_a}-{amino_acid_b}")
            for_matrix.append(per_row)
    return np.array(for_matrix)

def get_J_row(position1, position2, seq_length):
    """
    Identify in which row parameters for a given pair of positions is saved.
    :param position1: <0,seq_length>
    :param position2: <0,seq_length>
    :param seq_length: sequence length of the msa
    :return: Row of matrix, int, <0,number_rows_matrix>
    """
    assert position1 != position2, "positions must differ"
    # only top triangle is filled
    position1, position2 = sorted([position1, position2])[0], sorted([position1, position2])[1]
    assert position2 < seq_length, "out of bounds"
    row_1 = 0
    if position1 > 0:
        for i in range(0, position1):
            row_1 += (seq_length - i - 1)
    row_2 = 0
    if 0 < position2 < seq_length:
        row_2 = position2 - position1 - 1
    row_num = row_1 + row_2
    return row_num


def get_J_column(aa1, aa2, aa_count):
    """
    Identify in which row parameters for a given pair of amino acids is saved.
    :param aa1: <0,aa_count>
    :param aa2: <0,aa_count>
    :param aa_count: default 21
    :return: Column of J_matrix, int, <0,number_columns_matrix>
    """
    base = 2
    column_1 = aa_count * aa1
    column_2 = aa2
    return base + column_1 + column_2


def calculate_energy(h_matrix, J_matrix, sequence, aa_count, T):
    """
    Given the matrices with parameters and a sequence to be evaluated,
    return energy of the sequence.
    :param h_matrix:
    :param J_matrix:
    :param sequence: Sequence of amino acids with length equal to all the others in the msa
    :param aa_count: default 21
    :param T: temperature, as described by "An evolution-based model for designing chorismate mutase enzymes"
    :return: total energy of a given sequence, float
    """
    sum_h = 0
    sum_j = 0
    N = len(sequence)
    for i, aa_i in enumerate(sequence):
        # access location, aa is missing +1 because the conversion function starts at 1 instead of 0
        sum_h += h_matrix[i, aa_i]
        for j, aa_j in enumerate(sequence):
            # upper triangle only
            if i < j:
                # - 1 because the dictionary of aa starts at 1
                sum_j += J_matrix[get_J_row(i, j, N), get_J_column(aa_i - 1, aa_j - 1, aa_count)]
    return (-sum_h - sum_j)/T

