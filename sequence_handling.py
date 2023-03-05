import parameter_conversion as pc
import numpy as np

def conversion(sequence):
    """
    Transformation of letter encoded protein sequence into a int coding, which is required
    by bmDCA
    :param sequence: sequence of aa letters of length seq_length
    :return: list of ints representing the sequence
    """
    sequence = sequence.upper()
    alphabet = {"A":1, "C":2, "D":3, "E":4, "F":5, "G":6, "H":7, "I":8, "K":9, "L":10, "M":11, "N": 12, "P":13, "Q":14, "R":15, "S":16, "T":17, "V":18, "W":19, "Y": 20, "-":21}
    # check for non-protein letters
    assert not any(letter not in alphabet.keys() for letter in sequence), "Invalid sequence"
    for letter in sequence:
        if letter not in alphabet.keys():
            print(letter)
    numerical_sequence = []
    for letter in sequence:
        numerical_sequence.append(alphabet[letter])
    # convert sequence into numerical array
    return numerical_sequence


def score_sequences(fasta_sequences, h_matrix, J_matrix, aa_count, dictionary_path, T, reference):
    """
    Return dictionary of score of all the sequences in the MSA.
    :param fasta_sequences: SeqIO List of sequences (containing .id and .seq)
    :param h_matrix: Numpy array
    :param J_matrix: Numpy array
    :param aa_count: default 21
    :param dictionary_path: where to save
    :param T: temperature needed for calculate_energy, see bmDCA paper for more details
    :param reference: int, energy to normalize to
    :return: dictionary of id:energy
    """
    scores = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        scores[name] = pc.calculate_energy(h_matrix, J_matrix, conversion(sequence), aa_count, T) - reference
    np.save(dictionary_path, scores)
    return scores


def load_sequences(dictionary_path):
    return np.load(dictionary_path, allow_pickle=True)