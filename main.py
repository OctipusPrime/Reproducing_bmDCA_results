import csv
import parameter_conversion as pc
import matplotlib.pyplot as plt
import sequence_handling as sh
from Bio import SeqIO

file_J = "/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/paper_output/parameters_J_final.txt"
file_h = "/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/paper_output/parameters_h_final.txt"
path_J = "/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/J_matrix.npy"
path_h = "/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/h_matrix.npy"
dictionary_path = "/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/scores_dic.npy"
dictionary_csv = "/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/scores_dic.csv"
# Use SeqIO to read the .fasta file
fasta_sequences = SeqIO.parse(open('/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/msa.fasta'), 'fasta')
seq_length = len(list(fasta_sequences))
aa_count = 21
T_0 = 1
T_1 = 0.33
T_2 = 0.66
# Sequence that acts as a reference, it's energy set as 0
ref = "-TSENPLLALREKISALDEKLLALLAERRELAVEVGKAKLLSHRPVRDIDRERDLLERLITLGK-AHHLDAHYITRLFQLIIEDSVLTQQALLQQH"
create_new_matrix_from_file = True
calculate_enery_for_sequences = True

if __name__ == '__main__':
    if create_new_matrix_from_file:
        h_matrix, J_matrix = pc.create_matrix(file_J, file_h)
        pc.save_matrix(J_matrix, h_matrix, path_J, path_h)
    else:
        h_matrix, J_matrix = pc.load_matrix(path_h, path_J)

    ref_sequence = pc.calculate_energy(h_matrix, J_matrix, sh.conversion(ref), aa_count, T_0)

    # Extract the sequences from the .fasta file and print them
    if calculate_enery_for_sequences:
        scores = sh.score_sequences(fasta_sequences, h_matrix, J_matrix, aa_count, dictionary_path, T_0, ref_sequence)
    else:
        scores = sh.load_sequences(dictionary_path)

    # Write gathered energies into a csv file
    with open(dictionary_csv, 'w') as f:
        w = csv.writer(f)
        w.writerow(['species', 'value'])
        for key, value in scores.items():
            w.writerow([key, value])
