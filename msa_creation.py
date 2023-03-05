sequences, names = open("/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/sequences.txt", "r"), open("/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/seq_names.txt", "r")
msa = open("/Users/yanbarta/MUNI/Thesis/bmDCA_reproduce/msa.fasta", "w")
counter = 0
for seq, name in zip(sequences.readlines(), names.readlines()):
    counter += 1
    msa.write(f">{counter}_{name}{seq}")
sequences.close()
names.close()
msa.close()

