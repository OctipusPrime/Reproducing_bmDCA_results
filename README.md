# Reproducing_bmDCA_results
 An attempt to reproduce results shown in "An evolution-based model for designing chorismate mutase enzymes" that was made in order to evaluate the viability of using bmDCA to predict the melting temperature and solubility of a protein sequence. 

## Process

1. Sequences were gathered from the provided supplementary table "aba3304_table_s1.xlsx" using the msa_creation.py script.
2. Energies of sequences were calculated using the main.py, parameter_conversion.py and sequence_handling.py scripts. 
3. Energies were calculated using the arDCA package.
4. Results were compared using a pairplot.

## Results

As shown in the pariplot below, no correlation was found between values provided by the authors, energies calculated with my scripts and energies given by arDCA on the same dataset. I have not been able to identify what could be causing this apparent unrelatedness of these values, which should theoretically capture the same information. Please let me know if you see where the problem lies. 

This approach and project were subsequently terminated. 