import numpy as np

def load_file(filepath):
    fasta = {}
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                if active_sequence_name not in fasta:
                    fasta[active_sequence_name] = []
                continue
            sequence = line
            fasta[active_sequence_name].append(sequence)
    return fasta

def get_sequences(file):
    seqarr = []
    for key,value in enumerate(file):
            seq = ''
            for i in file[value]:
                seq = seq + i
            seqarr.append(seq)
    return seqarr

def get_weight_matrix(sequences):
    # Find the length of the longest sequence
    longest = 0
    for seq in sequences:
        if len(seq) > longest:
            longest = len(seq)
    # Initalize the weight matrix, A=0, C=1, G=2, T=3
    weight_matrix = np.zeros(shape=(longest,4))
    # Populate the weight matrix
    for seq in sequences:
        # For each nucleotide position
        for n in range(len(seq)):
            if(seq[n] == 'A'): weight_matrix[n][0]+=1
            if(seq[n] == 'C'): weight_matrix[n][1]+=1
            if(seq[n] == 'G'): weight_matrix[n][2]+=1
            if(seq[n] == 'T'): weight_matrix[n][3]+=1

    for n in range(len(weight_matrix)):
        sum_of_position = 0
        for x in range(len(weight_matrix[n])):
            sum_of_position +=weight_matrix[n][x]
        for x in range(len(weight_matrix[n])):
            weight_matrix[n][x] =weight_matrix[n][x]/sum_of_position
    print weight_matrix

def main():
    file = load_file("GATA2_chr1.fa")
    sequences = get_sequences(file)
    weight_matrix = get_weight_matrix(sequences)
    


main()