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
    seqarr = {}
    for key,value in file.iteritems():
        seq = ''
        for i in file[key]:
            seq = seq + i
        seqarr[key] = seq
    return seqarr

def get_largest_mixed(sequences):
    mixed = {}
    for protein,seq in sequences.iteritems():
        max = 
        count = 0
        for i in range(len(seq)):
            if seq[i] == 'M':
                count+=1
        mixed[protein]=float(count)/len(seq)
        if mixed[protein] > mixed[max]:
            max = mixed[protein]
    print max
    return max(mixed)



def main():
    file = load_file("q1_out.txt")
    sequences = get_sequences(file)
    
    largest_mixed = get_largest_mixed(sequences)

main()