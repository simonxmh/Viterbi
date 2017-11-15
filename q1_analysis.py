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
    max_mixed_protein = ''
    key = ''
    for protein,seq in sequences.iteritems():
        count = 0
        for i in range(len(seq)):
            if seq[i] == 'M':
                count+=1
        mixed[protein]=float(count)/len(seq)
    
    for protein,count in mixed.iteritems():
        print protein,count
        if mixed[protein] >= mixed[max]:
            max_mixed_protein = mixed[protein]
            key = protein
    print key,max
    return max(mixed)

def get_length_distributions(sequences):
    lengths = {}
    
    protein_id = 0
    for protein,sequence in sequences.iteritems(): 
        num_ho = 0
        num_hi = 0
        num_mi = 0  
        for el in sequence:
            if el == 'L':
                num_hi+=1
            elif el == 'B':
                num_ho+=1
            elif el =='M':
                num_mi+=1
        lengths[protein_id]=float(num_ho)/len(sequence),float(num_hi)/len(sequence),float(num_mi)/len(sequence)
        protein_id+=1
    return lengths


def get_aa_distributions(sequences,og_sequences):
    dist_l={}
    dist_b={}
    dist_m={}
    for protein,sequence in sequences.iteritems():
        protein = protein[:-1]
        for i in range(len(sequence)):
            if sequence[i] == 'L':
                if og_sequences[protein][i] in dist_l:
                    dist_l[og_sequences[protein][i]] +=1
                else:
                    dist_l[og_sequences[protein][i]] =1
            elif sequence[i] == 'B':
                if og_sequences[protein][i] in dist_b:
                    dist_b[og_sequences[protein][i]] +=1
                else:
                    dist_b[og_sequences[protein][i]] =1
            elif sequence[i] =='M':
                if og_sequences[protein][i] in dist_m:
                    dist_m[og_sequences[protein][i]] +=1
                else:
                    dist_m[og_sequences[protein][i]] =1
    # print dist_b
    return dist_b,dist_l,dist_m


def main():
    file = load_file("q1_out.txt")
    og_file = load_file("hw3_proteins.fa")
    sequences = get_sequences(file)
    og_sequences = get_sequences(og_file)
    # largest_mixed = get_largest_mixed(sequences)
    
    distributions_length = get_length_distributions(sequences)
    print distributions_length

    # dist_b,dist_l,dist_m = get_aa_distributions(sequences,og_sequences)
    # print dist_b,dist_l,dist_m

    

    outfile = open('dict.txt', 'w' )
    for key, value in sorted(distributions_length.items() ):
        outfile.write( str(key) + '\t' + str(value) + '\n' )


main()