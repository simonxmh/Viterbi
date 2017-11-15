import time

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

def get_permutations(alphabet):
    permutations = []
    count = 0
    print 'finding all possible combinations of consensus'
    for p1 in alphabet:
        for p2 in alphabet:
            for p3 in alphabet:
                for p4 in alphabet:
                    for p5 in alphabet:
                        for p6 in alphabet:
                            permutation = (p1,p2,p3,p4,p5,p6)
                            length = 0
                            for p in permutation:
                                length += len(p)
                            if(length <= 12):
                                permutations.append(permutation)
                                count +=1
    print 'found all possible combinations: ' + str(count)

    return permutations

def calculate_n_of_c(permuation,sequences):
    score = 0
    for seq in sequences:
        is_located = located(seq,permuation)
        if is_located == True:
            score += 1
    return score

def located(sequence,pc):
    all_possible_consensus = []
    seq = ''
    for pc1 in pc[0]:
        for pc2 in pc[1]:
            for pc3 in pc[2]:
                for pc4 in pc[3]:
                    for pc5 in pc[4]:
                        for pc6 in pc[5]:
                            seq = pc1+pc2+pc3+pc4+pc5+pc6
                            all_possible_consensus.append(seq)
    
    for possible_sequence in all_possible_consensus:
        for s in range(0,len(sequence)-6):
            if sequence[s]== possible_sequence[0]:
                if sequence[s+1]== possible_sequence[1]:
                    if sequence[s+2]== possible_sequence[2]:
                        if sequence[s+3]== possible_sequence[3]:
                            if sequence[s+4]== possible_sequence[4]:
                                if sequence[s+5]== possible_sequence[5]:
                                    return True
    return False




def main():
    gata_file = load_file("GATA2_chr1.fa")
    not_gata_file = load_file("not_GATA2_chr1.fa")
    sequences = get_sequences(gata_file)
    random_sequences = get_sequences(not_gata_file)

    alphabet = [['A'],['C'],['G'],['T'],['A','G'],['C','T'],['A','C','G','T']]
    permutations = get_permutations(alphabet)
    max = []
    scores = {}
    random_scores = {}

    time1 = time.time()
    for p in permutations:
        key =''.join(str(v) for v in p)
        scores[key] = calculate_n_of_c(p,sequences)
        random_scores[key] = calculate_n_of_c(p,random_sequences)
    time2 = time.time()
    
    final = {}

    for key,value in scores.iteritems():
        # Z-score = N of C - E of C over sqrt(E of C)
        final[key] = (score[key]-random_scores[key])/math.sqrt(random_scores[key]) 

    max = final[list(final.keys())[0]]
    max_key = ''
    for key,value in final.iteritems():
        if value > max:
            max = final[key]
            max_key = key

    return key
    
main()

    # print located("ACGTGCCCA",[['A'],['A','C','G','T'],['G'],['A','G'],['G'],['C'],['C']])
