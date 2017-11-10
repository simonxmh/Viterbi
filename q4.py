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
    all_possible_sequences = []
    seq = ''
    for pc1 in pc[0]:
        for pc2 in pc[1]:
            for pc3 in pc[2]:
                for pc4 in pc[3]:
                    for pc5 in pc[4]:
                        for pc6 in pc[5]:
                            seq = pc1+pc2+pc3+pc4+pc5+pc6
                            all_possible_sequences.append(seq)
    
    for possible_sequence in all_possible_sequences:
        for s in range(0,len(sequence)-6):
            if sequence[s:s+6] == possible_sequence:
                print "Found sequence match: " + possible_sequence
                return True
                
    return False




def main():
    file = load_file("GATA-sample.fa")
    sequences = get_sequences(file)
    alphabet = [['A'],['C'],['G'],['T'],['A','G'],['C','T'],['A','C','G','T']]
    permutations = get_permutations(alphabet)
    scores = {}

    time1 = time.time()
    # for p in permutations:
    #     key =''.join(str(v) for v in p)
    #     scores[key] = calculate_n_of_c(p,sequences)
    time2 = time.time()
    print (time2-time1)

    
main()

    # print located("ACGTGCCCA",[['A'],['A','C','G','T'],['G'],['A','G'],['G'],['C'],['C']])
