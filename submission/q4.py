import time
import math

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

def get_possible_consensus(alphabet):
    all_possible_consensus = []
    count = 0
    print 'finding all possible combinations of consensus'
    for p1 in alphabet:
        for p2 in alphabet:
            for p3 in alphabet:
                for p4 in alphabet:
                    for p5 in alphabet:
                        for p6 in alphabet:
                            possible_consensus= (p1,p2,p3,p4,p5,p6)

                            # no-filter out general sequences
                            # all_possible_consensus.append(possible_consensus)
                            # count +=1

                            # filter out general sequences
                            length = 0
                            for p in possible_consensus:
                                length += len(p)
                            if(length <= 12):
                                all_possible_consensus.append(possible_consensus)
                                count +=1
    print 'found all possible consensus: ' + str(count)

    return all_possible_consensus

def get_prob(consensus,p_a,p_c,p_g,p_t):
    probability= 1
    for element in consensus:
        element_prob= 0
        for i in element:
            if i=='A':
                element_prob+=p_a
            if i=='C':
                element_prob+=p_c
            if i=='G':
                element_prob+=p_g
            if i=='T':
                element_prob+=p_t
        probability*=element_prob
    # print probability
    return probability

def get_all_combinations(all_consensus):
    possible = {}
    for consensus in all_consensus:
        for nt1 in consensus[0]:
            for nt2 in consensus[1]:
                for nt3 in consensus[2]:
                    for nt4 in consensus[3]:
                        for nt5 in consensus[4]:
                            for nt6 in consensus[5]:
                                nt_tuple = (nt1,nt2,nt3,nt4,nt5,nt6)
                                if nt_tuple in possible:
                                    possible[nt_tuple].append(consensus)
                                else:
                                    possible[nt_tuple] = [consensus]
    return possible



def calculate_z_of_c(scores, total_random_sequence_length,p_a,p_c,p_g,p_t):
    score = scores
    e_score = {}

    for key,value in enumerate(score):
        e_score[value] = total_random_sequence_length*get_prob(value,p_a,p_c,p_g,p_t)
        # print e_score[value]
  

    for key,value in enumerate(e_score):
        if value in score:
            score[value] = float(score[value])-e_score[value]/math.sqrt(e_score[value])
            # print score[value]

    
    return score

def calculate_n_of_c(sequence, scores, all_possible_combinations):
    score = scores
    for s in range(0,len(sequence)-6):

        char_tuple = ()
        for el in sequence[s:s+6]:
            char_tuple += tuple(el)

        all_combinations = all_possible_combinations[char_tuple]

        if char_tuple in score:
            for combo in all_combinations:
                score[combo] += 1
        else:
            for combo in all_combinations:
                score[combo] = 1
    return score
    


def main():
    gata_file = load_file("GATA2_chr1.fa")
    not_gata_file = load_file("not_GATA2_chr1.fa")
    sequences = get_sequences(gata_file)
    random_sequences = get_sequences(not_gata_file)

    alphabet = [('A'),('C'),('G'),('T'),('A','G'),('C','T'),('A','C','G','T')]

    all_possible_consensus = get_possible_consensus(alphabet)
    all_possible_combinations = get_all_combinations(all_possible_consensus)

    scores = dict.fromkeys(all_possible_consensus,0)
    total_random_sequence_length = 0

    p_a=0
    p_c=0
    p_g=0
    p_t=0
    for i in random_sequences:
        for letter in i:
            if letter =='A':
                p_a +=1
            if letter =='C':
                p_c +=1
            if letter =='G':
                p_g +=1
            if letter =='T':
                p_t +=1
        total_random_sequence_length += len(i)-6+1
    p_a = float(p_a)/total_random_sequence_length
    p_c = float(p_c)/total_random_sequence_length
    p_g = float(p_g)/total_random_sequence_length
    p_t = float(p_t)/total_random_sequence_length

    

    time1 = time.time()
    num_seq = 0
    for s in sequences:
        scores = calculate_n_of_c(s,scores,all_possible_combinations)
        num_seq+=1
        print 'number of sequences' + str(num_seq)
    scores = calculate_z_of_c(scores,total_random_sequence_length,p_a,p_c,p_g,p_t)
 
    time2 = time.time()
    print time2-time1
    # print scores


    max = scores[list(scores.keys())[0]]
    max_key = ''
    for key,value in scores.iteritems():
        if value < max:
            max = scores[key]
            max_key = key

    print max_key, max
    
main()

    # print located("ACGTGCCCA",[['A'],['A','C','G','T'],['G'],['A','G'],['G'],['C'],['C']])
