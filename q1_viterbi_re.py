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

def write_file(data,filepath):
    with open(filepath, "w") as file:
        for key,value in data.iteritems():
            seq = ''
            file.write('>' + key + ':\n')
            # print value
            for i in range(len(value)):
                seq = seq+value[i]
            file.write(seq + '\n')
    file.close

def get_sequences(file):
    seqarr = {}
    for key,value in file.iteritems():
        seq = ''
        for i in file[key]:
            seq = seq + i
        seqarr[key] = seq
    return seqarr

def init_prob_matrix():
    hydrophobic_aa= ['A','V','I','L','M','F','Y','W']
    hydrophilic_aa= ['R','H','K','D','E','S','T','N','Q','C','G','P']

    # Initialize emission probability matrix
    emission_prob={}
    emission_prob['hydrophobic']={}
    emission_prob['hydrophilic']={}
    emission_prob['mixed']={}
    for b in hydrophobic_aa:
        for l in hydrophilic_aa:
            emission_prob['hydrophobic'][b] = float(3)/40
            emission_prob['hydrophobic'][l] = float(1)/30
            emission_prob['hydrophilic'][b] = float(1)/40
            emission_prob['hydrophilic'][l] = float(1)/15
    for m in hydrophobic_aa + hydrophilic_aa:
        emission_prob['mixed'][m]=float(1)/20

    # Initialize transition probability matrix
    transition_prob = {}
    transition_prob['hydrophobic']={}
    transition_prob['hydrophilic']={}
    transition_prob['mixed']={}

    transition_prob['hydrophobic']['hydrophobic'] = 0.80
    transition_prob['hydrophobic']['hydrophilic'] = 0.04
    transition_prob['hydrophobic']['mixed'] = 0.16
    transition_prob['hydrophilic']['hydrophilic'] = float(7)/8
    transition_prob['hydrophilic']['hydrophobic'] = 0.0375
    transition_prob['hydrophilic']['mixed'] = 0.0875
    transition_prob['mixed']['mixed']=float(6)/7
    transition_prob['mixed']['hydrophobic']=float(1)/14
    transition_prob['mixed']['hydrophilic']=float(1)/14

    return emission_prob, transition_prob

def viterbi_algo(emission,transition,seq):
    # Initialize viterbi matrix
    viterbi = {}
    viterbi[0]={}
    
    viterbi[0]['hydrophobic']=(math.log(emission['hydrophobic'][seq[0]])+math.log(float(1)/3),None)
    viterbi[0]['hydrophilic']=(math.log(emission['hydrophilic'][seq[0]])+math.log(float(1)/3),None)
    viterbi[0]['mixed']=(math.log(emission['mixed'][seq[0]])+math.log(float(1)/3),None)

    output = []

    for i in range(1,len(seq)):
        viterbi[i] = {}
        max_hydrophobic=get_max(viterbi,emission,transition,i,'hydrophobic')
        max_hydrophilic=get_max(viterbi,emission,transition,i,'hydrophilic')
        max_mixed=get_max(viterbi,emission,transition,i,'mixed')
        viterbi[i]['hydrophobic']=(max_hydrophobic[0][0]+math.log(emission['hydrophobic'][seq[i]]),max_hydrophobic[0][1])
        viterbi[i]['hydrophilic']=(max_hydrophilic[0][0]+math.log(emission['hydrophilic'][seq[i]]),max_hydrophilic[0][1])
        viterbi[i]['mixed']=(max_mixed[0][0]+math.log(emission['mixed'][seq[i]]),max_mixed[0][1])


    return viterbi



def get_max(viterbi,emission,transition,index,current):
    max_hydrophobic_prev = viterbi[index-1]['hydrophobic'][0]+math.log(transition['hydrophobic'][current])
    max_hydrophilic_prev = viterbi[index-1]['hydrophilic'][0]+math.log(transition['hydrophilic'][current])
    max_mixed_prev = viterbi[index-1]['mixed'][0]+math.log(transition['mixed'][current])
    # the smallest log is actually the largest number
    max_state = max(max_hydrophobic_prev,max_hydrophilic_prev,max_mixed_prev)
    if max_state == max_hydrophobic_prev:
        return ((max_hydrophobic_prev,'B'),viterbi[index-1]['hydrophobic'])
    if max_state == max_hydrophilic_prev:
        return ((max_hydrophilic_prev,'L'),viterbi[index-1]['hydrophilic'])
    if max_state == max_mixed_prev:
        return ((max_mixed_prev,'M'),viterbi[index-1]['mixed'])

def get_length_distributions(sequences):
    length_dist = {}
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
        lengths[protein]=((num_ho,num_hi,num_mi),len(sequence))
    return lengths

def traceback(viterbi,seq):
    output = []
    i = len(seq)-1
    max_final = max(viterbi[i]['hydrophobic'],viterbi[i]['hydrophilic'],viterbi[i]['mixed'])
    
    if max_final == viterbi[i]['hydrophobic']:
        while i>0:
            output.append(viterbi[i]['hydrophobic'][1])
            i-=1
    elif max_final == viterbi[i]['hydrophilic']:
        while i>0:
            output.append(viterbi[i]['hydrophilic'][1])
            i-=1
    elif max_final == viterbi[i]['mixed']:
        while i>0:
            output.append(viterbi[i]['mixed'][1])
            i-=1
        
    return output[::-1]

def main():
    file = load_file("test.fa")
    sequences = get_sequences(file)
    emission_prob, transition_prob = init_prob_matrix()

    output = {}

    for protein,seq in sequences.iteritems():
        viterbi_matrix = viterbi_algo(emission_prob,transition_prob,seq)
        output[protein] = traceback(viterbi_matrix,seq)
    # print viterbi_matrix

    length_dist = get_length_distributions

    write_file(output,'q1_re_out.txt')
main()