import random
from Bio import motifs


######################## Functions ###########################
def log_odd(pwmMatrix, back, s):
    scoreModelo = 1
    for k, i in enumerate (s):
        scoreModelo = scoreModelo * (pwmMatrix[i][k] + 0.01) / back[i]  
            # plus pseudo-count against wicked zeroes.
            # Divided by its background probability
    return scoreModelo

def actualizar_pwm(alineamientos, subseq, row, column):
    alineamientos[row] = subseq # change row
    pwm = motifs.create(alineamientos).pwm
    return alineamientos, pwm

####################### Parameters ############################

subseq_length = 4
epsilon = 0.001
sequencesFileName = "sequences.txt"

###################################################
seqstmp = open(sequencesFileName).readlines()
sequences = [x.strip() for x in seqstmp]

T = len(sequences)
N = len(sequences[1])

print " \n \n ", T, " sequences of size ", N, " loaded:\n", sequences, "\n"

# mt = motifs.create(sequences)
# print "> Motivo:\n", mt
background = {"A":0.25 ,"C":0.25 ,"G": 0.25 ,"T":0.25 }
alineamientos = []
for seq in sequences:
    seq_length = len(seq)
    subseq_start = random.randint(0, seq_length - subseq_length)
    subseq = seq[subseq_start:subseq_start+subseq_length]
    alineamientos.append(subseq)

mt_aln = motifs.create(alineamientos)
print "> Initial random motifs: \n", mt_aln

pwm = mt_aln.pwm
print "> Position weight matrix: \n", pwm

random_index = (range(0, len(sequences)))
random.shuffle(random_index)

print "random indexes: \n", random_index

maxScore = -1
score = maxScore

class AlgorithmConvergenceReached(Exception): pass
for row in random_index:
    list_sh = []
    random_seq = sequences[row]
    N = len(random_seq)
    
    try:
        for column in range(0, N - subseq_length):
            subseq = random_seq[column:column+subseq_length]

            if (subseq != alineamientos[row]):
                list_sh.append(subseq)
                new_score = log_odd(pwm, background, subseq)
                #  print "The new score is ", new_score
                if new_score > maxScore:
                    # maxIndice = j
                    
                    print "> Updated max score is: ", new_score, " replacing ", maxScore
                    maxScore = new_score
                    # maxSubseq = subseq
                    alineamientos, pwm = actualizar_pwm(alineamientos, subseq, row, column)
                    print "> Better motif found changing row # ", row, "with the sequence #", column
                    print "> New motifs: \n", alineamientos
                    print "> New PWM: \n", pwm
                    delta = abs(new_score - score)
                    print "> Delta error = ", delta , "\n"
                    if abs(new_score - score) < epsilon:
                        raise AlgorithmConvergenceReached
                    else:
                        score = new_score
            else:
                print "The sub-sequence ", subseq, " was already evaluated for row ", row
    except AlgorithmConvergenceReached:
        print "> Algorithm converged!"
        exit()

print "Algorithm never converged"




        

                


        
        

