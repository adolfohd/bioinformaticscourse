
from Bio import motifs
import random

w = 4

seqstmp = open("sequences.txt").readlines()
sequences = [x.strip() for x in seqstmp]

T = len(sequences)
N = len(se)

# mt = motifs.create(sequences)

print "> Motivo:\n", mt

background = {"A":0.25 ,"C":0.25 ,"G": 0.25 ,"T":0.25 }
alineamientos = []
for seq in sequences:
    n = len(seq)
    i = random.radint(0, n - w)
    subseq = seq[i:i+w]
    alineamientos.append(subseq)

mt_aln = motif.create(alineamientos)
print ">Initial random motifs: \n", mt_aln

random_index = randrange(0,len(sequences))

for i in random_index:
    list_sh = []    
    random_seq = sequences[i]
    N = len(random_seq)
    for j in range(0, N - w):
        subseq = random_seq[j:j+w]
        list_sh.append(subseq)
        score = log_odd(pwm, background, subseq)
        if score > maxScore:
            maxIndice = j
            maxScore = score
            maxSubseq = subseq

        
        

