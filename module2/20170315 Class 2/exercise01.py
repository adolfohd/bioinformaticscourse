#!/usr/bin/python

from Bio import motifs
# seqsTmp = open ()
sequences = []

for seq in open("sequences.txt"):
    sequences.append(seq.strip())

mt = motifs.create(sequences)

# print mt
# mt.counts
# mt.counts.normalize
# mt.pssm
# mt.consensus

background = {"A":0.25 ,"C":0.25 ,"G": 0.25 ,"T":0.25 }
# mtn = mt.counts.normalize
pwm = mt.pwm



S = "GGAGTAAAC"
score = 1

for k,i in enumerate(S):
    score = score * ( pwm [i] [k] + 0.01 ) # pseudo-count against the killer zeroes in the pwm matrix

print "The score for ", S , " is " , score

scoreBackground = 0.25 ** 9

log_odd = score / scoreBackground