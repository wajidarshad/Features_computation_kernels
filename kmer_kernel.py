# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 16:22:58 2015

@author: wajidarshad
"""
from Bio.Data import IUPACData 
import numpy as np
import itertools
import string
import random
import time
import matplotlib.pyplot as plt
def random_prot_seq_generator(size):
    prot_dic = dict((k, 0) for k in IUPACData.protein_letters) 
    amino_acids=prot_dic.keys()
    return ''.join(random.choice(amino_acids) for _ in range(size))
def gappy_triplet_kernel(sequences,m,n) :
    t0=time.time()
    data = []
    count=0.0
    for seq in sequences :
        motif_composition = {}
        for pos1 in range(len(seq) - 1) :
            #for offset1 in range(1, min(m, len(seq) - pos1 - 1) + 1) :
                #for offset2 in range(1, min(n, len(seq) - (pos1+offset1+2) - 1) + 1) :
            try:
                motif =  seq[pos1] + seq[pos1 + m+1]+ seq[pos1 + m+n+2] + str(pos1)
            except IndexError:
                break
            try:
                motif_composition[motif] += 1.0
            except KeyError:
                motif_composition[motif] = 1.0

        data.append(motif_composition)
    for key in data[0]:
        try:
            count+=data[0][key]*data[1][key]
        except KeyError:
            count+=0
    t1=time.time()
    print "time taken:",t1-t0
    return count,t1-t0

def pd_kspectrum_kernel(sequences, k) :
    t0=time.time()
    data = []
    count=0.0
    for s in sequences :
        kmers = {}
        for i in range(len(s) - k + 1) :
            kmer = s[i:i+k]
            key = kmer + str(i)
            try:
                kmers[kmer] += 1.0
            except KeyError:
               kmers[kmer] = 1.0
        data.append(kmers)
    for key in data[0]:
        try:
            data[1][key]
            count+=1
        except KeyError:
            count+=0
    t1=time.time()
    print "time taken:",t1-t0
    return count,t1-t0
def kspectrum_kernel(sequences, k) :
    t0=time.time()
    data = []
    count=0
    for s in sequences :
        kmers = {}
        for i in range(len(s) - k + 1) :
            kmer = s[i:i+k]
            try:
                kmers[kmer] += 1.0
            except KeyError:
               kmers[kmer] = 1.0
        data.append(kmers)
    for key in data[0]:
        try:
            count+=data[0][key]*data[1][key]
        except KeyError:
            count+=0
    t1=time.time()
    tt=t1-t0
    print "time taken:",tt
    return count,tt
def count_common_mers(seq1,seq2):
    t0=time.time()
    count = 0
    for letter in set(seq1):
        count += seq2.count(letter)
    t1=time.time()
    print "time taken:",t1-t0
    return count
def count_amino_acids(seq,k): 
     prot_dic = dict((k, 0) for k in IUPACData.protein_letters) 
     amino_acids=prot_dic.keys()
     keys=map(''.join, itertools.product(amino_acids, repeat=k))
     final_dict=dict.fromkeys(keys, 0)
     for aa in final_dict: 
         prot_dic[aa] = seq.count(aa) 
     return prot_dic.values() 

def kmers_kernel(seq1,seq2,k):
    t0=time.time()
    total_counts=0
    for i in range(len(seq1) - k + 1):
        count=seq2.count(seq1[i:i+k])
        total_counts+=count
    t1=time.time()
    print "time taken:",t1-t0
    return total_counts,t1-t0
    
seq =random_prot_seq_generator(100)
seq2=random_prot_seq_generator(10)
#f=count_amino_acids(seq,1)
#g=count_amino_acids(seq2,1)
#print np.inner(f,g)
k_value=kmers_kernel(seq,seq2,1)
print k_value
#d=count_common_mers(seq,seq2)
d,t=kspectrum_kernel([seq,seq2], 1)
print d
print pd_kspectrum_kernel([seq,seq2], 1)
tim=[]
length=[]

for i in range(10000,1000000,100000):
    seq =random_prot_seq_generator(i)
    seq2=random_prot_seq_generator(i)
    #d,t=kspectrum_kernel([seq,seq2], 3)
    #d,t=kmers_kernel(seq,seq2,1)
    d,t=pd_kspectrum_kernel([seq,seq2], 1)
    #d,t=gappy_triplet_kernel([seq,seq2],2,2)
    tim.append(t)
    length.append(i)
    
plt.plot(length,tim)
