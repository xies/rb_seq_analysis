#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 19:53:23 2017

@author: mimi
"""

import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from collections import Counter
from scipy import stats


def rec2dict(fasta):
    # convert a list of records from SeqIO into a dictionary using .id
    seqs = {}
    for f in fasta:
        seqs[f.id] = f
        
    return seqs
    
####

def randomize_position(pos):
    # Find non indel species
    I = np.where( pos != '-' )
    Iperm = np.permute(I)
    return pos

def mutual_information(seqs,i,j):
     
    seqs = [s for s in seqs if not '-' in [s[i],s[j]]]
    # Get position frequency table, pseudocount (+1), and normalize
    Pi = Counter(s[i] for s in seqs)
    Pi = add_all_AAs(Pi)
    Pi = add_pseudocount(Pi)
    Pi = normalize(Pi)

    # Other position
    Pj = Counter(s[j] for s in seqs)
    Pj = add_all_AAs(Pj)
    Pj = add_pseudocount(Pj)
    Pj = normalize(Pj)
    
    # Joint probability
    Pij = Counter((s[i],s[j]) for s in seqs)
    Pij = add_all_AAs_2d(Pij)
    Pij = add_pseudocount(Pij)
    Pij = normalize(Pij)
    
    return sum(Pij[(x,y)] * \
                  log(  np.float(Pij[(x,y)]) / np.float(Pi[x]*Pj[y]) ) \
                  for x,y in Pij)

def add_all_AAs(Pi):
    AAs = 'GAVLIPFYWSTCMNQKRHDE'
    for a in AAs:
        if not Pi.has_key(a):
            Pi[a] = 0
    return Pi
    
def add_pseudocount(Pi):
    for x in Pi:
        Pi[x] += 1
    return Pi
    
def normalize(Pi):
    total = sum(Pi.values())
    for x in Pi:
        Pi[x] = np.float(Pi[x]) / total
    return Pi
    
def add_all_AAs_2d(Pij):
    AAs = 'GAVLIPFYWSTCMNQKRHDE'
    for a in AAs:
        for b in AAs:
            if not Pij.has_key((a,b,)):
                Pij[(a,b)] = 0
    return Pij
    