#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 14:44:36 2017

@author: mimi
"""


import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from collections import Counter
from scipy import stats
from conkit import conkit


msa_file = '/Users/mimi/Box Sync/Bioinformatics/RB helix/liban_rb.txt'

RB = AlignIO.read(msa_file, "fasta")
RB_array = np.array([list(rec) for rec in RB], np.character)

helix_motif = RB[:,965:-1]
AlignIO.write(helix_motif,'/Users/mimi/Box Sync/Bioinformatics/RB helix/helix_motif.txt','fasta')

# Extract species names
names = [rec.name for rec in RB]
names_csvfile = open('/Users/mimi/Box Sync/Bioinformatics/RB helix/names.txt','w')
print >> names_csvfile, '\n'.join(names)
names_csvfile.close()
