#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:01:53 2017

@author: mimi
"""

import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from collections import Counter
from scipy import stats

msa_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/MSA/rb_filtered_1910_1940.fasta'
rb = conkit.io.read(msa_file,'fasta')

# Read in the filtered RB .aln clustal file
msa_file = '/Users/mimi/Box Sync/Mol Bio/DNA sequences/RB family/MSA/rb_filtered.aln'
RB = AlignIO.read(msa_file, "clustal")
RB_array = np.array([list(rec) for rec in RB], np.character)

# Filtered Cyclin .aln
msa_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/Cyclin family/cyc_filtered.aln'
Cyc = AlignIO.read(msa_file, "clustal")
Cyc_array = np.array([list(rec) for rec in Cyc], np.character)

# Compute pairwise mutual information

#MI = np.zeros( (len_align,len_align) )
#for i in xrange(len_align):
#    for j in xrange(len_align):
#        if i > j:
#            MI[i,j] = mutual_information(alignments_array,i,j)
#            print "Done with ", i,j
#
## Save immediately: expensive to compute!
#np.savetxt("MI.csv", MI, delimiter=",")

# Read MI from file
#MI = np.genfromtxt('/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/MI.csv', delimiter=',')

# Visualize mutual information
MI_flat = MI.flatten()
plt.hist(MI_flat,100)

# Sort and find the highest correlated residues
I = np.argsort(MI_flat)
I = I[::-1] # Reverse order so largest is first
I = np.column_stack(np.unravel_index(I, MI.shape))

sortedMI = np.sort(MI_flat)
sortedMI = sortedMI[::-1]

st = 10
for i,score in enumerate(sortedMI[st:st+10]):
    print 'Top ', st+i+1, 'th score: ', score, ' at alignment position: ', I[i], '\n'

# Save RB motif region
motif_beg = 1910
motif_end = 1940

motif_aligned = rb_cc[:,motif_beg:motif_end]
AlignIO.write([motif_aligned],
              '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/MSA/rb_cc_1910_1940.fasta','fasta')

# Save Cyc motif only
motif_beg = 1500
motif_end = 1590

motif_aligned = cyc_cc[:,motif_beg:motif_end]
AlignIO.write([motif_aligned],
              '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/Cyclin family/cyclin_cc_1500_1590.fasta','fasta')

