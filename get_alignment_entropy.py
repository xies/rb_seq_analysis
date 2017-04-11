#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:01:53 2017

@author: mimi
"""

from Bio import SeqIO, AlignIO, Phylo
from collections import Counter
from scipy import stats

fasta_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.fasta'

fasta = SeqIO.parse(open(fasta_file),'fasta')
# Read in FASTA files
seqs = {}
for f in fasta:
    seqs[f.id] = f.seq

# Read in the .aln clustal file
msa_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.aln'
alignments = AlignIO.read(msa_file, "clustal")
alignments_array = np.array([list(rec) for rec in alignments], np.character)

[num_seq,len_align] = alignments_array.shape

# Read in the .dnd phylogenetic tree file
dnd_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.dnd'
tree = Phylo.read(dnd_file,'newick')

# Compute sequence frequencies
aa_frequencies = [ Counter(alignments_array[:,i]) for i in xrange(len_align) ]

# Throw out '-' and 'X'
for f in aa_frequencies:
    try:
        f.pop('-')
    except KeyError:
        pass
    
    try:
        f.pop('X')
    except KeyError:
        pass
                  
aa_entropies = [ stats.entropy(f.values()) for f in aa_frequencies ]

# RB motif
motifOI = seqs['HsapRB1'][894:914]
motif_beg = 1909
motif_end = 1929

motif_aligned = alignments[:,motif_beg:motif_end]
AlignIO.write([motif_aligned],
              '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_motif_aligned.aln','clustal')


