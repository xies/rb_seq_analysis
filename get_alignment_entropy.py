#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:01:53 2017

@author: mimi
"""

from Bio import SeqIO, AlignIO, Phylo

fasta_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.fasta'

fasta = SeqIO.parse(open(fasta_file),'fasta')
# Read in FASTA files
seqs = {}
for f in fasta:
    seqs[f.id] = f.seq

# Read in the .aln clustal file
msa_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.aln'
alignments = AlignIO.read(msa_file, "clustal")
alignment_array = np.array([list(rec) for rec in alignments], np.character)


# Read in the .dnd phylogenetic tree file
dnd_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.dnd'
tree = Phylo.read(dnd_file,'newick')

