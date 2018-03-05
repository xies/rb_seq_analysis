t#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 19:58:27 2017

@author: mimi
"""

import numpy as np
from Bio import SeqIO, AlignIO, Align
from collections import Counter
from scipy import stats

# Read in the RB .aln clustal file
msa_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/MSA/rb_fastas.aln'
RB = AlignIO.read(msa_file, "clustal")
rb_seqs = rec2dict(RB.get_all_seqs())

# Cyclin .aln
msa_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/Cyclin family/cyclins.aln'
Cyc = AlignIO.read(msa_file, "clustal")
cyc_seqs = rec2dict(Cyc.get_all_seqs())

#Find all Cyclin D only
cyc_seqs = {k: v for (k, v) in cyc_seqs.items() if k.find('CycD') > 0}

# Extract species name as first 4 letters
rb_species = np.array([s[0:4] for s in rb_seqs.keys()])
cyc_species = np.array([s[0:4] for s in cyc_seqs.keys()])

# Find intersection (could vectorize but fast enough)
I = np.zeros((len(rb_species))) == 1
J = np.zeros((len(cyc_species))) == 1
for (i,s) in enumerate(rb_species):
    for (j,t) in enumerate(cyc_species):
        if s == t:
            I[i] = True
            J[j] = True
            
# Create new list of SeqIO records from the filtered protein list
rb_filt = Align.MultipleSeqAlignment([])
for k in np.array(rb_seqs.keys())[I]:
    rb_filt.append(rb_seqs[k])
cyc_filt = Align.MultipleSeqAlignment([])
for k in np.array(cyc_seqs.keys())[J]:
    cyc_filt.append(cyc_seqs[k])

AlignIO.write(rb_filt,'/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/MSA/rb_filtered.aln','clustal')
AlignIO.write(cyc_filt,'/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/Cyclin family/cyc_filtered.aln','clustal')

# Rewrite FASTA file with PROTEIN | PSECIES
rb_cc = Align.MultipleSeqAlignment([])
rb_keys = np.array([s for s in np.array(rb_seqs.keys())[I]])
for k in rb_keys.argsort():
    rec = rb_seqs[ rb_keys[k] ]
    rec.id = ' '. join([ rec.id[4:], '|', rec.id[0:4]])
    rec.description=''
    rb_cc.append(rec)

cyc_cc = Align.MultipleSeqAlignment([])
cyc_keys = np.array([s for s in np.array(cyc_seqs.keys())[J]])
for k in cyc_keys.argsort():
    rec = cyc_seqs[ cyc_keys[k] ]
    rec.id = ' '.join([ rec.id[4:] ,'|', rec.id[0:4]])   
    rec.description=''
    cyc_cc.append(rec)

AlignIO.write(rb_filt,'/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/MSA/rb_cc.fasta','fasta')
AlignIO.write(cyc_filt,'/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/Cyclin family/cyc_cc.fasta','fasta')


