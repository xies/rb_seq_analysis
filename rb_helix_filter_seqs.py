#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 17:45:43 2018

@author: mimi
"""

import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from numpy import genfromtxt
import re


# Load trimmed FASTA files
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed.mfa'
hmmer3_mafft_trimmed = AlignIO.read(filename,'fasta')
# Fix the record names
for rec in hmmer3_mafft_trimmed:
    name = rec.name
    # Negative lookbehind regex
    new_name = re.search('(?<!\/)\w+', name)
    rec.name = new_name.group(0)



##--- filter by occupancy ---
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/trimmed_occupancy.txt'
occupancy = genfromtxt(filename, delimiter=',')

filtered = []
I = occupancy > 35 #threshold = at least 7 sequences
for rec in hmmer3_mafft_trimmed:
    s = np.array(rec.seq)[I]
    rec.seq = Seq(''.join(s))
    filtered.append(rec)

SeqIO.write(filtered,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_occupancy>35.mfa','fasta')


##--- Filter by phyllum ---
# metazoa
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/metazoa/metazoa.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
metazoa = filter_by_gene_name(filtered,names)
SeqIO.write(metazoa,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_metazoa.fasta','fasta')

# viridiplantae
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/viridiplantae/viridiplantae.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
viridiplantae = filter_by_gene_name(filtered,names)
SeqIO.write(viridiplantae,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/viridiplantae/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_viridiplantae.fasta','fasta')

# oomycetes
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/oomycetes/oomycetes.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
oomycetes = filter_by_gene_name(filtered,names)
SeqIO.write(oomycetes,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oomycetes/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oomycetes.fasta','fasta')



def filter_by_gene_name(seqs,target):
    # Filter an AlignIO object by the records that are found within list of record names   
    filtered = []
    to_be_filtered = np.array([rec.name for rec in seqs])
    for n in target:
        # Find oomycete seq name and pull up the corresponding
        hits = find(n == to_be_filtered)
        if sum(hits.shape) > 0:        
            filtered.append(seqs[hits[0]])
            
    return filtered
        
