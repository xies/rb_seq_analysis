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
import weblogolib
from weblogolib import *


# Load trimmed FASTA files
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed.mfa'
hmmer3_mafft_trimmed = AlignIO.read(filename,'fasta')
# Fix the record names
for rec in hmmer3_mafft_trimmed:
    name = rec.name
    # Negative lookbehind regex
    new_name = re.search('(?<!\/)\w+', name)
    rec.name = new_name.group(0)

opt = LogoOptions()
opt.logo_start = 9
opt.logo_end = 35

##--- filter by occupancy ---
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/trimmed_occupancy.txt'
occupancy = genfromtxt(filename, delimiter=',')

filtered = []
I = occupancy > 40 #threshold = at least 7 sequences
for rec in hmmer3_mafft_trimmed:
    s = np.array(rec.seq)[I]
    rec.seq = Seq(''.join(s))
    filtered.append(rec)

opt.logo_title = 'Eukaryotic RB family'
SeqIO.write(filtered,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_occupancy>40.mfa','fasta')
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_occupancy>40.mfa',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_occupancy>40.eps',
           opt)

##--- Filter by phyllum ---
# metazoa
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/metazoa/metazoa.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
metazoa = filter_by_gene_name(filtered,names)
SeqIO.write(metazoa,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_metazoa.fasta','fasta')

opt.logo_title = 'Metazoan RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_metazoa.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_metazoa.eps',
           opt)

# viridiplantae
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/viridiplantae/viridiplantae.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
viridiplantae = filter_by_gene_name(filtered,names)
SeqIO.write(viridiplantae,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/viridiplantae/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_viridiplantae.fasta','fasta')

opt.logo_title = 'Viridiplantae RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/viridiplantae/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_viridiplantae.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/viridiplantae/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_viridiplantae.eps',
           opt)


# oomycetes
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/oomycetes/oomycetes.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
oomycetes = filter_by_gene_name(filtered,names)
SeqIO.write(oomycetes,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oomycetes/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oomycetes.fasta','fasta')

opt.logo_title = 'Oomycetes RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oomycetes/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oomycetes.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oomycetes/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oomycetes.eps',
           opt)

# fungi
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/fungi/fungi.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
fungi = filter_by_gene_name(filtered,names)
SeqIO.write(fungi,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/fungi/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_fungi.fasta','fasta')

opt.logo_title = 'Fungi RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/fungi/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_fungi.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/fungi/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_fungi.eps',
           opt)

# oligohymenophorea
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/oligohymenophorea/oligohymenophorea.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
oligohymenophorea = filter_by_gene_name(filtered,names)
SeqIO.write(oligohymenophorea,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oligohymenophorea/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oligohymenophorea.fasta','fasta')

opt.logo_title = 'Oligohymenophorea RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oligohymenophorea/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oligohymenophorea.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oligohymenophorea/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oligohymenophorea.eps',
           opt)

# bacillariophyta
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/bacillariophyta/bacillariophyta.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
bacillariophyta = filter_by_gene_name(filtered,names)
SeqIO.write(bacillariophyta,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/bacillariophyta/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_bacillariophyta.fasta','fasta')

opt.logo_title = 'Bacillariophyta RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/bacillariophyta/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_bacillariophyta.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/bacillariophyta/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_bacillariophyta.eps',
           opt)


# choanoflagellida
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/choanoflagellida/choanoflagellida.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
choanoflagellida = filter_by_gene_name(filtered,names)
SeqIO.write(bacillariophyta,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanoflagellida/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_choanoflagellida.fasta','fasta')

opt.logo_title = 'Choanoflagellida RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanoflagellida/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_choanoflagellida.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanoflagellida/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_choanoflagellida.eps',
           opt)

# dictyosteliales
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/dictyosteliales/dictyosteliales.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name']
dictyosteliales = filter_by_gene_name(filtered,names)
SeqIO.write(dictyosteliales,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/dictyosteliales/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_dictyosteliales.fasta','fasta')

opt.logo_title = 'Dictyosteliales RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/dictyosteliales/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_dictyosteliales.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/dictyosteliales/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_dictyosteliales.eps',
           opt)


# SAR supergroup
sar = oomycetes + oligohymenophorea + bacillariophyta + choanoflagellida
SeqIO.write(sar,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/sar/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_sar.fasta','fasta')
opt.logo_title = 'SAR RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/sar/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_sar.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/sar/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_sar.eps',
           opt)



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
        

def write_logo(filename,out_name,options):
    fin = open(filename)
    seqs = read_seq_data(fin)
    data = LogoData.from_seqs(seqs)
    form = LogoFormat(data,options)
    eps = eps_formatter(data,form)
    fout = open(out_name,'w')
    fout.write(eps)
    fout.close()
    fin.close()
    