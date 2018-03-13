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
from weblogolib import *


# Load trimmed FASTA files
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/readable_names/readable_names.fasta'
readable = AlignIO.read(filename,'fasta')

# Load trimmed FASTA files
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/readable_names/readable_names_trimmed.fasta'
readable_trimmed = AlignIO.read(filename,'fasta')
# Fix the record names
for rec in readable_trimmed:
    name = rec.name
    # Negative lookbehind regex
    new_name = re.search('(?<!\/)\w+', name)
    new_name = '_'.join(re.split('_',new_name.group(0))[-2:])
    rec.name = new_name
#    new_id = re.search('(?<!\/)\w+',rec.id)
#    rec.id = new_id.group(0)
    rec.description = rec.name

opt = LogoOptions()
opt.logo_start = 1
opt.logo_end = 22

#--- filter by occupancy ---
#filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/readable_names/occupancy.csv'
#occupancy = genfromtxt(filename, delimiter=',')
#filtered = filter_by_occupancy(readable_trimmed,occupancy,15) # thresh = 15

# Filter by metazoan occupancy
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/metazoan_occupancy.csv'
met_occ = genfromtxt(filename, delimiter=',')
filtered = filter_by_occupancy(readable_trimmed,met_occ,20) # thresh = 20

opt.logo_title = 'Eukaryotic RB family (filtered by metazoan occupancy)'
SeqIO.write(filtered,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/readable_names/trimmed_met_occ>20.mfa','fasta')
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/readable_names/trimmed_met_occ>20.mfa',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/readable_names/trimmed_met_occ>20.eps',
           opt)

##--- Filter by phyllum ---
# metazoa
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/metazoa/metazoa.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Accession'].values
metazoa = filter_by_gene_name(readable_trimmed,names)
SeqIO.write(metazoa,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_names_trimmed.fasta','fasta')

opt.logo_title = 'Metazoan RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_names_trimmed.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_names_trimmed.eps',
           opt)

# viridiplantae
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/viridiplantae/viridiplantae.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name'].values
viridiplantae = filter_by_gene_name(filtered,names)
SeqIO.write(viridiplantae,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/viridiplantae/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_viridiplantae.fasta','fasta')

opt.logo_title = 'Viridiplantae RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/viridiplantae/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_viridiplantae.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/viridiplantae/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_viridiplantae.eps',
           opt)

# oomycetes
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/oomycetes/oomycetes.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name'].values
oomycetes = filter_by_gene_name(filtered,names)
SeqIO.write(oomycetes,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oomycetes/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oomycetes.fasta','fasta')

opt.logo_title = 'Oomycetes RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oomycetes/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oomycetes.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oomycetes/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oomycetes.eps',
           opt)

# fungi
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/fungi/fungi.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name'].values
fungi = filter_by_gene_name(filtered,names)
SeqIO.write(fungi,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/fungi/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_fungi.fasta','fasta')

opt.logo_title = 'Fungi RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/fungi/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_fungi.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/fungi/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_fungi.eps',
           opt)

# oligohymenophorea
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/oligohymenophorea/oligohymenophorea.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name'].values
oligohymenophorea = filter_by_gene_name(filtered,names)
SeqIO.write(oligohymenophorea,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oligohymenophorea/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oligohymenophorea.fasta','fasta')

opt.logo_title = 'Oligohymenophorea RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oligohymenophorea/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oligohymenophorea.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/oligohymenophorea/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_oligohymenophorea.eps',
           opt)

# bacillariophyta
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/bacillariophyta/bacillariophyta.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name'].values
bacillariophyta = filter_by_gene_name(filtered,names)
SeqIO.write(bacillariophyta,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/bacillariophyta/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_bacillariophyta.fasta','fasta')

opt.logo_title = 'Bacillariophyta RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/bacillariophyta/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_bacillariophyta.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/bacillariophyta/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_bacillariophyta.eps',
           opt)


# choanoflagellida
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/choanoflagellida/choanoflagellida.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name'].values
choanoflagellida = filter_by_gene_name(filtered,names)
SeqIO.write(bacillariophyta,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanoflagellida/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_choanoflagellida.fasta','fasta')

opt.logo_title = 'Choanoflagellida RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanoflagellida/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_choanoflagellida.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanoflagellida/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_choanoflagellida.eps',
           opt)

# dictyosteliales
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/dictyosteliales/dictyosteliales.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Name'].values
dictyosteliales = filter_by_gene_name(filtered,names)
SeqIO.write(dictyosteliales,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/dictyosteliales/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_dictyosteliales.fasta','fasta')

opt.logo_title = 'Dictyosteliales RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/dictyosteliales/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_dictyosteliales.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/dictyosteliales/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_dictyosteliales.eps',
           opt)

# SAR supergroup
sar = oomycetes + oligohymenophorea + bacillariophyta + choanoflagellida
SeqIO.write(sar,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/sar/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_sar.fasta','fasta')
opt.logo_title = 'SAR RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/sar/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_sar.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/sar/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed_sar.eps',
           opt)

# LOAD manually curated fasta sets

#zoosporic (basal fungi)
zoosporic = AlignIO.read('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/zoosporic/zoosporia.fasta','fasta')
zoosporic = filter_by_occupancy(zoosporic,met_occ,20)
SeqIO.write(zoosporic,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/zoosporic/zoosporia_met_occ>20.fasta','fasta')
opt.logo_title = 'Zoosporia (basal fungi) RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/zoosporic/zoosporia_met_occ>20.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/zoosporic/zoosporia_met_occ>20.eps',
           opt)

#apusazoa
apusozoa = AlignIO.read('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/apusozoa/apusozoa.fasta','fasta')
apusozoa = filter_by_occupancy(apusozoa,met_occ,20)
SeqIO.write(apusozoa,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/apusozoa/apusozoa_met_occ>20.fasta','fasta')
opt.logo_title = 'Apusozoa RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/apusozoa/apusozoa_met_occ>20.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/apusozoa/apusozoa_met_occ>20.eps',
           opt)

#amoebazoa (basal fungi)
amoebozoa = AlignIO.read('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/amoebozoa/amoebozoa.fasta','fasta')
amoebozoa = filter_by_occupancy(amoebozoa,met_occ,20)
SeqIO.write(amoebozoa,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/amoebozoa/amoebozoa_met_occ>20.fasta','fasta')
opt.logo_title = 'Amoebozoa (slime molds) RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/amoebozoa/amoebozoa_met_occ>20.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/amoebozoa/amoebozoa_met_occ>20.eps',
           opt)

#choanomonada
choanomonada = AlignIO.read('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanomonada/choanomonada.fasta','fasta')
#choanomonada = filter_by_occupancy(choanomonada,met_occ,20)
SeqIO.write(choanomonada,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanomonada/choanomonada_met_occ>20.fasta','fasta')
opt.logo_title = 'Amoebozoa (slime molds) RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanomonada/choanomonada_met_occ>20.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/choanomonada/choanomonada_met_occ>20.eps',
           opt)

def filter_by_occupancy(seqs,occupancy,thresh):
    """ Filter by an occupancy array according to provided threshold; will also
    filter out completely empty sequences
    """
    I = occupancy > thresh
    filtered =[]
    for rec in seqs:
        s = np.array(rec.seq)[I]
        rec.seq = Seq(''.join(s))
        if rec.seq.ungap('-') != '':
            filtered.append(rec)
    return filtered
    

def filter_by_gene_name(seqs,target_list):
    # Filter an AlignIO object by the records that are found within list of record names   
    filtered = []
    for rec in seqs:
        # Find if seq name is within target_list
       if rec.name in target_list:
           filtered.append(rec)
    return filtered
        

def write_logo(filename,out_name,options):
    # Wrapper for weblogo
    fin = open(filename)
    seqs = read_seq_data(fin)
    data = LogoData.from_seqs(seqs)
    form = LogoFormat(data,options)
    eps = eps_formatter(data,form)
    fout = open(out_name,'w')
    fout.write(eps)
    fout.close()
    fin.close()
    
