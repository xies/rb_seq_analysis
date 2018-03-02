#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 20:16:13 2018

@author: mimi
"""

import numpy as np
from Bio import SeqIO, AlignIO, Phylo
from collections import Counter
from scipy import stats
from os import path
import re

msa_file = '/Users/mimi/Box Sync/Bioinformatics/RB helix/MSA/metazoa.aln'
RB = AlignIO.read(msa_file, "clustal")
RB_array = np.array([list(rec) for rec in RB], np.character)

# Filter
RB_filt = RB[:,1300:]
AlignIO.write(RB_filt,''.join((path.dirname(msa_file),'/metazoa_aln_filtered.fasta')),'fasta')


#-------

# Load RB metazoan helix file
rb_helix = '/Users/mimi/Box Sync/Bioinformatics/RB helix/1__metazoa_phmmer_profile/metazoa_aln_filtered.fasta'
RB_helix = AlignIO.read(rb_helix,'fasta')

# Load cyclinD for mutual info analysis
# Filtered Cyclin fasta
cyc_file = '/Users/mimi/Box Sync/Bioinformatics/RB helix/cyclinD_mutual_info/cyclins.fasta'
# Filter for D-type cyclins and write to new file
cycD = [rec for rec in SeqIO.parse(cyc_file,'fasta') if rec.id.find('CycD') > 0]
SeqIO.write(cycD, ''.join((path.dirname(cyc_file),'/cyclinD.fasta')),'fasta')

# Load Cyclin .aln
cyc_aln =  '/Users/mimi/Box Sync/Bioinformatics/RB helix/cyclinD_mutual_info/cyclind.aln'
cycD_alignment = AlignIO.read(cyc_aln,'clustal')
AlignIO.write(cycD_alignment,''.join((path.dirname(cyc_aln),'/cycd_aln.fasta')),'fasta')

# Filter specific region
cycD_mrk_helix = cycD_alignment[:,436:464]
AlignIO.write(cycD_mrk_helix,''.join((path.dirname(cyc_aln),'/cycd_aln_mrk_helix.fasta')),'fasta')


# --- 
# Load Liban file
rb_liban = '/Users/mimi/Box Sync/Bioinformatics/RB helix/liban_rb.txt'
RB_liban = AlignIO.read(rb_liban,'fasta')
RB_liban_filtered = []
for rec in RB_liban[:,965:]:
    if ~all(unique(rec.seq) == '-'):
        RB_liban_filtered.append(rec)
#        SeqIO.write(rec,''.join((path.dirname(rb_liban),'/liban_trimmed_helix.fasta')),'fasta')
        
SeqIO.write(RB_liban_filtered,''.join((path.dirname(rb_liban),'/liban_trimmed_helix.fasta')),'fasta')

#Elife_RB/72_eukaryotes
# --- Clean up hmmer3 server outputs ----

filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/3__mafft_align/hmmer3_e-10_input_clustalo_mafft_server.fasta'
hmmer3_mafft = AlignIO.read(filename,'fasta')

# Find and print duplicate names for manual deletion
names = np.array([rec.name for rec in hmmer3_mafft])
all_aligned = []
for rec in hmmer3_mafft:
    if rec.name in names:
        print 'Duplicate ', rec.name

# Load trimmed FASTA files
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed.mfa'
hmmer3_mafft_trimmed = AlignIO.read(filename,'fasta')

# Fix the record names
for rec in hmmer3_mafft_trimmed:
    name = rec.name
    # Negative lookbehind regex
    new_name = re.search('(?<!\/)\w+', name)
    rec.name = new_name.group(0)

    
#
## Load metazoa
#filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/2__hmmsearch/metazoa/hmmer3_e-10_input_clustalo_metazoa.fasta'
#metazoa = []
#names = np.array([rec.name for rec in hmmer3_mafft_trimmed])
#for rec in SeqIO.parse(filename,'fasta'):
#    # Find metazoan seq name and pull up the corresponding
#    hits = find(rec.name == names)
#    if sum(hits.shape) > 0:        
#        metazoa.append(hmmer3_mafft_trimmed[find(rec.name == names)[0]])
#    
#SeqIO.write(metazoa,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/3__mafft_align/metazoa/hmmer3_e-10_input_clustalo_mafft_server_trimmed_metazoa.fasta','fasta')
#
#
## Load viridiplantae
#filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/2__hmmsearch/viridiplantae/viridiplantae.fa'
#plants = []
#names = np.array([rec.name for rec in hmmer3_mafft_trimmed])
#for rec in SeqIO.parse(filename,'fasta'):
#    # Find plants seq name and pull up the corresponding
#    hits = find(rec.name == names)
#    if sum(hits.shape) > 0:        
#        plants.append(hmmer3_mafft_trimmed[find(rec.name == names)[0]])
#    
#SeqIO.write(plants,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/3__mafft_align/viridiplantae/hmmer3_e-10_input_clustalo_mafft_server_trimmed_plants.fasta','fasta')
#
#
## Load oomycetes
#filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/2__hmmsearch/oomycetes/oomycetes.fa'
#oomycetes = []
#names = np.array([rec.name for rec in hmmer3_mafft_trimmed])
#for rec in SeqIO.parse(filename,'fasta'):
#    # Find oomycete seq name and pull up the corresponding
#    hits = find(rec.name == names)
#    if sum(hits.shape) > 0:        
#        oomycetes.append(hmmer3_mafft_trimmed[find(rec.name == names)[0]])
#    
#SeqIO.write(oomycetes,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/3__mafft_align/oomycetes/hmmer3_e-10_input_clustalo_mafft_server_trimmed_oomycetes.fasta','fasta')
#
#
#
## Load heterotrichea
#filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/2__hmmsearch/heterotrichea/heterotrichea.fasta'
#heterotrichea = []
#names = np.array([rec.name for rec in hmmer3_mafft_trimmed])
#for rec in SeqIO.parse(filename,'fasta'):
#    # Find oomycete seq name and pull up the corresponding
#    hits = find(rec.name == names)
#    if sum(hits.shape) > 0:        
#        heterotrichea.append(hmmer3_mafft_trimmed[find(rec.name == names)[0]])
#    
#SeqIO.write(heterotrichea,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/3__mafft_align/heterotrichea/hmmer3_e-10_input_clustalo_mafft_server_trimmed_heterotrichea.fasta','fasta')
#




    
    