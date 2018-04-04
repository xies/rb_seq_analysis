#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 20:16:13 2018

@author: mimi
"""

import numpy as np
from Bio import SeqIO, AlignIO, Phylo, Seq
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
    if ~all(np.unique(rec.seq) == '-'):
        RB_liban_filtered.append(rec)
#        SeqIO.write(rec,''.join((path.dirname(rb_liban),'/liban_trimmed_helix.fasta')),'fasta')
        
SeqIO.write(RB_liban_filtered,''.join((path.dirname(rb_liban),'/liban_trimmed_helix.fasta')),'fasta')

#Elife_RB/72_eukaryotes
# --- Clean up hmmer3 server outputs ----
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_align.fasta'
hmmer3_mafft = AlignIO.read(filename,'fasta')

# Fix the record names
for rec in hmmer3_aligned:
    name = rec.name
    # Negative lookbehind regex
    new_name = re.search('(?<!\/)\w+', name)
    rec.name = new_name.group(0)

nonredundant = []
redundant = []
names_already_seen = []
for rec in hmmer3_aligned:
    if rec.name in names_already_seen:
        redundant.append(rec)
    else:
        names_already_seen.append(rec.name)
        nonredundant.append(rec)

        
##--- Load hit information from .csv file ---
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/hmmer3_e-10_input_mafft_jackhmmer.csv'
metadata = pd.read_csv(filename,delimiter='\t')

# Load trimmed FASTA files
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_trimmed.mfa'
hmmer3_mafft_trimmed = AlignIO.read(filename,'fasta')

# Load nonredundant entries
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/hmmer3_e-10_input_mafft_jackhmmer_mafft_align_nonredundant.fasta'
hmmer3_mafft = AlignIO.read(filename,'fasta')

# Load HMMER3-aligned sequence and look through it directly
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/3__jackhmmer3_aligned.fasta'
hmmer3_aligned = AlignIO.read(filename,'fasta')
# fix the . into - in gapped sequences
for rec in hmmer3_aligned:
    s = rec.seq.tostring().replace('.','-')
    rec.seq = Seq.Seq(s).upper()
SeqIO.write(hmmer3_aligned,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/3__jackhmmer3_aligned_fixed.fasta','fasta')


#----------
# Use Uniprot to figure out organism name, rewrite FASTA entries with organism + protein name + UnitProt_entrykey
u = UniProt()
new_recs = []
no_hits = []
for rec in nonredundant:
    if rec.seq.ungap('-') != '':
        d = u.quick_search(rec.name,limit=1)
        for k in d.keys():
            latin_name = re.search('(\w+)\s(\w+)',d[k]['Organism'])
            if latin_name != None:
                latin_name = '_'.join((latin_name.group(1),latin_name.group(2)))
                protein_names = '_'.join(re.findall('\w+',d[k]['Protein names']))
                entry_name = re.search('(?<!\/)\w+', rec.name).group(0)
                rec.id = '_'.join((latin_name,protein_names,entry_name))
                new_recs.append(rec)
            else:
                no_hits.append(rec)
    else:
        no_hits.append(rec)

SeqIO.write(new_recs,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/readable_names/readable_names.fasta','fasta')

SeqIO.write(nonredundant,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/non_redundant/nonredundant.fasta','fasta')



#-------

def find_redundant(seqs):
    "Return a list of redundant SeqIO .description entries"
    nonredundant = []
    redundant = []
    names_already_seen = []
    for rec in seqs:
        if rec.description in names_already_seen:
            redundant.append(rec)
        else:
            names_already_seen.append(rec.description)
            nonredundant.append(rec)
            
    return redundant,nonredundant


def split_into_individual_fastas(filename):
    #Split a .fasta file with multiple entries into multiple
    #individual .fasta files
    import os
    from Bio import SeqIO, AlignIO
     
    # Get file and directory names for unified output
    indir = os.path.dirname(filename)
    base_filename = os.path.basename(filename)
    base_noext = os.path.splitext(base_filename)[0]
    base_noext = ''.join((base_noext,'_individuals'))
    out_dirname = os.path.join(indir,base_noext)
    # make output directory if it DNE
    try:
        os.stat(out_dirname)
    except:
        os.mkdir(out_dirname)

    seqs = AlignIO.read(filename,'fasta')
    
    delete_trailing_aas_from_rec_id(seqs)
    delete_trailing_aas_from_rec_name(seqs)
    for rec in seqs:
       rec_basename = ''.join((rec.name,'.fasta'))
       rec_filename = path.join(out_dirname,rec_basename)
       rec.seq = rec.seq.ungap('-')
       if rec.seq == '':
           continue
       print rec_filename
       with open(rec_filename,'w') as output_handle:
           SeqIO.write(rec,output_handle,'fasta')
       
