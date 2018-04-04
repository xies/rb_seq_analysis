#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:45:10 2018

@author: xies
"""
import numpy as np
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from numpy import genfromtxt
import re, os
from weblogolib import *

BASE_DIR = '/home/xies/Desktop/Elife_RB/72_eukaryotes/mafft/4__mafft_align/'

# Load FASTA files
filename = os.path.join(BASE_DIR,'readable_names/readable_names.fasta')
readable = AlignIO.read(filename,'fasta')

# Trim 
left_bound = 9237
right_bound = 9655
for rec in readable:
    rec.seq = rec.seq[left_bound:right_bound]
filename = os.path.join(BASE_DIR,'readable_names/readable_names_trimmed_python.fasta')
SeqIO.write(readable,filename,'fasta')

# Load trimmed FASTA files
filename = os.path.join(BASE_DIR,'readable_names/readable_names_trimmed_python.fasta')
readable = AlignIO.read(filename,'fasta')

# Fix the record names
for rec in readable:
    name = rec.id
    # Negative lookbehind regex
    new_name = re.search('(?<!\/)\w+', name)
    new_name = '_'.join(re.split('_',new_name.group(0))[-2:])
    rec.name = new_name
#    new_id = re.search('(?<!\/)\w+',rec.id)
#    rec.id = new_id.group(0)

# 2----- Load metazoa phylla info
# a) Filter from overall eukaryotic info
filename = '/home/xies/Desktop/Elife_RB/72_eukaryotes/mafft/3__jackhammer/metazoa/metazoa.csv'
df = pd.read_csv(filename,delimiter='\t')
names = df['Target Accession'].values
[Imet,metazoa] = get_phylla_membership(readable,names,'metazoa')
# b) OR: manually load .mfa file
metazoa = AlignIO.read(os.path.join(BASE_DIR,'metazoa/manually_filtered/mafft.fasta'),'fasta')
delete_trailing_aas_from_rec_name(metazoa)
delete_trailing_aas_from_rec_id(metazoa)

# c) Get metazoan occupancy
occMet = get_occupancy(metazoa)
threshold = 50
readable_metocc = filter_by_occupancy(readable,occMet,threshold)
filename = os.path.join(BASE_DIR,''.join(('readable_names/readable_names_metocc>',str(threshold),'.fasta')))
SeqIO.write(readable_metocc,filename,'fasta')

# d) generate phylla specific logos ------
filename = os.path.join(BASE_DIR,''.join(('metazoa/readable_names_metocc>',str(threshold),'.fasta')))
SeqIO.write(metazoa,filename,'fasta')


opt = LogoOptions()
opt.logo_start = 22
opt.logo_end = 42

opt.logo_title = 'Manually filtered -- metazoan pocket proteins'
write_logo(os.path.join(BASE_DIR,'metazoa/manually_filtered/mafft_trimmed'),
           opt)

# ------ Generate subfamily logos -----
metazoa_names = [rec.name for rec in metazoa]

# Load the subtrees associated with pRB family and record trimmed FASTA alignment
rb_treefiles = ['metazoa/rb/tree.tree']
rb_names = []
for filename in rb_treefiles:
    rb_names.append(get_tree_leaves(''.join((BASE_DIR,filename))))
rb_names = set([x for tree in rb_names for x in tree])
#find alignment from trimmed metazoa sequences
rb = []
for name in rb_names:
    if name in metazoa_names:
        rb.append(metazoa[metazoa_names.index(name)])
    else:
        print "Oh no!", name
#        break
print "There are ", str(len(rb)), " pRB subfamily sequences"
SeqIO.write(rb,os.path.join(BASE_DIR,'metazoa/manually_filtered/rb/rb.fasta'),'fasta')
opt.logo_title = 'pRB subfamily'
write_logo(os.path.join(BASE_DIR,'metazoa/manually_filtered/rb/rb'),
           opt)

# Load the subtrees associated with RBL1 family and record trimmed FASTA alignment
dirname = '/home/xies/Desktop/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl1/'
rbl1_treefiles = ['tree.tree','tree2.tree','tree3.tree','tree4.tree','tree5.tree']
rbl1_names = []
for filename in rbl1_treefiles:
    rbl1_names.append(get_tree_leaves(''.join((dirname,filename))))
rbl1_names = set([x for tree in rbl1_names for x in tree])
#find alignment from trimmed metazoa sequences
rbl1 = []
for name in rbl1_names:
    if name in metazoa_names:
        rbl1.append(metazoa[metazoa_names.index(name)])
    else:
        print "Oh no!", name
#        break
print "There are ", str(len(rbl1)), " RB-like 1 subfamily sequences"
SeqIO.write(rbl1,os.path.join(BASE_DIR,'metazoa/manually_filtered/rbl1/rbl1.fasta'),'fasta')
opt.logo_title = 'RB-like 1 subfamily'
write_logo(os.path.join(BASE_DIR,'metazoa/manually_filtered/rbl1/rbl1'),
           opt)

# Load the subtrees associated with RBL2 (p130) family and record trimmed FASTA alignment
dirname = '/home/xies/Desktop/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl2/'
rbl2_treefiles = ['tree.tree','tree2.tree','tree3.tree','tree4.tree','tree5.tree']
rbl2_names = []
for filename in rbl2_treefiles:
    rbl2_names.append(get_tree_leaves(''.join((dirname,filename))))
rbl2_names = set([x for tree in rbl2_names for x in tree])
#rbl2_names = delete_leading_space(rbl2_names)
#find alignment from trimmed metazoa sequences
rbl2 = []
for name in rbl2_names:
    if name in metazoa_names:
        rbl2.append(metazoa[metazoa_names.index(name)])
    else:
        print "Oh no!", name
#        break
print "There are ", str(len(rbl2)), " RB-like 2 subfamily sequences"
SeqIO.write(rbl1,os.path.join(BASE_DIR,'metazoa/manually_filtered/rbl2/rbl2.fasta'),'fasta')
opt.logo_title = 'RB-like 1 subfamily'
write_logo(os.path.join(BASE_DIR,'metazoa/manually_filtered/rbl2/rbl2'),
           opt)



#------------------------delete_trailing_aas_from_rec_id-----------------


def get_occupancy(seqs):
    # Get occupancy
    Npos = len(seqs[0].seq)
    occupancy = np.zeros(Npos)
    for rec in seqs:
        s = np.array(rec.seq)
        I = s != '-'
        occupancy += I        
    return occupancy

def get_phylla_membership(seqs,target_list,phylla_name):
    # Filter an AlignIO object by the records that are found within list of record names   
    I = np.zeros(len(seqs))
# Load the subtrees associated with pRB family and record trimmed FASTA alignment
rb_treefiles = ['metazoa/rb/tree.tree']
rb_names = []
for filename in rb_treefiles:
    rb_names.append(get_tree_leaves(''.join((BASE_DIR,filename))))
rb_names = set([x for tree in rb_names for x in tree])
#find alignment from trimmed metazoa sequences
rb = []
for name in rb_names:
    if name in metazoa_names:
        rb.append(metazoa[metazoa_names.index(name)])
    else:
        print "Oh no!", name
#        break
    filtered = []
    for (i,rec) in enumerate(seqs):
       # Find if seq name is within target_list
       if rec.name in target_list:
           I[i] = 1
           filtered.append(rec)
    return I,filtered

def finter_by_single_s_occupancy(seqs,gene):
    filtered = []
    I = np.array( gene.seq ) != '-'
    for rec in seqs:
        s = np.array(rec.delete_trailing_aas_from_rec_idseq)[I]
        rec.seq = Seq(''.join(s))
        if rec.seq.ungap('-') != '':
            filtered.append(rec)
    return filtered

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
        

def write_logo(base_name,options):
    # Wrapper for weblogo
    filename = ''.join((base_name,'.fasta'))
    out_name = ''.join((base_name,'.p'))
    fin = open(filename)
    seqs = read_seq_data(fin)
    data = LogoData.from_seqs(seqs)
    form = LogoFormat(data,options)
    eps = eps_formatter(data,form)
    fout = open(out_name,'w')
    fout.write(eps)
    fout.close()
    fin.close()
    

def get_tree_leaves(tree_filename):
    # Read Newick subtrees and trim away everything except leaf names
    fin = open(tree_filename); txt = fin.read(); fin.close();
    txt = re.split('[:\(\),\n]+',txt)
    names = []
    for x in txt:
        # Search for any alphabet characters
        if re.search('[a-z]',x) != None:
            names.append(x.lstrip())
    return names

def delete_trailing_aas_from_rec_name(seqs):
    # Fix the record names
    for rec in seqs:
        name = rec.name
        # Negative lookbehind regex
        new_name = re.search('(?<!\/)\w+', name)
#        new_name = '_'.join(re.split('_',new_name.group(0))[-2:])
        rec.name = new_name.group(0)
        

def delete_trailing_aas_from_rec_id(seqs):
    # Fix the record id to delete anything after a /
    for rec in seqs:
        name = rec.id
        # Negative lookbehind regex
        new_name = re.search('(?<!\/)\w+', name)
        new_name = '_'.join(re.split('_',new_name.group(0))[-2:])
        rec.id = new_name
        


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
       
