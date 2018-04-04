#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 17:52:13 2018

@author: mimi
"""

import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from numpy import genfromtxt
import re
from weblogolib import *
from ete3 import Tree, TreeStyle
from bioservices import UniProt

# Load occupancy-filtered FASTA files
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_names_trimmed.fasta'
metazoa = AlignIO.read(filename,'fasta')

# Fix the record names from the Jalview trimming process
for rec in metazoa:
    name = rec.name
    # Negative lookbehind regex
    new_name = re.search('(?<!\/)\w+', name)
    rec.name = new_name.group(0)
    
##--- Load hit information from .csv file ---
filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/3__jackhammer/hmmer3_e-10_input_mafft_jackhmmer.csv'
metadata = pd.read_csv(filename,delimiter='\t')
metazoa = AlignIO.read('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_names_trimmed.fasta','fasta')

# Find and print duplicate names for manual deletion
#names = np.array([rec.name for rec in metazoa])
names = []
non_redundant = []
counter = 0
for rec in metazoa:
    if rec.name in names:
        print 'Duplicate ', rec.name
        counter += 1
        rec = None
    else:
        names.append(rec.name)
        

#descriptions = []
#for rec in metazoa:
#    I = metadata['Target Name'] == rec.name
#    if sum(I) > 1:
#        J = np.argmax(metadata[I]['Target Length'])
#        descriptions.append(metadata.iloc[J,]['Description'])
#        rec.description = metadata.iloc[J,]['Description']
#    elif sum(I) == 1:
#        # To get the raw string need .values[0] to convert to np array and grab first item
#        descriptions.append(metadata[I]['Description'].values[0])
#        rec.description = metadata[I]['Description'].values[0]
#        

# USE REGEX TO FIND RB-like proteins
# Look for RBL-1 (p107) and RBL-2 (p130)
rb = []
rb_l1 = []
rbl = []
rb_l2 = []
for rec in metazoa:
    d = rec.description
    name = re.search('like.*(\d)',d)
    if name != None:
        rbl.append(rec)
#        if name.group(1) == '1':
#            rb_l1.append(rec)
#            print d
#        elif name.group(1) == '2':
#            rb_l2.append(rec)
#            print d
    else:
        name = re.search('Uncharacterized',d)
        if name == None:
            rb.append(rec)
            print d

opt = LogoOptions()
opt.logo_start = 1
opt.logo_end = 21

SeqIO.write(rb,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_trimmed_rb.mfa','fasta')
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_trimmed_rb.mfa',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_trimmed_rb.eps',
           opt)
SeqIO.write(rbl,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_trimmed_rbl.mfa','fasta')
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_trimmed_rbl.mfa',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_trimmed_rbl.eps',
           opt)


#### USING ete3 Tree visualizer ####

filename = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/metazoa.tree'
tree = Tree(filename,format=1)
tree.show()
#
## Fix the record names
#for rec in metazoa:
#    name = rec.name
#    # Negative lookbehind regex
#    new_name = re.search('(?<!\/)\w+', name)
#    rec.id = new_name

metazoa_names = [rec.name for rec in metazoa]
# Load the subtrees associated with RBL1 family and record trimmed FASTA alignment
dirname = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl1/'
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
SeqIO.write(rbl1,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl1/rbl1_tree_humanRB.mfa','fasta')
opt.logo_title = 'RB-like 1 family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl1/rbl1_tree_humanRB.mfa',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl1/rbl1_tree_humanRB.eps',
           opt)

# Load the subtrees associated with RBL2 (p130) family and record trimmed FASTA alignment
dirname = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl2/'
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
SeqIO.write(rbl2,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl2/rbl2_tree_humanRB.mfa','fasta')
opt.logo_title = 'RB-like 2 family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl2/rbl2_tree_humanRB.mfa',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rbl2/rbl2_tree_humanRB.eps',
           opt)


# Load the subtrees associated with pRB family and record trimmed FASTA alignment
dirname = '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rb/'
rb_treefiles = ['tree.tree']
rb_names = []
for filename in rb_treefiles:
    rb_names.append(get_tree_leaves(''.join((dirname,filename))))
rb_names = set([x for tree in rb_names for x in tree])
#rbl2_names = delete_leading_space(rbl2_names)
#find alignment from trimmed metazoa sequences
rb = []
for name in rb_names:
    if name in metazoa_names:
        rb.append(metazoa[metazoa_names.index(name)])
    else:
        print "Oh no!", name
#        break
SeqIO.write(rb,'/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rb/rb_tree_humanRB.mfa','fasta')
opt.logo_title = 'pRB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rb/rb_tree_humanRB.mfa',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/rb/rb_tree_humanRB.eps',
           opt)


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

#
#def delete_leading_space(list_of_string):
#    return [s.lstrip() for s in list_of_string]