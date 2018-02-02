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

# Filter specific region
cycD_mrk_helix = cycD_alignment[:,436:464]
AlignIO.write(cycD_mrk_helix,''.join((path.dirname(cyc_aln),'/cycd_aln_mrk_helix.fasta')),'fasta')

