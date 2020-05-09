#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 14:42:50 2020

@author: xies
"""

from os import path
from skbio import Protein, TabularMSA


#%% Load initial alignment used to seed JackHMMER
 
basedir = '/Users/xies/Box/Bioinformatics/Cyclin family/'

filename = '1__medina_elife_cycd_aln.fasta'

cycD_init = TabularMSA.read(path.join(basedir,filename),constructor=Protein)
# Reindex using FASTA IDs
seqnames = [seq.metadata['id'] for seq in cycD_init]
cycD_init.index = seqnames

#%% Load jackHMMER sequences

# EMBL-EBI server
# hmmsearch
# -E 1e-50 --domE 1e-50 --incE 1e-50 --incdomE 1e-50 --seqdb uniprotkb

#  9DB88C1A-918F-11EA-ABBF-2463F75AEC3D.1

filename = '3__cycd_sb2.fasta'
sb2 = TabularMSA.read(path.join(basedir,filename),constructor=Protein)
