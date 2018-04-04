#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 19:44:25 2018

@author: mimi
"""


import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from numpy import genfromtxt
import re
from weblogolib import *



opt = LogoOptions()
opt.logo_start = 1
opt.logo_end = 19

# Metazoa
opt.logo_title = 'Metazoan RB family'
write_logo('/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_names_trimmed_humanRB.fasta',
           '/Users/mimi/Box Sync/Bioinformatics/RB helix/Elife_RB/72_eukaryotes/mafft/4__mafft_align/metazoa/readable_names_trimmed_humanRB.eps',
           opt)





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
    
