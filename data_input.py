#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:01:53 2017

@author: mimi
"""

from Bio import SeqIO

fasta_sequences = SeqIO.parse(open(input_file),'fasta')

with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        new_sequence = some_function(sequence)
        write_fasta(out_file)