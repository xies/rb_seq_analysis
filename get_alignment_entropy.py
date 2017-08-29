#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 16:01:53 2017

@author: mimi
"""

from Bio import SeqIO, AlignIO, Phylo
from collections import Counter
from scipy import stats

fasta_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.fasta'

fasta = SeqIO.parse(open(fasta_file),'fasta')
# Read in FASTA files
seqs = {}
for f in fasta:
    seqs[f.id] = f.seq

# Read in the .aln clustal file
msa_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.aln'
alignments = AlignIO.read(msa_file, "clustal")
alignments_array = np.array([list(rec) for rec in alignments], np.character)

[num_seq,len_align] = alignments_array.shape

# Read in the .dnd phylogenetic tree file
dnd_file = '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_fastas.dnd'
tree = Phylo.read(dnd_file,'newick')

# Compute sequence frequencies
aa_frequencies = [ Counter(alignments_array[:,i]) for i in xrange(len_align) ]

# Throw out '-' and 'X'
for f in aa_frequencies:
    try:
        f.pop('-')
    except KeyError:
        pass
    
    try:
        f.pop('X')
    except KeyError:
        pass
                  
aa_entropies = [ stats.entropy(f.values()) for f in aa_frequencies ]

# Compute pairwise mutual information

#MI = np.zeros( (len_align,len_align) )
#for i in xrange(len_align):
#    for j in xrange(len_align):
#        if i > j:
#            MI[i,j] = mutual_information(alignments_array,i,j)
#            print "Done with ", i,j
#
## Save immediately: expensive to compute!
#np.savetxt("MI.csv", MI, delimiter=",")

# Read MI from file
MI = np.genfromtxt('/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/MI.csv', delimiter=',')

# Visualize mutual information
MI_flat = MI.flatten()
plt.hist(MI_flat,100)

# Sort and find the highest correlated residues
I = np.argsort(MI_flat)
I = I[::-1] # Reverse order so largest is first
I = np.column_stack(np.unravel_index(I, MI.shape))

sortedMI = np.sort(MI_flat)
sortedMI = sortedMI[::-1]

for i,score in enumerate(sortedMI[0:10]):
    print 'Top ', i+1, 'th score: ', score, ' at alignment position: ', I[i], '\n'

# RB motif
motifOI = seqs['HsapRB1'][894:914]
motif_beg = 1909
motif_end = 1929

motif_aligned = alignments[:,motif_beg:motif_end]
AlignIO.write([motif_aligned],
              '/Users/mimi/Dropbox (Personal)/Mol Bio/DNA sequences/RB family/rb_motif_aligned.aln','clustal')



####


def mutual_information(seqs,i,j):
     
    seqs = [s for s in seqs if not '-' in [s[i],s[j]]]
    # Get position frequency table, pseudocount (+1), and normalize
    Pi = Counter(s[i] for s in seqs)
    Pi = add_all_AAs(Pi)
    Pi = add_pseudocount(Pi)
    Pi = normalize(Pi)
    
    # Other position
    Pj = Counter(s[j] for s in seqs)
    Pj = add_all_AAs(Pj)
    Pj = add_pseudocount(Pj)
    Pj = normalize(Pj)
    
    # Joint probability
    Pij = Counter((s[i],s[j]) for s in seqs)
    Pij = add_all_AAs_2d(Pij)
    Pij = add_pseudocount(Pij)
    Pij = normalize(Pij)
    
    return sum(Pij[(x,y)] * \
                  log(  np.float(Pij[(x,y)]) / np.float(Pi[x]*Pj[y]) ) \
                  for x,y in Pij)

def add_all_AAs(Pi):
    AAs = 'GAVLIPFYWSTCMNQKRHDE'
    for a in AAs:
        if not Pi.has_key(a):
            Pi[a] = 0
    return Pi
    
def add_pseudocount(Pi):
    for x in Pi:
        Pi[x] += 1
    return Pi
    
def normalize(Pi):
    total = sum(Pi.values())
    for x in Pi:
        Pi[x] = np.float(Pi[x]) / total
    return Pi
    
def add_all_AAs_2d(Pij):
    AAs = 'GAVLIPFYWSTCMNQKRHDE'
    for a in AAs:
        for b in AAs:
            if not Pij.has_key((a,b,)):
                Pij[(a,b)] = 0
    return Pij
    