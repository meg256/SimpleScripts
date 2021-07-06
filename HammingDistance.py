#takes in a multifasta with seqs of equal length.
#command should be in format hamming_distance_matrix.py filename.fasta
#Returns distance matrix in tsv format

import re
import pandas as pd
import os
import sys


#open multifasta
filepath = os.getcwd()+'//'
fasta = open(filepath+sys.argv[1])
sequences = fasta.read()
sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
fasta.close()

#first item in sequences list will be blank output
del sequences[0]

#make dict of fastas
fastadict={}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
    sequence = sequence.replace("\n","") # Take line breaks out of sequence.
    fastadict[header]=sequence

def hamming_distance(s1, s2) -> int:
    if len(s1) != len(s2):
        raise ValueError("Sequences not same length.")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

#make dataframe
df = pd.DataFrame(columns=fastadict.keys(),index=fastadict.keys())

#iterate through df and put distance metrics into cells
for i,r in df.iterrows():
    for col in df.columns:
        query=fastadict[r.name]
        seqlength=len(query)
        ref=fastadict[col]
        hamming = hamming_distance(query,ref)
        hamdist=(seqlength-hamming)/seqlength
        df.at[r.name,col] = hamdist

#print to tsv
df.to_csv('HammingMatrix.txt', index=True)
