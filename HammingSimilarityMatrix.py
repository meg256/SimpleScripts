#takes in a multifasta with seqs of equal length, where multifasta is saved in working directory.
#command should be in format HammingSimilarityMatrix.py filename.fasta
#Returns similarity matrix in .txt format

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

#make dict of each fasta in form of [header]=seq
fastadict={}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
    sequence = sequence.replace("\n","") # Take line breaks out of sequence.
    fastadict[header]=sequence

def hamming_distance(s1, s2) -> int:
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

#make dataframe with key names as index
df = pd.DataFrame(columns=fastadict.keys(),index=fastadict.keys())

#iterate through df and put distance metrics into cells
for i,r in df.iterrows():
    for col in df.columns:
        query=fastadict[r.name]
        seqlength=len(query)
        ref=fastadict[col]
        hamming = hamming_distance(query,ref)
        hamdist=(seqlength-hamming)/seqlength #proportion of same-character positions to total sequence length
        df.at[r.name,col] = hamdist

#print to tsv
df.to_csv('matrix.txt', index=True)
