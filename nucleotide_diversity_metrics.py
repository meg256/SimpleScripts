#This script takes in a multifasta, aligned and trimmed to equal lengths, and
#returns both pi and also the average nucleotide diversity per 100 nts.


import re
import os
import sys

#READ IN INPUT FASTA
filepath = os.getcwd()+'//'
fasta = open(filepath+sys.argv[1])
sequences = fasta.read()
sequences = re.split("^>", sequences, flags=re.MULTILINE) # Only splits string at the start of a line.
fasta.close()
#the first in this list will be blank, so delete it
del sequences[0]
#read in dict of header:seq
fastadict={}
for fasta in sequences:
    header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
    header = ">" + header + "\n" # Replace ">" lost in ">" split, Replace "\n" lost in split directly above.
    sequence = sequence.replace("\n","") # Take line breaks out of sequence.
    fastadict[header]=sequence

#RETURN ERROR IF SEQS ARE DIFFERENT LENGTHS
addlist=0
for el in fastadict.keys():
    seqlen=len(fastadict[el])
    addlist=addlist+seqlen
avglen= addlist/len(fastadict.keys())
boollist=[]
for el in fastadict.keys():
    if len(fastadict[el])==avglen:
        boollist.append("True")
    else:
        boollist.append("False")
if "False" in boollist:
    raise ValueError("Error: seqs are not all same length.")

#function to calculate total number of different characters between two seqs (assumed same length)
def hamming_distance(s1, s2) -> int:
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

#CALCULATE HAMMING DISTANCE FOR EVERY SEQUENCE PAIR
hamminglist=[]
donelist=[]
#iterate through dict and calculate all nonredundant pairwise comparisons
for key in fastadict.keys():
    for el in fastadict.keys():
        if key != el and (key+el) not in donelist and (el+key) not in donelist:
            query=fastadict[key]
            ref=fastadict[el]
            hamming = hamming_distance(query,ref)
            hamminglist.append(hamming)
            donelist.append(key+el)

#CALCULATE AVERAGE HAMMING DISTANCE
#add all hamming distances divided by number of comparisons
avgham = sum(hamminglist)/len(hamminglist)

#PER 100 nt
#divide by num of 100nt sections in alignment
ref=sorted(fastadict.keys())[0]
refseq=fastadict[ref]
nt100 = len(refseq)/100

#return avg subs/100nt
print("Subs per 100nt: ") 
print(avgham/nt100)

#return pi
print("Pi: )
print(avgham/len(refseq))
