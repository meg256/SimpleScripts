#This script is designed to produce a labeled database of isolates with ATs assigned.
#It takes in a redundant fasta and outputs a csv file with three columns: IsolateName,AT,sequence
#It also outputs a redundant fasta with all headers changed to AT_IsolateName and a nonredundant fasta with AT names as headers.

#Your redundant fasta should be a multifasta, with your isolate IDs and their associated typing region sequence.
#Bases should be in all caps, although this will only make a difference in two situations: 1) if you have an N ambiguity code, or 2) you have
#two copies of the same sequence, where one is in lowercase and the other is in uppercase. This script would read them as different sequences
#and consider them unique ATs.

import re
import os
import sys


#read in the two files
filepath = os.getcwd()
redundantfasta = sys.argv[1]

#redundant fasta
redf = open(filepath+'//'+redundantfasta)
redlines = redf.readlines()
redf.close()


#make a dictionary in format header:sequence.
fastadict={}
for i in range(0,len(redlines)):
	if re.match('>',redlines[i]):
		key=redlines[i].strip('\n')
		val=redlines[i+1].strip('\n')
		fastadict[key]=val

#get a list of unique sequences in this dictionary.
uniqueseqs=list(set(fastadict.values()))


#Assign an AT to each unique sequence.
#Make a new dictionary with format uniqueseq:AT
atdict = {}
for i in range(0,len(uniqueseqs)):
	key=uniqueseqs[i]
	val=i+1
	atdict[key]=str(val)


#Now iterate through each header in the fastadict dictionary. Write each IsolateName,AT,sequence to a line in a list of lines.
csvlines = []
csvlines.append('IsolateName,AT,sequence\n')
for key in fastadict.keys():
	line=''
	name=key.strip('>')
	seq = fastadict[key]
	atnum = atdict[seq]
	line = name+','+atnum+','+seq+'\n'
	csvlines.append(line)

#first output file: a csv in the format IsolateName,AT,sequence named "database.csv"
with open('ATdatabase.csv','w') as databasefile:
	for line in csvlines:
		databasefile.write(line)

#second output file: a redundant fasta with header format AT_IsolateName
fastalines=[]
for key in fastadict.keys():
	line=''
	name=key.strip('>')
	seq = fastadict[key]
	atnum = 'AT'+str(atdict[seq])
	header='>'+atnum+'_'+name
	line = header+'\n'+seq+'\n'
	fastalines.append(line)

with open('AT_labeled_fasta.fasta','w') as atlabelfile:
	for line in fastalines:
		atlabelfile.write(line)

#third output file: nonredundant fasta with AT labels
nonredlines=[]
for key in atdict.keys():
	seq=key
	header='>AT'+atdict[key]
	nonredlines.append(header+'\n'+seq+'\n')

with open('ATlabel_nonredundant.fasta','w') as nonredfile:
	for line in nonredlines:
		nonredfile.write(line)

databasefile.close()
atlabelfile.close()
nonredfile.close()