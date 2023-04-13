#!/usr/bin/python

import sys #utilize arguments from command line
import os #utilize operating system functions
import re #use regular expressions, RegEx

# The goal of this script is to generate a best reciprocal blast hits for two fasta files,
# which are protein sequences.

# On the command line write: python rbh.py [query.fasta] [subject.fasta] [output.txt]

# Save the arguments from the command line, starts at 0 with the python script; ignore 0
query = sys.argv[1]
subject = sys.argv[2]

# Check
print ("Query sequence is:",query)
print ("Database is built from:",subject)

# Create a directory to store results

wd_path = os.getcwd() #get the current directory
path_db = os.path.join(wd_path,"db_results") #add the db folder name

os.mkdir (path_db)

# Check
print ("Created database folder:",path_db)

# Create database
database = query + "_db"

cmd_mkdb = "makeblastdb -in " + subject + " -dbtype 'prot' -out ./db_results/" + database
os.system(cmd_mkdb)

# Database location
database_file = os.path.join(path_db,database)

# Run blastp with first query 

out1 = query + "_vs_" + subject + ".txt"

cmd_query = "blastp -db " + database_file + " -query " + query + " -out " + out1 + " -outfmt 7"

os.system(cmd_query)

####################

# Now do the reciprocal blast run
query1 = sys.argv[2]
subject1 = sys.argv[1]

# Make the second database
database1 = query1 + "_db"

cmd_mkdb1 = "makeblastdb -in " + subject1 + " -dbtype 'prot' -out ./db_results/" + database1
os.system(cmd_mkdb1)

# Second database location
database_file1 = os.path.join(path_db,database1)

# Run blastp with second query 
out2 = query1 + "_vs_" + subject1 + ".txt"

cmd_query1 = "blastp -db " + database_file1 + " -query " + query1 + " -out " + out2 + " -outfmt 7"

os.system(cmd_query1)

####################

# Parse the first BLAST results
blast1 = open(out1,'r') #read ('r') tab-delimited file
D1 = {} #create a dictionary for the first blast output file

for Line in blast1:
	if (Line[0] != '#'): #ignore comments
		Line.strip() #strip each line
		Elements = re.split('\t',Line) #strip by tab-delimited columns, and save as elements
		queryID = Elements[0]
		subjectID = Elements[1]
		if (not (queryID in D1.keys())): #add the queryID to the dictionary
			D1[queryID] = subjectID #pick the first subjectID from the database, which is the unidirectional best hit

#####################

# Parse the second BLAST results; repeat above but with blast2
blast2 = open(out2,'r')
D2 = {}

for Line in blast2:
	if (Line[0] != '#'):
		Line.strip()
		Elements = re.split('\t',Line)
		queryID = Elements[0]
		subjectID = Elements[1]
		if ( not(queryID in D2.keys())):
			D2[queryID] = subjectID

#####################

# Identify the shared pairs, i.e. the reciprocal best blast hits

sharedPairs = {} #new dictionary to identify shared ids

for id1 in D1.keys():
	value1 = D1[id1]
	if (value1 in D2.keys()):
		if (id1 == D2[value1]):
			sharedPairs[value1] = id1
			
saved_file = sys.argv[3]
			
outfile = open(saved_file,'w') #write RBHs

for k1 in sharedPairs.keys():
	line = k1 + '\t' + sharedPairs[k1] + '\n'
	outfile.write(line)
	
outfile.close()

print ("Done! RBH from ",query,"and ",subject,"have been saved to ",saved_file,".")