#!/usr/bin/env python
# coding: utf-8


"""
Shreyas Krishnan; ArcherDX 
1) Recursively find all FASTQ files in a directory and report each file name and the percent of sequences 
in that file that are greater than 30 nucleotides long.
"""



import os, re, glob
import numpy as np
import pandas as pd
from tqdm import tqdm
import argparse



"""joining the path obtained above with the intermediate unknown path and the user given extension"""
path2 = os.path.split(os.path.abspath('fastq'))[0]

final_path = os.path.join(path2 + "\\**\\*." + 'fastq') # joining the user given path and extension



""" function to count number of fastq entries"""

def num_lines(x):
    
    counter2 = 0
    
    for line in x:
        if re.search(r'^@', line):
            counter2 += 1
    
    if counter2 == 1:
        return(counter2)
    else: return(counter2-2)



"""function to count the number of sequences of length > 30 nt in each file"""

def seq_len(y, q1):
    
    list_len = []
    
    # make an array out of the file
    y_array = np.array(y.split("\n"))
    
    y_array2 = np.reshape(y_array, (-1,4))   # the -1 allow the function to determine the second dimension automatically
    
    for m,n in np.ndenumerate(y_array2.T[1]):  # transposing the new array puts all sequences in one column
        
        if len(n)>30:
            list_len.append(len(n))
    
    return(len(list_len))



"""
iterating through the directory tree searching for the files with extension
splitting the discovered file paths into dir name and base name and populating the array
Then open the array and count number of @ symbols and append that to a new axis in the array as the number of Fastq statements
call a function to count string length of sequence
"""

file_list = []

file_array = np.array(file_list, dtype=str) # an array to contain the dirname, basename for each file

line_counts = np.zeros((0), dtype = int)

num_files = 0

for i in tqdm(glob.iglob(final_path, recursive=True), desc = "parsing files (File number, Elapsed time, Time per file)"):   # iteratively going through folders to find files with the extension in the path
    
    num_files += 1
    
    print(num_files, "\t", os.path.basename(i), "\t")
    
    line_list = []   # a lit that will eventually contain the number of reads for a file
    
    seq_length = []  # a list that will eventually contain the number of seq > 30 nt for each read for a file
    
    file_array = np.append(file_array, os.path.split(i), axis = 0) # populating the array with dirname and basename 
                                                                   # of all files that meet the extension condition

    with open(i) as file:   # open each file in path one at time
        
        file_open = file.read()    # read the file contents 
        
        q = num_lines(file_open)  # calling function to count number of reads and returning the value into q
        
        line_list.append(q)       # appending q number of reads into the list
        
        seq_length.append(seq_len(file_open, q))  # calling seq_len function, sending q number of lines value and 
                                                  # returning and appending the returned length of each read into the list
        file.close() # close the file
		
    line_counts = np.append(line_counts, np.array(line_list).T, axis = 0) # 
    line_counts = np.append(line_counts, np.array(seq_length).T, axis = 0) # 



"""output for below line. top row is number of lines. lower row is number of reads that meet the condition."""

line_counts2 = np.transpose(line_counts.reshape((-1,2), order='C'))



file_array2 = file_array.reshape((-1,2), order='C')



"""calculating the percentage of read > 30 nt"""

line_counts2[1] = line_counts2[1]/line_counts2[0] *100



"""creating a dataframe to output results"""
df = pd.DataFrame(index = range(num_files), columns = ["Filename", "Number of reads in file", "Percentage of reads > 30 nt"])


"""assigning output to dataframe"""
df['Filename'] = file_array2.T[1]
df['Number of reads in file'] = line_counts2[0]
df['Percentage of reads > 30 nt'] = line_counts2[1]

"""output file to the following path"""
fname = os.path.join(path2 + "\\fastq list-reads grtr 30nt.csv")

df.to_csv(fname, sep = "\t", index = False, header = True)
print("\nOutput saved to fastq list, reads grtr 30nt.csv at", fname)