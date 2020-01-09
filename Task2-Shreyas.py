#!/usr/bin/env python
# coding: utf-8

# In[ ]:


"""
Given a FASTA file with DNA sequences, find 10 most frequent sequences and return the sequence and their counts 
in the file. 
"""


# In[2]:


"""importing libraries"""

import numpy as np, pandas as pd
import os, glob
from io import StringIO
import re
from collections import Counter 
import argparse



"""Command line argument parser: """

parser = argparse.ArgumentParser(description='Look up 10 most frequent sequences', prog='python Task2.py --in[sample.fasta]')
parser.add_argument('--in', help='input the fasta file name with .extension, example: sample.fasta', dest='input', type=str, required=True)
args = parser.parse_args()


"""getting input file data from user"""

directory1 = args.input


# In[8]:


#file_content1 = open('C:\\Users\\shreyas\\Downloads\\data analysis datasets\\ArcherDX\\coding_tasks\\sample_files\\fasta\\sample.fasta', mode='r').read()



"""getting input file data from user"""

#directory1 = input("Enter file name for fasta file, followed by extension:\t\t")



"""evaluate that the enterd file name and extension are valid and can be found"""

try:
    """finding the input files paths"""

    path1 = os.path.split(os.path.abspath(directory1))
    path1_full = os.path.join(path1[0] + "\\**\\" + path1[1])
    final_path1 = glob.glob(path1_full, recursive = True)

    
    """opening to read contents of the files"""
    
    with open(final_path1[0], mode='r') as f:
        file_content1 = f.read()
    
    print('The file name you provided is valid')
    
except:
    print('One of the file names you provided is NOT VALID')



"""reading file into a dataframe"""

df1 = pd.read_csv(StringIO(file_content1), sep = '\t', header=None)



"""converting the dataframe into a numpy array"""

df1_array = df1.to_numpy()



"""reshaping the numpy array"""

df1_array2 = df1_array.reshape(-1,2)



"""creating a counter object for the sequences"""

x = Counter(df1_array2.T[1].tolist())



"""identifying the most common sequences"""

x_freq = np.array(x.most_common(10), dtype = object)



"""building a Treu False matrix of matches between the most frequent sequences and those in the original input"""

seqID_TF = df1_array2.T[1][:,np.newaxis] == x_freq.T[0]



"""getting the coordinates in the counter object and the original input array of the positions/sequences that matched"""

seq_ID_coords = np.where(seqID_TF == True)



"""creating an output array. The list will be populated with sequences, counts, and seq IDs"""

out_array = np.zeros((10,2), dtype = object)



"""assigning the sequences and their counts to the new output array"""

out_array.T[0] = x_freq.T[0]
out_array.T[1] = x_freq.T[1]



"""output the results to a dataframe and a csv file"""

new_df = pd.DataFrame.from_records(out_array, index=None, exclude=None, columns=['Sequence', 'Count'])
new_df.to_csv("task2_out.csv", sep="\t", index=False)
print(new_df.head(n=5))