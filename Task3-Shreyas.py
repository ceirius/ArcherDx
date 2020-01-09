#!/usr/bin/env python
# coding: utf-8

# 
# Given a chromosome and coordinates, write a program for looking up its annotation.  Keep in mind you'll be doing 
# this annotation millions of times.
# 
# -Input: 
# o   Tab-delimited file: Chr<tab>Position
# o   GTF formatted file with genome annotations.
# 
# -Output: 
# o   Annotated file of gene name that input position overlaps.
# 
# -Hint: Most of the sequence reads come from a small portion of the genome. Try to use this information to 
# improve performance, if possible.
# 
# 



# python v3.7.4

import argparse
import numpy as np  # numpy Version: 1.16.5
import pandas as pd # pandas Version: 0.25.1
import os, re, glob
from tqdm import tqdm
from io import StringIO
from datetime import datetime



"""Command line argument parser: """

parser = argparse.ArgumentParser(description='Look up SNPs against gtf annotations. Processes the gtf file in chunks of 30000.', prog='python Task3.py --gtf [sample.gtf] --in [sample.txt]')
parser.add_argument('--gtf', help='input the gtf file name with .extension, example: sample1.gtf', dest='input1', type=str, required=True)
parser.add_argument('--in', help='input the SNP file name with .extension, example: sample2.txt', dest='input2', type=str, required=True)
args = parser.parse_args()



"""getting input file data from user"""

directory1 = args.input1
directory2 = args.input2



"""getting input file data from user"""

"""directory1 = input("Enter file name for gtf file, followed by extension\t\t")
directory2 = input("Enter file name for annotate, followed by extension file\t\t")"""




"""Evaluate the path for the specific file names and then evaluate that the files can be opened and read. 
Does not evaluate if the correct file was entered in correct order, or if the files are of type gtf and annotations."""

try:
    """finding the input files paths"""

    path1 = os.path.split(os.path.abspath(directory1))
    path1_full = os.path.join(path1[0] + "\\**\\" + path1[1])
    final_path1 = glob.glob(path1_full, recursive = True)

    path2 = os.path.split(os.path.abspath(directory2))
    path2_full = os.path.join(path2[0] + "\\**\\" + path2[1])
    final_path2 = glob.glob(path2_full, recursive = True)
    
    """opening to read contents of the files"""

    file_content1 = open(final_path1[0], mode='r').read()
    file_content2 = open(final_path2[0], mode='r').read()
    
    print('\nBoth file names you provided are valid\n\n')
    
except:
    print('\nOne or both of the file names you provided are NOT VALID\n\n')




"""Creating dataframes out of the input files. By having a default chunk size, we can manage memory usage.""" 
df1 = pd.read_csv(StringIO(file_content1), sep = '\s+', header=None, chunksize = 30000)
df2 = pd.read_csv(StringIO(file_content2), sep = '\s+', header=None)



df2['annotate_Index'] = df2.index
#df2.drop_duplicates(subset=[1,0], keep='first', inplace=True)



"""create a new dataframe with labelled columns and the contained output"""
new_df = pd.DataFrame(columns=['Chromosome', 'SNP', 'Chr copy', 'Gene', 'Region', 'Exon#', "Transcript ID", "Exon ID", 'Orientation', 'Start', 'End'])



def append_True(coords_list, x,y,z,w, final_array):
    
    """idx 4 is 2D so that one dimension gives the 0 - 1, and the second dimension == len(coords). 
    So the condition takes item4 as index position of df2 and df1"""

    for idx4, item4 in enumerate(coords_list):
        
        temp_list = []
        
        
        if idx4 < len(coords_list)//2:
            
            temp_list.append(w.T[0][coords_list[idx4]]) 
            temp_list.append(z.T[0][coords_list[idx4]])
            temp_list.append(y.T[5][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(y.T[0][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(y.T[4][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(y.T[1][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(y.T[2][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(y.T[6][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(y.T[3][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(x.T[0][coords_list[idx4 + len(coords_list)//2]])
            temp_list.append(x.T[1][coords_list[idx4 + len(coords_list)//2]])  ### these two lines are to test that indices match up
            #temp_list.append(z.T[1][coords_list[idx4]])
        
        
        """This condition executes for each item in the coords_list3; so it's (array p) shape should be 12,1
        final_array will keep growing for every mapped SNP"""
        if len(temp_list) > 0:
            p =  np.reshape(np.array(temp_list), (1,11))
            
            final_array = np.append(final_array, p.T, axis=1) # append the array p to the final_array.
            
        del(temp_list)

        
    return(final_array)
        


def append_False(coords_list3, z,w, final_array5, new_df):
    
    """corods_list3 is 1D so that one dimension idx4 merely gives index position of unmapped SNPs in the list
    item4 gives the index position of the unmapped SNPs in the input df2"""
    
    for idx4, item4 in enumerate(coords_list3):
        
        temp_list = []
        
        """if value at index item4 z.T array is not in df2.SNP """
        if z.T[0][item4] not in new_df['SNP']:
            
            """append valued of df2 at that position into a temp list"""
            temp_list.append(w.T[0][item4]) 
            temp_list.append(z.T[0][item4])
            temp_list.append(".")
            temp_list.append(".")
            temp_list.append(".")
            temp_list.append(".")
            temp_list.append(".")
            temp_list.append(".")
            temp_list.append(".")
            temp_list.append(".")
            temp_list.append(".")
            #temp_list.append(z.T[1][item4])
        
        
        """This condition executes for each item in the coords_list3; so it's shape should be 12,1"""
        if len(temp_list) > 0:
            p2 =  np.reshape(np.array(temp_list), (1,11))
            final_array5 = np.append(final_array5, p2.T, axis=1)
            
    
        del(temp_list)          # to make sure the list is clean in the next iteration
    
    return(final_array5)



def update_DF(new_df, final_array4):
        
    """create a new dataframe with labelled columns and append with the final_array the contained output"""

    new_df1 = pd.DataFrame.from_records(final_array4.T, index=None, exclude=None, columns=['Chromosome', 'SNP', 'Chr copy', 'Gene', 'Region', 'Exon#', "Transcript ID", "Exon ID", 'Orientation', 'Start', 'End'])
    
    new_df = new_df.append(new_df1)
    
    return(new_df)



def update_unmapped_SNPs(new_df, z, w):
    
    
    """identifying items in df2 that were not mapped to the gtf"""
    false_array = (new_df.SNP.astype("int64")[np.newaxis, :] == df2[1][:, np.newaxis]) & (new_df.Chromosome[np.newaxis, :] == df2[0][:, np.newaxis])
    
    """evaluate items in the input coordiantes file that have and have not been mapped (boolean matrix)"""
    coords2 =  false_array.sum(axis = 1)
    
    """extract only those coordinates where the axis sums to 0, i.e. there is not match and has not yet mapped"""
    coords_list2 = np.ravel([np.where(coords2 == 0)])
    
    """creating a new array to populate with unmapped SNPs"""
    final_array5 = np.zeros((11,0), dtype = object) 
    
    """calling function to map unmapped SNPs"""
    final_array5 = append_False(coords_list2, z,w, final_array5, new_df)
    
    """update the results dataframe"""
    new_df = update_DF(new_df, final_array5)
        
    return(new_df)



for chunk in tqdm(df1, desc="Mapping SNPs to chromosomal coordinates. ~2 sec per iteration.(Iteration#, Elapsed time, Time per iteration)"):
    
    """creating a final_array of type object to populate with values from arrays x and y
    the shape of array equals len(df2) and 6 columns that I am extracting. We can adjust this later as needed"""

    final_array = np.zeros((11,0), dtype = object)


    """assigning the index value to a column and then to the downstream array will help track SNP-gtf matches"""
    chunk['GTF_Index'] = chunk.index
    df1a = chunk                   # this could be deprecated
    
    
    """creating np arrays from df1 and df2"""
    x = df1a[[3, 4, 'GTF_Index']].to_numpy()               # integer array
    y = df1a[[9, 13, 11, 6, 2, 0, 15]].to_numpy()              # object array
    z = df2[[1, 'annotate_Index']].to_numpy()              # integer array
    w = df2[[0]].to_numpy()                                # object array
    
    
    """creating a boolean matrix"""
    n = (y.T[5] == w.T[0][:,np.newaxis]) & (x.T[0] <= z.T[0][:,np.newaxis]) & (x.T[1] >= z.T[0][:,np.newaxis])
    
     
    """creating an array for coordinates. the coordinates respectively represent the index positions in df2 and df1"""
    coords1 = np.where(n==True)  # creating an array of the coordinates where array n is True
    
    
    """creating a list from the unraveled array coords"""
    coords_list1 = list(np.ravel(coords1))
    
    
    """calling function to append SNPs that matched gtf to new_df"""
    final_array = append_True(coords_list1, x,y,z,w, final_array)
    
    new_df = update_DF(new_df, final_array)
    

"""calling function to append unmapped SNPs that did not match gtf to new_df"""
new_df = update_unmapped_SNPs(new_df, z, w)



"""sort the output data by chromosome and SNP"""

new_df.sort_values(by = ['Chromosome','SNP', 'Exon#', 'Transcript ID'], axis=0, inplace=True)



"""clean up functions"""

def changeType(val):
	"""convert numeric columns to interger type"""
	
	val = val.astype("int")
	return(val)

def cleanup(valu):
	"""stripping the ';' from column values"""
	
	return valu.rstrip(';')


"""cleaning up output data and converting numbers to dtype int"""

for value in new_df.columns:
    
    new_df[value] = new_df[value].apply(cleanup)
    
    if re.match(r'^\d+$', str(new_df[value][0])):
        new_df[value] = changeType(new_df[value])


"""each output is saved with a date and time stamp"""

filename = "SNP-Chr_map_" + str(datetime.now().strftime("%Y%m%d_%H%M")) + ".csv"

filename =os.path.join(path2[0] + "\\" + filename)



"""Conclusion - save the SNP-GTF map to a csv file, output a sample to screen, and END"""

new_df.to_csv(filename, sep = "\t", index = False)

print("\n\nOutput file saved as (%s)\n" %filename, "\n\nSample view\n\n", new_df.sample(n=5))