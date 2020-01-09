# ArcherDx

The following modules are required to run these scripts: argparse, os, io, collections, datetime, glob, re, numpy, pandas

Task1: Task1-Shreyas.py
execute as python Task1-Shreyas.py from the terminal
output file: fastq list-reads grtr 30nt.csv

Task2: Task2-Shreyas.py
execute as python Task2-Shreyas.py --in sample.fasta from the terminal
output file: task2_out.csv

Task3: Task3-Shreyas.py
execute as python Task3.py --gtf [sample.gtf] --in [sample.txt] from the terminal
output file: Task3 - SNP-Chr_map_date_timeGMT.csv

Task3: 
Since the sample files are presumed to be dummy files (for all tasks) I have not presumed errors in the input files. i.e.  I have not dropped duplicate entries in both the coordinates/SNP file or in the gtf. 
Task3 accounts for SNPs that do not map anywhere in the gtf file (novel SNPs) in addition to those that map (to multiple transcripts and duplicates) to the gtf. 
The input file is taken in in chunks of 30000 lines and completes the mapping in under a minute. Chunk size can be altered in the code (line 86), not as an argument at the command line. 
This script completes in under a minute and takes about 700 Mb memory.
Task3 does not sort the input file to process the mapping by chromosome or region of chromosome as theoretically and simplistically using chromosomal coordinate ranges could include an entire chromosome anyway. Such a situation would reduce to the solution I have provided. However, given that you work with markers for specific conditions it is reasonable to further optimize the code to perform the mapping by chromosomal region.
