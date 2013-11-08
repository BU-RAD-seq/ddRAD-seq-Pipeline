#######################################################################
##
## fastq2qseq.py
##
## Version 2 -- 2012
##
## Created by Michael Sorenson
## Copyright (c) 2011-2013 Boston University. All rights reserved.
##
## This program is written for execution in the python (version 3)
## language. It is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
## This program is to convert fastq files (in Cassava v1.8 format) to
## a qseq file format, which is used in the BU-RAD-seq pipeline.
##
#######################################################################

import sys, os, re #imports necessary python modules

##command in terminal:
##python3 fastq2qseq.py infile.fastq

fastq_filename = sys.argv[1] #defines input file
fastq_file = open(fastq_filename, 'r') #opens input file for reading

#define output file name by replacing .fastq with .qseq
qseq_file = open(fastq_filename.replace('.fastq','.qseq'), 'w')

#this sets up a  counter to mark progress
#will print every time 100,000 reads have been processed
line_count=0
target = 100000

#starts loop to go through fastq file
for line in fastq_file:
    line_count += 1
    if line_count == target:
        print(line_count)
        target += 100000

    #splits readname at characters @ : # /  \n
    data = re.split('@|:|#|/| |\n',line)

#comments below define parts of split line
#    name = data[1]+'_'+data[2]
#    run = data[3]
#    coord = data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]
#    index = data[11]
#    read = data[8]
#    filt = data[9]
    
    #get read by reading next line, stripping return character
    seq = fastq_file.readline().strip()

    #skip line 3
    fastq_file.readline()

    #get read quality by reading next line, stripping return character
    qual = fastq_file.readline().strip()

    #write information to qseq output file
    qseq_file.write(data[1]+'_'+data[2]+'\t'+
                    data[3]+'\t'+
                    data[4]+'\t'+data[5]+'\t'+data[6]+'\t'+data[7]+'\t'+
                    data[11]+'\t'+
                    data[8]+'\t'+
                    seq+'\t'+
                    qual+'\t'+
                    data[9]+'\n')

#close input and output files
fastq_file.close()
qseq_file.close()

#print finished to screen so you know script completed
print('\nFinished!!\n')

