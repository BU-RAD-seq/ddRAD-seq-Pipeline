#######################################################################
##
## uCLust85toBLAST.py
##
## Version 1.00 -- 2011
##
## Created by Michael Sorenson
## Copyright (c) 2011-2013 Boston University. All rights reserved.
##
## Dependencies: usearch (version 5 in $PATH as "usearch")
##
## This program is written for execution in the python (version 3)
## language. It is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
#######################################################################

import sys
import os, math, subprocess

qseq_filename = sys.argv[1]
qseq_file = open(qseq_filename,'r')

min_depth = int(sys.argv[2])

#sort by quality
print('Sorting qseq file by quality')
tempfile1 = qseq_filename.replace('.qseq','temp1.qseq')
command=('sort -k13,13nr '+qseq_filename+' -o '+tempfile1)
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

#make fasta file
print('Making fasta file for uClust')
command = ("awk '//{print(\">\"$1\":\"$4\":\"$5\":\"$6\"\\n\"$9)}' "+tempfile1+" > temp.fas")
print(command)
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

#run uclust on fasta file
print('Running uclust on fasta file')
uc_results=(qseq_filename.replace('.qseq','')+'.uc85')
command = ('usearch --usersort --nofastalign --iddef 1 --cluster temp.fas --uc '+uc_results+' --id 0.85 --minlen 30')
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

#get clusterdepths
print('Getting cluster depths')
command=("awk '/C\t/{print$2,$3}' "+uc_results+" > clusterdepth")
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

#Adding Cluster # and depth to qseq-like file
print('Adding Cluster # and depth to qseq-like file')
qseq_file = open(tempfile1,'r')
ucresults_file = open(uc_results,'r')
tempfile2 = (qseq_filename.replace('.qseq','')+'temp2.qseq')
temp_qseq_file = open(tempfile2, 'w')

depths=[]
clusterdepth_file = open('clusterdepth','r')
for line in clusterdepth_file:
    clusterdata = line.split()
    depths.append(int(clusterdata[1]))
    
for i in range(8):
    ucresults_file.readline()
    
line_count=0
target_count=100000
for line in qseq_file:
    line_count += 1
    if line_count > target_count:
        print(target_count)
        target_count += 100000
        
    ucresult = ucresults_file.readline()
    ucresults = ucresult.split()
    d_index=int(ucresults[1])
    temp_qseq_file.write(line[:line.rfind('\n')]+'\t'+str(ucresults[1])+'\t'+str(depths[d_index])+'\n')

qseq_file.close()
ucresults_file.close()
temp_qseq_file.close()

#sort combined qseq file
print('Sorting combined qseq file')
tempfile3 = qseq_filename.replace('.qseq','')+'temp3.qseq'
command=('sort -k14,14n -k1,1d -k4,4n -k5,5n -k6,6n '+tempfile2+' -o '+tempfile3)
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]
command=('rm '+tempfile1+' '+tempfile2+' clusterdepth temp.fas')
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

#Select best seq per cluster for BLAST
print('Select best seq per cluster for BLAST')
inputfile = open(tempfile3,'r')
blastfile = open(qseq_filename.replace('.qseq','BLAST.fasta'),'w')

oldcluster=-1
maxscore=0
linecount=0
target_count = 100000
seq_count=0
for line in inputfile:
    linecount+=1
    if linecount > target_count:
        print(target_count)
        target_count += 100000
        
    data1 = line.split()
    newcluster=data1[13]
    score = float(data1[12])
        
    if newcluster != oldcluster:
        if seq_count >= min_depth:
            best_parts = bestseq.strip('\n').split('\t')
            blastfile.write('>' + best_parts[0] + ":" + best_parts[13] + ':'+str(maxscore)+'\n' + best_parts[8] + '\n')
        oldcluster=newcluster
        maxscore=float(data1[12])
        bestseq=line
        seq_count=int(data1[11])
        
    else:
        seq_count+=int(data1[11])
        if score > maxscore:
            maxscore = score
            bestseq = line 
        
if seq_count >= min_depth:
    best_parts = bestseq.strip('\n').split('\t')
    blastfile.write('>' + best_parts[0] + ":" + best_parts[13] + ':'+str(maxscore)+'\n' + best_parts[8] + '\n')

blastfile.close()

print('\nFinished!!\n')

