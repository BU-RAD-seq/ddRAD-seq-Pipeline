#######################################################################
##
## CondenseSequences.py
##
## Version 1.01 -- 8 November 2013
##
## Created by Michael Sorenson
## Copyright (c) 2011-2013 Boston University. All rights reserved.
##
## This program is written for execution in the python (version 3)
## language. It is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
#######################################################################

import sys
import os, subprocess

def complexity(seq):
    triplets=[]
    triplet_counts=[]
    for i in range(len(seq)-2):
        triplet = seq[i:i+3]
        if triplet in triplets:
            triplet_counts[triplets.index(triplet)]+=1
        else:
            triplets.append(seq[i:i+3])
            triplet_counts.append(1)

    S = 0.0
    if len(triplets) > 1:
        for i in range(len(triplets)):
            S += float(triplet_counts[i])*(triplet_counts[i]-1)/2
        S = round(S/(len(triplets)-1),2)
    else:
        S = 100

    return(S)

def collapse_sequences(filename):
    input_file = open(filename,'r')
    line_count=0
    seeds = []
    haps = []
    seed_counts = []
    hap_counts = []
    cluster_list1 = []

    interval = 100000
    target_count = 100000
    for line in input_file:
        line_count += 1
        if line_count > target_count:
            print(target_count)
            target_count += interval
        data = line.split()
        if data[8][13:17] not in seeds:
            seeds.append(data[8][13:17])
            seed_counts.append(1)
            haps.append([data[8]])
            hap_counts.append([1])
            cluster_list1.append([len(seeds)-1,0])
        else:
            list_number = seeds.index(data[8][13:17])
            seed_counts[list_number]+=1
            if data[8] not in haps[list_number]:
                haps[list_number].append(data[8])
                hap_counts[list_number].append(1)
                hap_number=len(hap_counts[list_number])-1
            else:
                hap_number = haps[list_number].index(data[8])
                hap_counts[list_number][hap_number] += 1
            cluster_list1.append([list_number,hap_number])

    count = 0
    clusters = []
    #construct sequentially numbered list of clusters
    print('Renumbering clusters')
    for i in range(len(seeds)):
        clusters.append([])
        for j in range(len(haps[i])):
            count += 1
            clusters[i].append(count)
    print(count,'clusters')

    #create list of new cluster numbers and depths for each sequence
    cluster_list2 = []
    for k in range(len(cluster_list1)):
        i = cluster_list1[k][0]
        j = cluster_list1[k][1]
        cluster_list2.append([clusters[i][j],hap_counts[i][j]])

    return (cluster_list2,count)


qseq_filename = sys.argv[1]
print(qseq_filename)

#collapse identical sequences
print('Collapsing identical sequences')
Clist,Ccount = collapse_sequences(qseq_filename)

#add cluster numbers and depth to qseq-like file
print('Adding depth to qseq-like file')
qseq_file = open(qseq_filename,'r')
temp_qseq_file = open(qseq_filename.replace('.qseq','temp.qseq'), 'w')
    
line_count=0
for line in qseq_file:
    line_count += 1
    temp_qseq_file.write(line[:line.rfind('\n')] + '\t'+str(Clist[line_count-1][1])+'\t'+str(Clist[line_count-1][0])+'\n')

qseq_file.close()
temp_qseq_file.close()

#sort combined qseq file
print('Sorting combined qseq file')
command=('sort -k13,13n -k1,1d -k4,4n -k5,5n -k6,6n '+qseq_filename.replace('.qseq','temp.qseq')+' -o '+qseq_filename.replace('.qseq','sort.qseq'))
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

print('Select best seq per cluster for downstream analysis')

fileinput = open(qseq_filename.replace('.qseq','sort.qseq'),'r')
qseqsort_file = open(qseq_filename.replace('.qseq','C.qseq'),'w')

oldcluster=-1
bestseq=['']
maxscore=0
for line in fileinput:           
    data1 = line.split()
    newcluster=int(data1[12])
        
    if newcluster != oldcluster:
        if len(bestseq)>2:
            score = 0 
            for a in maxquals: 
                score += ord(a)-33
            score = round(score/len(maxquals),2)
            bestseq[9]=''.join(maxquals)
            S = complexity(bestseq[8])
            qseqsort_file.write('\t'.join(bestseq[:12])+'\t'+str(score)+'\t'+str(S)+'\n')
        oldcluster=newcluster
        bestseq=line.split()
        maxquals = bestseq[9]    
    else:
        tempqual = []
        for a,b in zip(maxquals,data1[9]):
            if ord(b)>ord(a):
                tempqual.append(b)
            else:
                tempqual.append(a)
        maxquals = tempqual
                
score = 0 
for a in maxquals: 
    score += ord(a)-33
score = round(score/len(maxquals),2)
bestseq[9]=''.join(maxquals)
S = complexity(bestseq[8][int(sys.argv[2]):])
qseqsort_file.write('\t'.join(bestseq[:12])+'\t'+str(score)+'\t'+str(S)+'\n')
qseqsort_file.close()

sumfile = open(sys.argv[1].replace('.qseq','.sumtemp'),'w')
sumfile.write(qseq_filename+'\t'+str(line_count)+'\t'+str(Ccount)+'\n')
sumfile.close()

#cleanup
command=('rm '+qseq_filename.replace('.qseq','sort.qseq')+' '+qseq_filename.replace('.qseq','temp.qseq'))
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]


