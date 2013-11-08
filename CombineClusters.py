#######################################################################
##
## CombineClusters.py
##
## Version 1.01 -- 7 August 2013
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
import os, math, subprocess

qseq_filename = sys.argv[1]

#sort BLASTsummary file
BLASTsummary_filename =(qseq_filename.replace('.qseq','')+'BLASTsummary.out')
sortedBLAST_filename = (qseq_filename.replace('.qseq','')+'BLASTsorted.out')

command=('sort -k6,6n -k4,4d -k5,5n '+BLASTsummary_filename+' -o '+sortedBLAST_filename)
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

offsets = [0]*101
#construct list of "homologous" clusters
temp_file = open('temp.out', 'w')
sortedBLAST_file = open(sortedBLAST_filename, 'r')
chrom = ('')
position = 0
idlist=[]
for line in sortedBLAST_file:
    data=line.split()
    if data[3]==chrom:
        if data[3]!='.':
            if int(data[4])-int(position) <= 50:
                diff = int(data[4])-int(position)
                offsets[diff]+=1
                idlist.append(int(data[0]))
                position=data[4]
            else:
                if len(idlist)>1:
                    print(idlist)
                    idlist.sort()
                    for i in range(1,len(idlist)):
                        temp_file.write(str(idlist[i])+'\t'+str(idlist[0])+'\t'+str(diff)+'\n')
                        print(idlist[i],idlist[0])
                cluster=data[0]
                idlist=[int(cluster)]
                position=data[4]
                chrom = data[3]
    else:
        if len(idlist)>1:
            print(idlist)
            idlist.sort()
            for i in range(1,len(idlist)):
                temp_file.write(str(idlist[i])+'\t'+str(idlist[0])+'\t'+str(diff)+'\n')
                print(idlist[i],idlist[0])
        cluster=data[0]
        idlist=[int(cluster)]
        position=data[4]
        chrom = data[3]
temp_file.close()

clusterIDs_filename=(qseq_filename.replace('.qseq','')+'clusterIDs.out')    
command=('sort -k1,1n temp.out -o '+clusterIDs_filename)
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

print(offsets)

#update cluster numbers and depths in qseq file
qseq_filename=(qseq_filename.replace('.qseq','temp3.qseq'))
qseq_file=open(qseq_filename,'r')
clusterIDs_file=open(clusterIDs_filename,'r')
out_filename=(qseq_filename.replace('temp3.qseq','temp4.qseq'))
out_file=open(out_filename,'w')

clusterstofix=[]
newclusters=[]
diffs=[]
for line in clusterIDs_file:
    cluster_edit=line.split()
    clusterstofix.append(cluster_edit[0])
    newclusters.append(cluster_edit[1])
    diffs.append(cluster_edit[2]) 
line_count=0
target_count=100000
for line in qseq_file:
    line_count+=1
    if line_count > target_count:
        print(target_count)
        target_count+=100000
    data=line.split()
    if data[13] in clusterstofix:
        #print(data[13],diffs[clusterstofix.index(data[13])])
        data[13]=newclusters[clusterstofix.index(data[13])]
        out_file.write('\t'.join(data)+'\n')
    else:
        out_file.write(line)

out_file.close()

print('Fixed',len(clusterstofix))
                
#sort qseq file
print('Sorting combined qseq file')
out_filename2=(out_filename.replace('temp4.qseq','COMB.qseq'))
command=('sort -k14,14n -k1,1d -k4,4n -k5,5n -k6,6n '+out_filename+' -o '+out_filename2)
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

#cleanup
command=('rm '+qseq_filename+' '+out_filename+' '+sortedBLAST_filename+' temp.out '+clusterIDs_filename)
print(command)
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

print('\nFinished!!\n')
