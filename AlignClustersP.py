#######################################################################
##
## AlignClustersP.py
##
## Version 1.05 -- 8 July 2013
##
## Created by Michael Sorenson
## Copyright (c) 2011-2013 Boston University. All rights reserved.
##
## Dependencies: muscle (version 3 in $PATH as "muscle")
##               text file with list of nodes to use in analysis
##
## This program is written for execution in the python (version 3)
## language, and uses the multiprocessing module and ssh to run
## parallel jobs in a Diskless Remote Boot Linux (DRBL) environment.
## It is free and distributed WITHOUT warranty; without even the
## implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.
##
#######################################################################

import sys
import os, math, subprocess
import multiprocessing
from multiprocessing import Process, Lock

base_filename = str(sys.argv[1])
min_depth = int(sys.argv[2])
host_filename = sys.argv[3]

path = os.getcwd()

def collapse_cluster(fasfile):
    seqsin = open(fasfile,'r')
    seqs=[]
    for line in seqsin:
        if ">" not in line:
            seqs.append(line.strip())
    seqsin.close()
    
    seeds = []
    haps = []
    seed_counts = []
    hap_counts = []
    cluster_list1 = []

    for i in seqs: 
        if i[13:17] not in seeds:
            seeds.append(i[13:17])
            seed_counts.append(1)
            haps.append([i])
            hap_counts.append([1])
            cluster_list1.append([len(seeds)-1,0])
        else:
            list_number = seeds.index(i[13:17])
            seed_counts[list_number]+=1
            if i not in haps[list_number]:
                haps[list_number].append(i)
                hap_counts[list_number].append(1)
                hap_number=len(hap_counts[list_number])-1
            else:
                hap_number = haps[list_number].index(i)
                hap_counts[list_number][hap_number] += 1
            cluster_list1.append([list_number,hap_number])

    count = 0
    clusters = []
    #construct sequentially numbered list of clusters
    for i in range(len(seeds)):
        clusters.append([])
        for j in range(len(haps[i])):
            count += 1
            clusters[i].append(count)

    #create list of new cluster numbers and depths for each sequence
    cluster_list2 = []
    for k in range(len(cluster_list1)):
        i = cluster_list1[k][0]
        j = cluster_list1[k][1]
        cluster_list2.append([clusters[i][j],hap_counts[i][j],haps[i][j]])
    new_muscle=open(fasfile+'s','w')
    uniques=[]
    for i in cluster_list2:
        if i[0] not in uniques:
            uniques.append(i[0])
            new_muscle.write(">"+str(i[0]).zfill(5)+'\n'+i[2]+'\n')
    new_muscle.close()
   
    return (cluster_list2,count)
    #return(count)
        
def did_muscle_finish(filename):
    try:
        testfile = open(filename,'r')
        testfile.close()
        resultcode=1
    except IOError:
        resultcode=0
    return(resultcode)

def run_muscle(clustno,d,s,host,seq_ids):
    cmd = ('ssh '+host+' muscle -in '+path+'/muscle'+clustno+'.fas -quiet -out '+path+'/outalign.'+clustno)
    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)[1]

    testcode = did_muscle_finish('outalign.'+clustno)

    if testcode == 1:
        #input muscle results and sort
        muscle_output=open('outalign.'+clustno,'r')
        aligned_seqs = []
        for line2 in muscle_output:
            if line2[0]==('>'):
                parts2=line2.strip()
                seqID=parts2[1:]
                aligned_seqs.append([seqID,''])
            else:
                parts2=line2.strip('\n')
                aligned_seqs[-1][1]+=parts2
        muscle_output.close()
        
        aligned_seqs.sort()
        
        #merge reformatted muscle file with tempout file
        t_file=open('clustertemp.'+clustno,'r')
        tempout_file=open('tempout.'+clustno,'w')
        line_count=0
        for line3 in t_file:
            line_count+=1
            parts3=line3.split()
            parts3[8]=aligned_seqs[(seq_ids[line_count-1][0]-1)][1]
            found=parts3[8].find('-')
            while found > -1:
                if found == 0:
                    newqual=(parts3[9][0]+parts3[9])
                else:
                    newqual=(parts3[9][:found]+parts3[9][found-1]+parts3[9][found:])
                found = parts3[8].find('-',found+1)
                parts3[9]=newqual
            parts5 = '\t'.join(parts3)
            tempout_file.write(parts5[:parts5.rfind('\t')]+'\t'+str(d)+'\n')
        tempout_file.close()
        
        command=('rm muscle'+clustno+'.fa* outalign*.'+clustno+' clustertemp.'+clustno)
        p = subprocess.Popen(command, shell=True)
        sts = os.waitpid(p.pid, 0)[1]

    else:
        print (clustno,d,s,'muscle failed')

        command=('mv clustertemp.'+clustno+' tempout.'+clustno)
        p = subprocess.Popen(command, shell=True)
        sts = os.waitpid(p.pid, 0)[1]

        command=('rm muscle'+clustno+'.fa*')
        p = subprocess.Popen(command, shell=True)
        sts = os.waitpid(p.pid, 0)[1]

input_file=open(base_filename.replace('.qseq','COMB.qseq'),'r')
muscle_file = open('muscle0.fa','w')

bigclusters = open('BigClusters','w')
output_filename = base_filename.replace('.qseq','A.qseq')
output_file = open(base_filename.replace('.qseq','A.qseq'),'w')
                    
threads=[]
runs=[]
clust2align=[]
clust2write=[]

hosts = []
host_file = open(host_filename,'r')
for line in host_file:
    host = line.strip('\n')
    hosts.append(host)
nthreads = int(len(hosts))
    
line=input_file.readline()
parts = line.split()
oldcluster=parts[13]
temp_file=open('clustertemp.'+oldcluster,'w')
temp_file.write(line)
seqID = (parts[0]+':'+parts[2]+':'+parts[3]+':'+parts[4]+':'+parts[5])
seq = parts[8]
muscle_file.write('>'+seqID+'\n'+seq+'\n')
RAD_depth=int(parts[11])
seqs=1

for line in input_file:
    parts = line.split()
    cluster=parts[13]
    if cluster==oldcluster:
        temp_file.write(line)
        seqID = (parts[0]+':'+parts[2]+':'+parts[3]+':'+parts[4]+':'+parts[5])
        seq = parts[8]
        muscle_file.write('>'+seqID+'\n'+seq+'\n')
        RAD_depth+=int(parts[11])
        seqs+=1
    else:
        print(oldcluster)
        temp_file.close()
        muscle_file.close()
        
        seqids,collapsed = collapse_cluster('muscle'+oldcluster+'.fa')

        bigclusters.write(str(oldcluster)+'\t'+str(seqs)+'\t'+str(RAD_depth)+'\t'+str(collapsed)+'\n')

        if RAD_depth>=min_depth and seqs > 1 and collapsed < 6000:
            if __name__ == '__main__':
                lock = Lock()
                if len(threads) < nthreads:
                    p = multiprocessing.Process(target=run_muscle,args=[oldcluster,RAD_depth,seqs,hosts[len(threads)],seqids])
                    p.start()
                    threads.append(p)
                    clust2align.append(oldcluster)
                    
                else:
                    while len(threads) == nthreads:
                        tnum = -1
                        for thread in threads:
                            tnum += 1
                            if not thread.is_alive():
                                threads.remove(thread)
                                hosts.append(hosts.pop(tnum))
                                clust2write.append(clust2align.pop(tnum))

                    p = multiprocessing.Process(target=run_muscle,args=[oldcluster,RAD_depth,seqs,hosts[len(threads)],seqids])
                    p.start()
                    threads.append(p)
                    clust2align.append(oldcluster)
 
        else:
            if collapsed > 5999:
                print(oldcluster,RAD_depth,seqs,'too many seqs')
            else:
                print(oldcluster,RAD_depth,seqs)
            temp_filewc = open('clustertemp.'+oldcluster,'r')
            for lineT in temp_filewc:
                output_file.write(lineT[:lineT.rfind('\t')]+'\t'+str(RAD_depth)+'\n')
            temp_filewc.close()
     
            command=('rm clustertemp.'+oldcluster+' muscle'+oldcluster+'.fa*')
            p = subprocess.Popen(command, shell=True)
            sts = os.waitpid(p.pid, 0)[1]

            while len(clust2write) > 0:
                filenum = clust2write.pop(0)
                temp_filewc = open('tempout.'+filenum,'r')
                for lineT in temp_filewc:
                    output_file.write(lineT)
                temp_filewc.close()
                
                command=('rm tempout.'+filenum)
                p = subprocess.Popen(command, shell=True)
                sts = os.waitpid(p.pid, 0)[1]

        oldcluster=parts[13]
        temp_file=open('clustertemp.'+oldcluster,'w')
        muscle_file=open('muscle'+oldcluster+'.fa','w')
        temp_file.write(line)
        seqID = (parts[0]+':'+parts[2]+':'+parts[3]+':'+parts[4]+':'+parts[5])
        seq = parts[8]
        muscle_file.write('>'+seqID+'\n'+seq+'\n')
        RAD_depth=int(parts[11])
        seqs=1

temp_file.close()
muscle_file.close()

seqids,collapsed = collapse_cluster('muscle'+oldcluster+'.fa')

bigclusters.write(str(oldcluster)+'\t'+str(seqs)+'\t'+str(RAD_depth)+'\t'+str(collapsed)+'\n')
if RAD_depth>=min_depth and seqs > 1 and collapsed < 6000:
    if __name__ == '__main__':
        lock = Lock()
        if len(threads) < nthreads:
            p = multiprocessing.Process(target=run_muscle,args=[oldcluster,RAD_depth,seqs,hosts[len(threads)],seqids])
            p.start()
            threads.append(p)
            clust2align.append(oldcluster)
            
        else:
            while len(threads) == nthreads:
                tnum = -1
                for thread in threads:
                    tnum += 1
                    if not thread.is_alive():
                        threads.remove(thread)
                        hosts.append(hosts.pop(tnum))
                        clust2write.append(clust2align.pop(tnum))

            p = multiprocessing.Process(target=run_muscle,args=[oldcluster,RAD_depth,seqs,hosts[len(threads)],seqids])
            p.start()
            threads.append(p)
            clust2align.append(oldcluster)

else:                
    if collapsed > 5999:
        print(oldcluster,RAD_depth,seqs,'too many seqs')
    else:
        print(oldcluster,RAD_depth,seqs)
    temp_filewc = open('clustertemp.'+oldcluster,'r')
    for lineT in temp_filewc:
        output_file.write(lineT[:lineT.rfind('\t')]+'\t'+str(RAD_depth)+'\n')
    temp_filewc.close()

    command=('rm clustertemp.'+oldcluster+' muscle'+oldcluster+'.fa*')
    p = subprocess.Popen(command, shell=True)
    sts = os.waitpid(p.pid, 0)[1]

while len(threads) > 0:
    tnum = -1
    for thread in threads:
        tnum += 1
        if not thread.is_alive():
            threads.remove(thread)
            hosts.append(hosts.pop(tnum))
            clust2write.append(clust2align.pop(tnum))

while len(clust2write) > 0:
    filenum = clust2write.pop(0)
    temp_filewc = open('tempout.'+filenum,'r')
    for lineT in temp_filewc:
        output_file.write(lineT)
    temp_filewc.close()
        
    command=('rm tempout.'+filenum)
    p = subprocess.Popen(command, shell=True)
    sts = os.waitpid(p.pid, 0)[1]

input_file.close()
bigclusters.close()
output_file.close()

#sort qseq file
print('Sorting combined qseq file')
command=('sort -k14,14n -k1,1d -k4,4n -k5,5n -k6,6n '+output_filename+' -o '+base_filename.replace('.qseq','AS.qseq'))
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

##cleanup
command=('rm '+base_filename.replace('.qseq','A.qseq')+' '+base_filename.replace('.qseq','temp4.qseq')+' '+base_filename.replace('.qseq','temp3.qseq'))
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

print('\nFinished!!\n')


 

