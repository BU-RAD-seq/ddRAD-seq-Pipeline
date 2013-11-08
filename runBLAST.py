#######################################################################
## runBLAST.py
##
## Version 1.01 -- 7 August 2013
##
## Created by Michael Sorenson
## Copyright (c) 2011-2013 Boston University. All rights reserved.
##
## Dependencies: blastn (version 2 in $PATH as "blastn")
##               blastn database for reference genome
##
## This program is written for execution in the python (version 3)
## language. It is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
#######################################################################

import sys
import os, subprocess, heapq

basefilename = sys.argv[1].replace('.qseq','')
BLASTdatabase = sys.argv[2] #enter with path - e.g., /home/Zeeb/Zeeb324

infile = open(basefilename+'BLAST.fasta','r')

target_length = 40

#create input files with small number of seqs per file to speed BLAST searches
line_count = 0
file_count = 10001
target = target_length
outfile = open(basefilename+str(file_count)+'B.fas','w')
for line in infile:
    line_count+= 1
    outfile.write(line)
    if line_count == target_length:
        outfile.close()
        file_count+=1
        line_count=0
        outfile = open(basefilename+str(file_count)+'B.fas','w')
outfile.close()
infile.close()


#BLAST sequences against specified database
for i in range(10001,file_count+1):
    print(basefilename+str(i))
    command = ('blastn -query '+basefilename+str(i)+'B.fas -db '+sys.argv[2]+' -out '+basefilename+str(i)+'.Bout -evalue 0.0001 -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 1 -outfmt "7 std qlen" -dust yes -num_threads 1')
    p = subprocess.Popen(command, shell=True)
    sts = os.waitpid(p.pid, 0)[1]

command = ('cat *.Bout > '+basefilename+'BLAST.results')
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]

command = ('rm *.Bout *B.fas')
p = subprocess.Popen(command, shell=True)
sts = os.waitpid(p.pid, 0)[1]


#parse results file
infile = open(basefilename+'BLAST.results','r')
outfile = open(basefilename+'BLASTsummary.out','w')

for line in infile:
    if 'Query:' in line:
            parts=line.split(':')
            cluster=parts[2]
    else:
        if 'hits' in line:
            parts=line.split(' ')
            hits=int(parts[1])
            if hits==0:
                outfile.write(cluster+'\t0\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n')
            else:                                    
                evalues=[]
                hits_data = []
                for i in range(hits):
                    line = infile.readline()
                    if line[0]=='#':
                        break
                    else:
                        hits_data.append(line.split())
                        evalues.append(float(hits_data[i][10]))
                    
                best_hits=[]
                if hits == 1:
                    result_code = 'one'
                    best_hits.append(hits_data[0])
                else:
                    #find second best hit
                    evalue1,evalue2 = heapq.nsmallest(2,evalues)
                    
                    for line in hits_data:
                        if float(line[10]) <= evalue2:
                            best_hits.append(line)

                    if evalue2 > evalue1:
                        if best_hits[1][1].lower().find('random') > -1:
                            for i in range(2,len(best_hits)):
                                if best_hits[i][1].lower().find('random') < 0:
                                    best_hits.insert(1,best_hits.pop(i))
                                    break
                    else:
                        if best_hits[0][1].lower().find('random') > -1:
                            for i in range(1,len(best_hits)):
                                if best_hits[i][1].lower().find('random') < 0:
                                    best_hits.insert(0,best_hits.pop(i))
                                    break

                    if evalue2 > evalue1:
                        if best_hits[1][1].lower().find('chrun') > -1:
                            for i in range(2,len(best_hits)):
                                if best_hits[i][1].lower().find('chrun') < 0:
                                    best_hits.insert(1,best_hits.pop(i))
                                    break
                    else:
                        if best_hits[0][1].lower().find('chrun') > -1:
                            for i in range(1,len(best_hits)):
                                if best_hits[i][1].lower().find('chrun') < 0:
                                    best_hits.insert(0,best_hits.pop(i))
                                    break
                    
                    if evalue1==evalue2:
                        if evalues.count(evalue1) > 2:
                            result_code = 'multiple'
                        elif hits == 2:
                            result_code = 'tied'
                        else:
                            if len(best_hits) == 2:
                                evalue1,evalue2,evalue3 = heapq.nsmallest(3,evalues)
                                for line in hits_data:
                                    if float(line[10]) == evalue3:
                                        best_hits.append(line)
                            if (float(best_hits[2][2])/100*int(best_hits[2][3])/int(best_hits[2][12]))/(float(best_hits[0][2])/100*int(best_hits[0][3])/int(best_hits[0][12])) < 0.7:
                                result_code = 'tied+'
                            else:
                                result_code = 'multiple'
                    else:
                        if (float(best_hits[1][2])/100*int(best_hits[1][3])/int(best_hits[1][12]))/(float(best_hits[0][2])/100*int(best_hits[0][3])/int(best_hits[0][12])) < 0.7:
                            result_code = 'one+'
                        else:
                            result_code = 'best'

                if best_hits[0][8] > best_hits[0][9]: 
                    direction1=-1
                    approxstart1 = int(best_hits[0][8])+int(best_hits[0][6])-1
                else:
                    direction1=1
                    approxstart1 = int(best_hits[0][8])-(int(best_hits[0][6])-1)
                               

                if result_code=='one':

                    outfile.write(cluster+'\t'+str(hits)+'\t'+result_code+'\t'+best_hits[0][1]+'\t'+str(approxstart1)+'\t'+str(direction1)+'\t'+best_hits[0][12]+'\t'+best_hits[0][3]+'\t'+best_hits[0][2]+'\t'+best_hits[0][10]+'\t.\t.\t.\t.\t.\t.\t.\n')
                else:
                    if best_hits[1][8] > best_hits[1][9]:
                        direction2=-1
                        approxstart2 = int(best_hits[1][8])+int(best_hits[1][6])-1
                    else:
                        direction2=1
                        approxstart2 = int(best_hits[1][8])-(int(best_hits[1][6])-1)

                    outfile.write(cluster+'\t'+str(hits)+'\t'+result_code+'\t'+best_hits[0][1]+'\t'+str(approxstart1)+'\t'+str(direction1)+'\t'+best_hits[0][12]+'\t'+best_hits[0][3]+'\t'+best_hits[0][2]+'\t'+best_hits[0][10]+'\t'+best_hits[1][1]+'\t'+str(approxstart2)+'\t'+str(direction2)+'\t'+best_hits[1][12]+'\t'+best_hits[1][3]+'\t'+best_hits[1][2]+'\t'+best_hits[1][10]+'\n')
                    
                                   
outfile.close()
infile.close()
