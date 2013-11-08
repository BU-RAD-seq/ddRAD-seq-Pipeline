#######################################################################
##
## ParallelParser.py
##
## Version 1.00 -- 11 July 2013
##
## Created by Michael Sorenson
## Copyright (c) 2011-2013 Boston University. All rights reserved.
##
## Dependencies: specs file (see example.specs)
##               ddRADparser.py
##               CondenseSequences.py
##               FilterSequences.py
##
## This program is written for execution in the python (version 3)
## language, and uses the multiprocessing module and ssh to run
## parallel jobs in a Diskless Remote Boot Linux (DRBL) environment.
## It is free and distributed WITHOUT warranty; without
## even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.
##
#######################################################################

import sys
import os
import subprocess
import multiprocessing
from multiprocessing import Process, Lock

def start_run(cmd):
    p = subprocess.Popen(cmd, shell=True)
    sts = os.waitpid(p.pid, 0)[1]


RADspecs_file = open(sys.argv[1], 'r')

#get info from RADspecs file
line = RADspecs_file.readline()
data = line.split()
if data[1]=='2':
    Paired = True
else:
    Paired = False
    
line = RADspecs_file.readline()
data = line.split('\t')
path = data[1]

line = RADspecs_file.readline()

line = RADspecs_file.readline()
data = line.split()
basefilename1 = data[1]
qseq_file1 = open(path+data[1]+'.qseq','r')

line = RADspecs_file.readline()
if Paired == True:
    data = line.split()
    basefilename2 = data[1]
    qseq_file2 = open(path+data[1]+'.qseq','r')

line = RADspecs_file.readline()
data = line.split()
barcode_file = open(path+data[1],'r')

RADspecs_file.readline()
RADspecs_file.readline()
RADspecs_file.readline()

line = RADspecs_file.readline()
data = line.split()
r_site1_complete_length = len(data[1])

for i in range(5):
    line = RADspecs_file.readline()

data = line.split()
nruns = int(data[1])
print("N_processes: ",nruns)

line = RADspecs_file.readline()
hosts = []
for i in range(nruns):
    line = RADspecs_file.readline()
    data = line.split()
    hosts.append(data[0])

hosts3 = hosts*3

#length of input files
cmd = ('wc -l '+path+basefilename1+'.qseq')
p = subprocess.getoutput(cmd)
file_data = p.split()
file_length = int(file_data[0])
target_length = round((file_length/(nruns*3))+1)

if Paired == True:
    #create input files - including multiple specs files
    line_count = 0
    file_count = 101
    target = target_length
    outfile1 = open(path+basefilename1+str(file_count)+'.qseq','w')
    outfile2 = open(path+basefilename2+str(file_count)+'.qseq','w')
    for line1 in qseq_file1:
        line2 = qseq_file2.readline()
        line_count += 1
        outfile1.write(line1)
        outfile2.write(line2)
        if line_count == target:
            outfile1.close()
            outfile2.close()
            cmd = ('cp '+sys.argv[1]+' temp')
            p = subprocess.Popen(cmd, shell=True)        
            sts = os.waitpid(p.pid, 0)[1]
            cmd = ('sed s/'+basefilename1+'/'+basefilename1+str(file_count)+'/ <temp >temp2')
            p = subprocess.Popen(cmd, shell=True)        
            sts = os.waitpid(p.pid, 0)[1]
            cmd = ('sed s/'+basefilename2+'/'+basefilename2+str(file_count)+'/ <temp2 >'+path+sys.argv[1]+str(file_count))
            p = subprocess.Popen(cmd, shell=True)        
            sts = os.waitpid(p.pid, 0)[1]
            file_count+=1
            target += target_length
            outfile1 = open(path+basefilename1+str(file_count)+'.qseq','w')
            outfile2 = open(path+basefilename2+str(file_count)+'.qseq','w')
    outfile1.close()
    outfile2.close()
    cmd = ('cp '+sys.argv[1]+' temp')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('sed s/'+basefilename1+'/'+basefilename1+str(file_count)+'/ <temp >temp2')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('sed s/'+basefilename2+'/'+basefilename2+str(file_count)+'/ <temp2 >'+path+sys.argv[1]+str(file_count))
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm temp temp2')
    p = subprocess.Popen(cmd, shell=True) 
    sts = os.waitpid(p.pid, 0)[1]


    #launch runs
    threads=[]
    if __name__ == '__main__':
        lock = Lock()
        for i in range(nruns*3):
            if len(threads) < nruns:
                cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'ddRADparser.py '+path+sys.argv[1]+str(101+i)+' '+str(101+i))
                print(cmd)
                p = multiprocessing.Process(target=start_run,args=[cmd])
                p.start()
                threads.append(p)
            else:
                while len(threads) == nruns:
                    tnum = -1
                    for thread in threads:
                        tnum += 1
                        if not thread.is_alive():
                            threads.remove(thread)
                            hosts.append(hosts.pop(tnum))

                cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'ddRADparser.py '+path+sys.argv[1]+str(101+i)+' '+str(101+i))
                print(cmd)

                p = multiprocessing.Process(target=start_run,args=[cmd])
                p.start()
                threads.append(p)
                
        while len(threads) > 0:
            for thread in threads:
                if not thread.is_alive():
                    threads.remove(thread)

    #concatenate files
    samples=[]
    for line in barcode_file:
        data = line.split()
        samples.append(data[2])
    num_samples = len(samples)
    for i in range(num_samples):
        cmd = ('cat '+path+samples[i]+'R11*qseq > '+path+samples[i]+'R1.qseq')
        p = subprocess.Popen(cmd, shell=True)        
        sts = os.waitpid(p.pid, 0)[1]
        cmd = ('rm '+path+samples[i]+'R11*qseq')
        p = subprocess.Popen(cmd, shell=True)        
        sts = os.waitpid(p.pid, 0)[1]
        cmd = ('cat '+path+samples[i]+'R21*qseq > '+path+samples[i]+'R2.qseq')
        p = subprocess.Popen(cmd, shell=True)        
        sts = os.waitpid(p.pid, 0)[1]
        cmd = ('rm '+path+samples[i]+'R21*qseq')
        p = subprocess.Popen(cmd, shell=True)        
        sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'TinyR1seqs1* > '+path+'TinyR1seqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'TinyR1seqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'TinyR2seqs1* > '+path+'TinyR2seqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'TinyR2seqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'UnassignedR1seqs1* > '+path+'UnassignedR1seqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'UnassignedR1seqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'UnassignedR2seqs1* > '+path+'UnassignedR2seqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'UnassignedR2seqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    
    cmd = ('cat '+path+'ConcatamerR1seqs1* > '+path+'ConcatamerR1seqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'ConcatamerR1seqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'ConcatamerR2seqs1* > '+path+'ConcatamerR2seqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'ConcatamerR2seqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('rm '+path+basefilename1+'1*qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('rm '+path+basefilename2+'1*qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('rm '+path+sys.argv[1]+'1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    summ_data = []
    inputfile = open(path+basefilename1+str(101)+'summary.out')
    for line in inputfile:
        data = line.strip('\n').split('\t')
        summ_data.append(data)

    for i in range(nruns-1):
        inputfile = open(path+basefilename1+str(102+i)+'summary.out')
        inputfile.readline()
        line_count=0
        for line in inputfile:
            line_count+=1
            data = line.strip('\n').split('\t')
            for j in range(1,len(data)):
                summ_data[line_count][j] = int(summ_data[line_count][j])+int(data[j])

    outputfile = open(path+basefilename1+'summary.out','w')
    for i in range(len(summ_data)):
        outputfile.write(str(summ_data[i][0]))
        for j in range(1,len(summ_data[i])):
            outputfile.write('\t'+str(summ_data[i][j]))
        outputfile.write('\n')
    outputfile.close()

    cmd = ('rm '+path+basefilename1+'1*summary*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    bstart = []
    bstart_counts= []
    for i in range(nruns):
        inputfile = open(path+'Unassigned_start_seqs'+str(101+i)+'.out')
        inputfile.readline()
        for line in inputfile:
            data = line.strip('\n').split('\t')
            if data[0] in bstart:
                bstart_counts[bstart.index(data[0])] += int(data[1])
            else:
                bstart.append(data[0])
                bstart_counts.append(int(data[1]))

    outputfile = open(path+'Unassigned_start_seqs.out','w')
    outputfile.write('Sequence\tCount\n')
    for i in range(len(bstart)):
        outputfile.write(bstart[i]+'\t'+str(bstart_counts[i])+'\n')
    outputfile.close()

    cmd = ('rm '+path+'Unassigned_start_seqs1*.out')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
            
##
## starts here for single read processing
##

else:
    #create input files - including multiple specs files
    line_count = 0
    file_count = 101
    print("Generating Files")
    print(file_count)
    target = target_length
    outfile1 = open(path+basefilename1+str(file_count)+'.qseq','w')
    for line1 in qseq_file1:
        line_count += 1
        outfile1.write(line1)
        if line_count == target:
            outfile1.close()
            cmd = ('cp '+sys.argv[1]+' temp')
            p = subprocess.Popen(cmd, shell=True)        
            sts = os.waitpid(p.pid, 0)[1]
            cmd = ('sed s/'+basefilename1+'/'+basefilename1+str(file_count)+'/ <temp >'+path+sys.argv[1]+str(file_count))
            p = subprocess.Popen(cmd, shell=True)        
            sts = os.waitpid(p.pid, 0)[1]
            file_count+=1
            target += target_length
            outfile1 = open(path+basefilename1+str(file_count)+'.qseq','w')
            print(file_count)
    outfile1.close()
    cmd = ('cp '+sys.argv[1]+' temp')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('sed s/'+basefilename1+'/'+basefilename1+str(file_count)+'/ <temp >'+path+sys.argv[1]+str(file_count))
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm temp')
    p = subprocess.Popen(cmd, shell=True) 
    sts = os.waitpid(p.pid, 0)[1]
   

    #launch runs
    threads=[]
    if __name__ == '__main__':
        lock = Lock()
        for i in range(nruns*3):
            if len(threads) < nruns:
                cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'ddRADparser.py '+path+sys.argv[1]+str(101+i)+' '+str(101+i))
                p = multiprocessing.Process(target=start_run,args=[cmd])
                p.start()
                threads.append(p)
            else:
                while len(threads) == nruns:
                    tnum = -1
                    for thread in threads:
                        tnum += 1
                        if not thread.is_alive():
                            threads.remove(thread)
                            hosts.append(hosts.pop(tnum))

                cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'ddRADparser.py '+path+sys.argv[1]+str(101+i)+' '+str(101+i))
                p = multiprocessing.Process(target=start_run,args=[cmd])
                p.start()
                threads.append(p)
            
        while len(threads) > 0:
            for thread in threads:
                if not thread.is_alive():
                    threads.remove(thread)

    #concatenate files
    samples=[]
    for line in barcode_file:
        data = line.split()
        samples.append(data[2])
    num_samples = len(samples)
    for i in range(num_samples):
        cmd = ('cat '+path+samples[i]+'1*qseq > '+path+samples[i]+'.qseq')
        p = subprocess.Popen(cmd, shell=True)        
        sts = os.waitpid(p.pid, 0)[1]
        cmd = ('rm '+path+samples[i]+'1*qseq')
        p = subprocess.Popen(cmd, shell=True)        
        sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'TinySeqs1* > '+path+'TinySeqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'TinySeqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'UnassignedSeqs1* > '+path+'UnassignedSeqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'UnassignedSeqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('cat '+path+'ConcatamerSeqs1* > '+path+'ConcatamerSeqs.qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]
    cmd = ('rm '+path+'ConcatamerSeqs1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('rm '+path+basefilename1+'1*qseq')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    cmd = ('rm '+path+sys.argv[1]+'1*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    summ_data = []
    inputfile = open(path+basefilename1+str(101)+'summary.out')
    for line in inputfile:
        data = line.strip('\n').split('\t')
        summ_data.append(data)

    for i in range(nruns*3-1):
        inputfile = open(path+basefilename1+str(102+i)+'summary.out')
        inputfile.readline()
        line_count=0
        for line in inputfile:
            line_count+=1
            data = line.strip('\n').split('\t')
            for j in range(1,len(data)):
                summ_data[line_count][j] = int(summ_data[line_count][j])+int(data[j])

    outputfile = open(path+basefilename1+'summary.out','w')
    for i in range(len(summ_data)):
        outputfile.write(str(summ_data[i][0]))
        for j in range(1,len(summ_data[i])):
            outputfile.write('\t'+str(summ_data[i][j]))
        outputfile.write('\n')
    outputfile.close()

    cmd = ('rm '+path+basefilename1+'1*summary*')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    bstart = []
    bstart_counts= []
    for i in range(nruns*3):
        inputfile = open(path+'Unassigned_start_seqs'+str(101+i)+'.out')
        inputfile.readline()
        for line in inputfile:
            data = line.strip('\n').split('\t')
            if data[0] in bstart:
                bstart_counts[bstart.index(data[0])] += int(data[1])
            else:
                bstart.append(data[0])
                bstart_counts.append(int(data[1]))

    outputfile = open(path+'Unassigned_start_seqs.out','w')
    outputfile.write('Sequence\tCount\n')
    for i in range(len(bstart)):
        outputfile.write(bstart[i]+'\t'+str(bstart_counts[i])+'\n')
    outputfile.close()

    cmd = ('rm '+path+'Unassigned_start_seqs1*.out')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    #run CondenseSequences.py
    threads=[]
    for i in range(num_samples):
        if len(threads) < nruns:
            cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'CondenseSequences.py '+path+samples[i]+'.qseq '+str(r_site1_complete_length))
            p = multiprocessing.Process(target=start_run,args=[cmd])
            p.start()
            threads.append(p)
        
        else:
            while len(threads) == nruns:
                tnum = -1
                for thread in threads:
                    tnum += 1
                    if not thread.is_alive():
                        threads.remove(thread)
                        hosts.append(hosts.pop(tnum))

            cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'CondenseSequences.py '+path+samples[i]+'.qseq '+str(r_site1_complete_length))
            p = multiprocessing.Process(target=start_run,args=[cmd])
            p.start()
            threads.append(p)

                        
    while len(threads) > 0:
        for thread in threads:
            if not thread.is_alive():
                threads.remove(thread)

    inputfile = open(path+basefilename1+'summary.out','r')
    line = inputfile.readline()
    outputfile = open(path+basefilename1+'sumtemp.out','w')
    outputfile.write(line)
    for i in range(num_samples):
        inputfile2 = open(path+samples[i]+'.sumtemp','r')
        line2=inputfile2.readline()
        data2 = line2.split()
        line = inputfile.readline()
        data = line.split()
        data.insert(6,data2[2])
        outputfile.write('\t'.join(data)+'\n')
        cmd = ('rm '+path+samples[i]+'.sumtemp')
        p = subprocess.Popen(cmd, shell=True)        
        sts = os.waitpid(p.pid, 0)[1]

    for i in range(4):
        line = inputfile.readline()
        outputfile.write(line)

    outputfile.close()
    
    cmd = ('mv '+path+basefilename1+'sumtemp.out '+path+basefilename1+'summary.out')
    p = subprocess.Popen(cmd, shell=True)        
    sts = os.waitpid(p.pid, 0)[1]

    #run FilterSequences.py
    threads=[]
    for i in range(num_samples):
        if len(threads) < nruns:
            cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'FilterSequences.py '+path+samples[i]+'C.qseq')
            p = multiprocessing.Process(target=start_run,args=[cmd])
            p.start()
            threads.append(p)
        
        else:
            while len(threads) == nruns:
                tnum = -1
                for thread in threads:
                    tnum += 1
                    if not thread.is_alive():
                        threads.remove(thread)
                        hosts.append(hosts.pop(tnum))

            cmd = ('ssh '+hosts[len(threads)]+' python3 '+path+'FilterSequences.py '+path+samples[i]+'C.qseq')
            p = multiprocessing.Process(target=start_run,args=[cmd])
            p.start()
            threads.append(p)

                    
    while len(threads) > 0:
        for thread in threads:
            if not thread.is_alive():
                threads.remove(thread)
