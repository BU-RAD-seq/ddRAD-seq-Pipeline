#######################################################################
##
## ddRADparser.py
##
## Version 1.04 -- 14 September 2013
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
import os, math, subprocess, heapq

def reverse_complement(seq):
    seq = list(seq)
    seq.reverse()
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    complseq = [complement[base] for base in seq]
    return (''.join(complseq))

def Align(s1,s2): ###Corrected code inserted 14 Sept 2013, version 1.04###
    M=[]
    S=[]
    for i in range(len(s1)+1):
        M.append([1])
        S.append([i])
    #initialize first row of M matrix
    M[0] = [0]+[1 for x in range(len(s2))]
    #calculate M matrix
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i]==s2[j]:
                M[i+1].append(0)
            else:
                M[i+1].append(1)
    #initialize first row of S matrix
    S[0] = [x for x in range(len(s2)+1)]    
    #calculate rest of S matrix
    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):
            score = S[i-1][j-1]+M[i][j]
            if S[i][j-1]+1 < score:
                score = S[i][j-1]+1
            if S[i-1][j]+1 < score:
                score = S[i-1][j]+1
            S[i].append(score)            
    return (S)

def FindRsite(Smatrix):
    Rrow=Smatrix[-len(r_site1)]
    while Rrow.count(min(Rrow)) > 1:
        Rrow[Rrow.index(min(Rrow))]=min(Rrow)+1
    cutpoint = Rrow.index(min(Rrow))
    return(cutpoint)

def find_P2(inseq,matchseq,len_rsite):
    found = False
    code = 0
    location = -1
    if inseq.find(matchseq[0:17]) > 0:
        location = inseq.find(matchseq[0:17])
        found=True
        code=1
    else:
        locations = []
        i = inseq.find(matchseq[0:len_rsite])
        while i != -1:
            locations.append(i)
            i = inseq.find(matchseq[0:len_rsite],i+1)
        i = inseq.find(matchseq[len_rsite:len_rsite+5])
        while i != -1:
            locations.append(i-len_rsite)
            i = inseq.find(matchseq[len_rsite:len_rsite+5],i+1)
        locations = list(set(locations))

        for i in locations:            
            test_seq = inseq[i:i+21]
            if len(test_seq) > len_rsite:
                if test_seq == matchseq[:len(test_seq)]:
                    found=True
                    code=2
                    location = i
                    break
                else:     
                    if sum(c1 != c2 for c1,c2 in zip(test_seq,matchseq[:len(test_seq)]))/(len(test_seq)-5) < 0.3:
                        found=True
                        code=3
                        location = i
                        break
                    else:
                        AlignResult=Align(test_seq,matchseq[:len(test_seq)])
                        if AlignResult[-1][-1]/(len(test_seq)-5) < 0.3:
                            found=True
                            code=4
                            location = i
                            break
    return(found,location,code)

Nfix = 0
rsitefix1 = 0
rsitefix2 = 0

RADspecs_file = open(sys.argv[1],'r')

while len(sys.argv) < 3:
    sys.argv.append('')
    
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
data = line.split()
read_length = int(data[1])

line = RADspecs_file.readline()
data = line.split()
basefilename1 = data[1]

line = RADspecs_file.readline()
if Paired == True:
    data = line.split()
    basefilename2 = data[1]

line = RADspecs_file.readline()
data = line.split()
barcode_file = open(path+data[1],'r')

line = RADspecs_file.readline()
data = line.split()
if data[1]=='True':
    BCcorrection = True
else:
    BCcorrection = False

line = RADspecs_file.readline()
data = line.split()
if data[1]=='True':
    adjust_quals = True
else:
    adjust_quals = False
        
line = RADspecs_file.readline()
data = line.split()
r_site1 = data[1]

line = RADspecs_file.readline()
data = line.split()
r_site1_complete = data[1]

line = RADspecs_file.readline()
data = line.split()
r_site2 = data[1]

line = RADspecs_file.readline()
data = line.split()
r_site2_complete = data[1]

line = RADspecs_file.readline()
data = line.split()
P1adapter = data[1]

line = RADspecs_file.readline()
data = line.split()
P2adapter = data[1]

matchseq2 = r_site2+P2adapter
P2_codes = [0]*5

failed_align1 = 0
successful_align1 = 0
failed_align2 = 0
successful_align2 = 0

###########################################
## starts here for paired end processing ##
###########################################

if Paired == True:
    qseq_file1 = open(path+basefilename1+'.qseq','r')
    qseq_file2 = open(path+basefilename2+'.qseq','r')
    barcodes=[]
    barcodesBC=[]
    barcodesBCA=[]
    samples=[]
    for line in barcode_file:
        data = line.split()
        barcodes.append(data[1]+r_site1[:2])
        barcodesBC.append(data[1]+r_site1[:3])
        barcodesBCA.append(data[1]+r_site1)
        samples.append(data[2])
        vars()['filename1' + str(len(barcodes))] = open(path+data[2] + 'R1'+sys.argv[2]+'.qseq','w')
        vars()['filename2' + str(len(barcodes))] = open(path+data[2] + 'R2'+sys.argv[2]+'.qseq','w')
    num_samples = len(samples)

    r_site1_length = len(r_site1)
    r_site2_length = len(r_site2)
    r_site1_complete_length = len(r_site1_complete)
    r_site2_complete_length = len(r_site2_complete)
    barcode_length = len(barcodes[0])-2
    test_length = len(barcodes[0])
    test_lengthBC = len(barcodesBC[0])
    test_lengthBCA = len(barcodesBCA[0])
    
    R1addseq = r_site1_complete.replace(r_site1,'')
    revcomp_R1addseq = reverse_complement(R1addseq)
    R2addseq = r_site2_complete.replace(r_site2,'')
    revcomp_R2addseq = reverse_complement(R2addseq)
    R2addqual = 'H' * len(R2addseq)
    len_R1add = len(R1addseq)
    len_R2add = len(R2addseq)
    R1R2add = len_R1add+len_R2add
    
    R1length = read_length - barcode_length + len_R1add
    R2length = read_length + len_R2add

    S1Locations1=[0] * (R1length + 4)
    S2Locations1=[0] * (R1length + 4)
    S1Locations2=[0] * (R2length + 4)
    S2Locations2=[0] * (R2length + 4)
    STagLocations=[0] * (R1length + 4)

    seq_counts = [0] * num_samples
    tiny_seqs = [0] * num_samples
    concatamer1R1 = [0] * num_samples
    concatamer1R2 = [0] * num_samples
    concatamer2R1 = [0] * num_samples
    concatamer2R2 = [0] * num_samples
    short_tags=[]
    for i in range(num_samples):
        short_tags.append(STagLocations[0:R1length+3])

    unassigned_fileR1 = open(path+'UnassignedR1seqs'+sys.argv[2]+'.qseq', 'w')
    unassigned_fileR2 = open(path+'UnassignedR2seqs'+sys.argv[2]+'.qseq', 'w')
    ConcatamersR1 = open(path+'ConcatamerR1seqs'+sys.argv[2]+'.qseq', 'w')
    ConcatamersR2 = open(path+'ConcatamerR2seqs'+sys.argv[2]+'.qseq', 'w')
    tiny_fileR1 = open(path+'TinyR1seqs'+sys.argv[2]+'.qseq', 'w')
    tiny_fileR2 = open(path+'TinyR2seqs'+sys.argv[2]+'.qseq', 'w')
    start_listB=[]
    start_countB=[]
    bstart = open(path+'Unassigned_start_seqs'+sys.argv[2]+'.out', 'w')

    no_barcode=0

    line_count=0
    target = 10000
    for line1 in qseq_file1:
        line2 = qseq_file2.readline()
        line_count += 1
        if line_count == target:
            print(target)
            target += 10000
        data1 = line1.split()
        data2 = line2.split()
            
        sample_index = -1    
        #fix any single N in restriction site if adjacent bases are basically correct
        teststring = data1[8][barcode_length:test_length]
        if teststring.count("N") == 1:
            if sum(c1 != c2 for c1,c2 in zip(teststring,r_site1)) < 3:
                N_pos = teststring.find("N")
                data1[8] = data1[8][:barcode_length+N_pos]+r_site1[N_pos]+data1[8][barcode_length+N_pos+1:]
                Nfix+=1
            
        #fix any remaining single base error in overhang of restriction site
        if teststring[:4] != r_site1[:4]:
            if sum(c1 != c2 for c1,c2 in zip(teststring[:4],r_site1[:4])) == 1:
                data1[8] = data1[8][:barcode_length]+r_site1[:4]+data1[8][barcode_length+4:]
                rsitefix1+=1
        
        test_code = data1[8][:barcode_length]
        #check for r_site in correct position
        if sum(c1 != c2 for c1,c2 in zip(teststring[:4],r_site1[:4])) < 2:    
            # identify barcode or find best match to barcode
            if test_code in barcodes:
                sample_index = barcodes.index(test_code)
            elif BCcorrection == True:
                cost=[]
                for i in barcodes:
                    cost.append(sum(c1 != c2 for c1,c2 in zip(i,test_code)))
                min_cost = heapq.nsmallest(2,cost)
                if min_cost[0] < 2 and min_cost[1]-min_cost[0] > 0:
                    sample_index = cost.index(min_cost[0])
                    successful_align1 += 1
                else:
                    failed_align1 += 1
                                                        
        if sample_index < 0:
            unassigned_fileR1.write(line1)
            unassigned_fileR2.write(line2)
            no_barcode+=1
            if test_code in start_listB:
                start_countB[start_listB.index(test_code)]+=1
            else:
                start_listB.append(test_code)
                start_countB.append(1)

        else:
            seq_counts[sample_index] += 1
            data1[8] = data1[8].replace(data1[8][0:r_site_position], R1addseq, 1)
            data1[9] = data1[9].replace(data1[9][0:r_site_position-len_R1add], '', 1)

            # check for concatamers and short inserts
            P2found=False
            Concatamer=False
            # check for r_site1 concatamer
            position = data1[8].find(r_site1_complete,2)+r_site1_complete_length
            if position > r_site1_complete_length:    
                ConcatamersR1.write(str(position)+'\t'+line1)
                ConcatamersR2.write(line2)
                concatamer1R1[sample_index] += 1
                S1Locations1[position]+=1
            else:
                # check for short insert or r_site2 concatamer
                P2found,position,findcode = find_P2(data1[8],matchseq2,r_site2_length)
                P2_codes[findcode]+=1
                    
                if P2found==True:
                    data1[8]=data1[8][:position]+r_site2_complete
                    data1[9]=data1[9][:len(data1[8])]
                    short_tags[sample_index][len(data1[8])] += 1
                    data2[8]=data2[8][:len(data1[8])-R1R2add]+revcomp_R1addseq
                else:
                    position = data1[8].find(r_site2_complete)+r_site2_complete_length
                    if position > r_site2_complete_length:
                        concatamer2R1[sample_index] += 1
                        S2Locations1[position]+=1
                        ConcatamersR1.write(str(position)+'\t'+line1)
                        ConcatamersR2.write(line2)
                        Concatamer=True                    
                    else:
                        position = data2[8].find(r_site1_complete)+r_site1_complete_length
                        if position > r_site1_complete_length:
                            ConcatamersR1.write(line1)
                            ConcatamersR2.write(str(position)+'\t'+line2)
                            concatamer1R2[sample_index] += 1
                            S1Locations2[position]+=1
                            Concatamer=True
                        else:
                            position = data2[8].find(r_site2_complete) + r_site2_complete_length
                            if position > r_site2_complete_length:
                                concatamer2R2[sample_index] += 1
                                S2Locations2[position]+=1
                                ConcatamersR1.write(line1)
                                ConcatamersR2.write(str(position)+'\t'+line2)
                                Concatamer=True
                        
                if Concatamer == False:
                    if len(data1[8]) > 31:
                        #if 'N' in data1[8][:8]:
                        #    data1[8] = fix_rsite(data1[8])
                        #recode qualities
                        if adjust_quals == True:
                            new_qual1=list(data1[9])
                            for i in range(len(data1[9])):
                                new_qual1[i] = chr(round(10*math.log10(1+10**((ord(data1[9][i])-64)/10)))+33)
                            data1[9]= "".join(new_qual1)
    
                            new_qual2=list(data2[9])
                            for i in range(len(data2[9])):
                                new_qual2[i] = chr(round(10*math.log10(1+10**((ord(data2[9][i])-64)/10)))+33)
                            data2[9]= "".join(new_qual1)
    

                        data2[8]=revcomp_R2addseq+data2[8]
                        data2[9]=R2addqual+data2[9]
                        if P2found == True:
                            data2[9]=data2[9][:len(data2[8])]
                        outfile_name = 'filename1' + str(sample_index + 1)
                        vars()[outfile_name].write(samples[sample_index])
                        for i in range(1,11):
                            vars()[outfile_name].write('\t'+data1[i])
                        vars()[outfile_name].write('\n')
                            
                        outfile_name = 'filename2' + str(sample_index + 1)
                        vars()[outfile_name].write(samples[sample_index])
                        for i in range(1,11):
                            vars()[outfile_name].write('\t'+data2[i])
                        vars()[outfile_name].write('\n')   

                    else:
                        tiny_fileR1.write(line1)
                        tiny_fileR2.write(line2)
                        tiny_seqs[sample_index] +=1    
                
    for i in range(num_samples):
        vars()['filename1' + str(i + 1)].close()
        vars()['filename2' + str(i + 1)].close()
    unassigned_fileR1.close()
    unassigned_fileR2.close()
    ConcatamersR1.close()
    ConcatamersR2.close()
    tiny_fileR1.close()
    tiny_fileR2.close()

    bstart.write('Sequence\tCount\n')
    for i in range(len(start_listB)):
        bstart.write(start_listB[i]+'\t'+str(start_countB[i])+'\n')
    bstart.close()    

    for i in range(len(barcodes)):
        print(samples[i]+'\t'+str(seq_counts[i]))

    sum_stat_file=open(path+basefilename1+'summary.out','w')
    sum_stat_file.write('SampleID\tSeq_count\tTiny_seqs\tRS1.1_cat\tRS2.1_cat\tRS1.2_cat\tRS2.2_cat\tShort_tags\tUniqueSeqs\tShort_tag_locations\n')
    for i in range(num_samples):
        sum_stat_file.write(samples[i]+'\t'+str(seq_counts[i])+'\t'+str(tiny_seqs[i])+'\t'+str(concatamer1R1[i])+'\t'
                            +str(concatamer2R1[i])+'\t'+str(concatamer1R2[i])+'\t'+str(concatamer2R2[i])+'\t'+str(sum(short_tags[i])))
        for j in range(len(short_tags[i])):
            sum_stat_file.write('\t'+str(short_tags[i][j]))
        sum_stat_file.write('\n')
    sum_stat_file.write('Unassigned\t'+str(no_barcode)+'\n')
    sum_stat_file.write('\nConcatamer locations summed across all samples\n')
    sum_stat_file.write('Rsite1 concatamers in read 1')
    for i in range(len(S1Locations1)):
        sum_stat_file.write('\t'+str(S1Locations1[i]))
    sum_stat_file.write('\nRsite2 concatamers in read 1')
    for i in range(len(S2Locations1)):
        sum_stat_file.write('\t'+str(S2Locations1[i]))
    sum_stat_file.write('\nRsite1 concatamers in read 2')
    for i in range(len(S1Locations2)):
        sum_stat_file.write('\t'+str(S1Locations2[i]))
    sum_stat_file.write('\nRsite2 concatamers in read 2')
    for i in range(len(S2Locations2)):
        sum_stat_file.write('\t'+str(S2Locations2[i]))
    sum_stat_file.write('\n')    
    sum_stat_file.close()

############################################
## starts here for single read processing ##
############################################
    
else:
    qseq_file1 = open(path+basefilename1+'.qseq','r')
        
    barcodes=[]
    samples=[]
    for line in barcode_file:
        data = line.split()
        barcodes.append(data[1])
        samples.append(data[2])
        vars()['filename1' + str(len(barcodes))] = open(path+data[2]+sys.argv[2]+'.qseq','w')
    num_samples = len(samples)

    r_site1_length = len(r_site1)
    r_site2_length = len(r_site2)
    r_site1_complete_length = len(r_site1_complete)
    r_site2_complete_length = len(r_site2_complete)
    barcode_length = len(barcodes[0]) #-2
    test_length = barcode_length+r_site1_length

    R1addseq = r_site1_complete.replace(r_site1,'')
    revcomp_R1addseq = reverse_complement(R1addseq)
    len_R1add = len(R1addseq)
    
    R1length = read_length - barcode_length + len_R1add

    S1Locations1=[0] * (R1length + 4)
    S2Locations1=[0] * (R1length + 4)
    STagLocations=[0] * (R1length + 4)

    seq_counts = [0] * num_samples
    corrected_codes = [0] * num_samples
    tiny_seqs = [0] * num_samples
    concatamer1R1 = [0] * num_samples
    concatamer2R1 = [0] * num_samples
    short_tags=[]
    for i in range(num_samples):
        short_tags.append(STagLocations[0:R1length+3])

    unassigned_fileR1 = open(path+'UnassignedSeqs'+sys.argv[2]+'.qseq', 'w')
    ConcatamersR1 = open(path+'ConcatamerSeqs'+sys.argv[2]+'.qseq', 'w')
    tiny_fileR1 = open(path+'TinySeqs'+sys.argv[2]+'.qseq', 'w')
    start_listB=[]
    start_countB=[]
    bstart = open(path+'Unassigned_start_seqs'+sys.argv[2]+'.out', 'w')

    no_barcode=0

    line_count=0
    target = 10000
    for line1 in qseq_file1:
        line_count += 1
        if line_count == target:
            print(target)
            target += 10000
        data1 = line1.split()
        
        sample_index = -1    
        #fix any single N in restriction site if adjacent bases are basically correct
        teststring = data1[8][barcode_length:test_length]
        if teststring.count("N") == 1:
            if sum(c1 != c2 for c1,c2 in zip(teststring,r_site1)) < 3:
                N_pos = teststring.find("N")
                data1[8] = data1[8][:barcode_length+N_pos]+r_site1[N_pos]+data1[8][barcode_length+N_pos+1:]
                Nfix+=1
            
        #fix any remaining single base error in overhang of restriction site
        if teststring[:4] != r_site1[:4]:
            if sum(c1 != c2 for c1,c2 in zip(teststring[:4],r_site1[:4])) == 1:
                data1[8] = data1[8][:barcode_length]+r_site1[:4]+data1[8][barcode_length+4:]
                rsitefix1+=1
        
        test_code = data1[8][:barcode_length]
        #check for r_site in correct position
        if sum(c1 != c2 for c1,c2 in zip(teststring[:4],r_site1[:4])) < 2:    
            # identify barcode or find best match to barcode
            if test_code in barcodes:
                sample_index = barcodes.index(test_code)
            elif BCcorrection == True:
                cost=[]
                for i in barcodes:
                    cost.append(sum(c1 != c2 for c1,c2 in zip(i,test_code)))
                min_cost = heapq.nsmallest(2,cost)
                if min_cost[0] < 2 and min_cost[1]-min_cost[0] > 0:
                    sample_index = cost.index(min_cost[0])
                    corrected_codes[sample_index] += 1
                    successful_align1 += 1
                    data1[7] = '3'
                else:
                    failed_align1 += 1
                                
        if sample_index < 0:
            unassigned_fileR1.write(line1)
            no_barcode+=1
            if test_code in start_listB:
                start_countB[start_listB.index(test_code)]+=1
            else:
                start_listB.append(test_code)
                start_countB.append(1)
        else:
            seq_counts[sample_index] += 1
            data1[8] = R1addseq + data1[8][barcode_length:]
            data1[9] = data1[9][barcode_length-len_R1add:]
                            
            # check for short tags

            #check for r_site2 concatamer
            position = data1[8].find(r_site2_complete)+r_site2_complete_length
            if position > r_site2_complete_length:
                concatamer2R1[sample_index] += 1
                S2Locations1[position]+=1
                data1[8]=data1[8][:position]
                data1[9]=data1[9][:len(data1[8])]
                
            # check for r_site1 concatamer
            position = data1[8].find(r_site1_complete,r_site1_length)+r_site1_complete_length
            if position > r_site1_complete_length:    
                data1[8]=data1[8][:position]
                ConcatamersR1.write(str(position)+'\t'+line1)
                concatamer1R1[sample_index] += 1
                S1Locations1[position]+=1
            else:
                # check for short insert
                P2found,position,findcode = find_P2(data1[8],matchseq2,r_site2_length)
                P2_codes[findcode]+=1
                
                if P2found==True:
                    data1[8]=data1[8][:position]+r_site2_complete
                    data1[9]=data1[9][:len(data1[8])]
                    short_tags[sample_index][len(data1[8])] += 1
                                                    
                if len(data1[8]) > 31:
                    if adjust_quals == True:
                        new_qual1=list(data1[9])
                        for i in range(len(data1[9])):
                            new_qual1[i] = chr(round(10*math.log10(1+10**((ord(data1[9][i])-64)/10)))+33)
                        data1[9]= "".join(new_qual1)
    
                    outfile_name = 'filename1' + str(sample_index + 1)
                    vars()[outfile_name].write(samples[sample_index])
                    for i in range(1,11):
                        vars()[outfile_name].write('\t'+data1[i])
                    vars()[outfile_name].write('\n')
                        
                else:
                    tiny_fileR1.write(line1)
                    tiny_seqs[sample_index] +=1    
                
    for i in range(num_samples):
        vars()['filename1' + str(i + 1)].close()
    unassigned_fileR1.close()
    ConcatamersR1.close()
    tiny_fileR1.close()

    bstart.write('Sequence\tCount\n')
    for i in range(len(start_listB)):
        bstart.write(start_listB[i]+'\t'+str(start_countB[i])+'\n')
    bstart.close()    

    for i in range(len(barcodes)):
        print(samples[i]+'\t'+str(seq_counts[i]))

    sum_stat_file=open(path+basefilename1+'summary.out','w')
    sum_stat_file.write('SampleID\tSeq_count\tBarcodeCorrected\tTiny_seqs\tRS1.1_cat\tRS2.1_cat\tShort_tags\tUniqueSeqs\tShort_tag_locations\n')
    for i in range(num_samples):
        sum_stat_file.write(samples[i]+'\t'+str(seq_counts[i])+'\t'+str(corrected_codes[i])+'\t'+str(tiny_seqs[i])+'\t'+str(concatamer1R1[i])+'\t'
                            +str(concatamer2R1[i])+'\t'+str(sum(short_tags[i])))
        for j in range(len(short_tags[i])):
            sum_stat_file.write('\t'+str(short_tags[i][j]))
        sum_stat_file.write('\n')
    sum_stat_file.write('Unassigned\t'+str(no_barcode)+'\n')
    sum_stat_file.write('\nConcatamer locations summed across all samples\n')
    sum_stat_file.write('Rsite1 concatamers in read 1')
    for i in range(len(S1Locations1)):
        sum_stat_file.write('\t'+str(S1Locations1[i]))
    sum_stat_file.write('\nRsite2 concatamers in read 1')
    for i in range(len(S2Locations1)):
        sum_stat_file.write('\t'+str(S2Locations1[i]))
    sum_stat_file.write('\n')    
    sum_stat_file.close()


print("P2_codes: ",P2_codes)
print("Corrected barcodes (succeed/fail): ",successful_align1,"/",failed_align1)
print("Nfix: ",Nfix,"Rsitefix: ",rsitefix1)
