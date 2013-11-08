#######################################################################
##
## RADGenotypes.py
##
## Version 1.12 -- 30 October 2013
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

import sys, os, random, heapq, re
from math import factorial
#from decimal import *
#from scipy.stats import norm
#from scipy.stats import binom

def median(alist):
    srtd = sorted(alist)
    mid = int(len(alist)/2)  
    if len(alist) > 1:
        if len(alist) % 2 == 0:
            return (srtd[mid-1] + srtd[mid]) / 2.0
        else:
            return srtd[mid]
    else:
        return(srtd[0]) 
           
def binom_prob(n,k):
    combs = 0
    for i in range(k+1):    
        combs += factorial(n) / (factorial(i)*factorial(n-i))
    prob = 2*combs*0.5**n
    return prob

def binom_test(count1,count2):
    if count1 > 1:
        if float(count1)/(count1+count2) > 0.29:
            result = 1
        else:
            if count1 < 7: ##changed from "count1+count2 < 41" to remove the inconsistency of accepting 9,32 but not 9,30 or 9,31 as provisional heteros
                pvalue = binom_prob(count1+count2,count1)
                if pvalue > 0.05:
                    result = 1
                elif pvalue > 0.001:
                    result = 3
                else:
                    result = 4
            elif float(count1)/(count1+count2) > 0.19999:
                result = 3
            else:
                result = 4
    else:
        result = 3
    return(result)

def strip_alignment_gaps(con_seq,seq_array1,seq_array2,site_list):
    gap_sites = [m.start() for m in re.finditer('-', con_seq)]
    gap_sites.reverse()
    for i in range(len(gap_sites)):
        if gap_sites[i] not in site_list:
            for j in range(len(site_list)):
                if site_list[j] > gap_sites[i]:
                    site_list[j] -= 1
            con_seq = con_seq[:gap_sites[i]]+con_seq[gap_sites[i]+1:]
            for j in range(len(seq_array1)):
                if seq_array1[j][1] != '.':
                    seq_array1[j][1] = seq_array1[j][1][:gap_sites[i]]+seq_array1[j][1][gap_sites[i]+1:]
                if seq_array2[j][1] != '.':
                    seq_array2[j][1] = seq_array2[j][1][:gap_sites[i]]+seq_array2[j][1][gap_sites[i]+1:]
    return(con_seq,seq_array1,seq_array2,site_list)

#get input file name and open files
base_filename = sys.argv[1]
BLAST_depth = int(sys.argv[2])

input_file = open(base_filename.replace('.qseq','AS.qseq'),'r')
BLASTinput = open(base_filename.replace('.qseq','BLASTsummary.out'),'r')
barcode_file = open(base_filename.replace('.qseq','.index'), 'r')
sumfile = open(base_filename.replace('.qseq','clustersummary.out'),'w')

sumfile.write('Clstr\tHits\tChr\tPos\tDir\tCat\tLength\tPoly\tSNPS\tGaps\tHaps\tEdit\tInfSites\tHW-X2\thets+/-\tHetRatio\tHets\tNoData\tGood\tLowDepth\tBadRatio\tExtraReads\t3rdAllele\tConsensus')

#get sample names from sample list file
sample_list = []
for line in barcode_file:
    sample_data = line.strip('\n').split('\t')
    print(sample_data)
    sample_list.append(sample_data[2])
num_samples = len(sample_list)

for i in range(num_samples):
    sumfile.write('\t'+str(i+1))
for i in range(num_samples):
    sumfile.write('\t'+str(i+1))
sumfile.write('\n')

bigclusters = open('BigClusters','r')

vars()[base_filename + 'out']=open(base_filename.replace('.qseq','.out'),'w')

resultsCL = [[0]*6 for i in range(num_samples)]
resultsVL = [[0]*6 for i in range(num_samples)]
resultsVLX2 = [[0]*6 for i in range(num_samples)]
resultsVLXX = [[0]*6 for i in range(num_samples)]
resultsC = [[0]*6 for i in range(num_samples)]
resultsV = [[0]*6 for i in range(num_samples)]
resultsVX2 = [[0]*6 for i in range(num_samples)]
resultsVXX = [[0]*6 for i in range(num_samples)]

#create arrays for individual sample data
four_bases = [0,0,0,0,0,0]
ACGT=['A','C','G','T','-','N']
rsite = sys.argv[3]

RAD_sample_bases = []
RAD_sample_qual = []
weighted_sample_bases = []
sample_depth = [0]*num_samples
sample_weight = [0]*num_samples

Loci_recovered = [0]*(num_samples+1)
Loci_with_depth = [0]*(num_samples+1)
Constant_loci=[0]*(num_samples+1)
Passed_variable_loci=[0]*(num_samples+1)
Imperfect_variable_loci=[0]*(num_samples+1)
Total_depth=[0]*(num_samples+1)

for i in range(num_samples):
    RAD_sample_bases.append(four_bases)
    RAD_sample_qual.append(four_bases)
    weighted_sample_bases.append(four_bases)

failed = []

BLASTparts=[-1]

#read first line from qseq file
line = input_file.readline()
parts = line.split()
cluster = parts[13]

if parts[0] in sample_list:
    RAD_sample = [parts[0]]
    RAD_seqs = [parts[8]]
    RAD_qual = [parts[9]]
    weight = int(parts[11])
    RAD_weight = [weight]
    sample_depth[sample_list.index(parts[0])] += weight
else:
    RAD_sample = []
    RAD_seqs = []
    RAD_qual = []
    RAD_weight = []

hom_test1,hom_test2,hom_test3,het_test1,het_test2,het_test3,fail_test1,fail_test2,fail_test3,fail_test4,fail_test5,no_data = 0,0,0,0,0,0,0,0,0,0,0,0

EOF = False
while EOF == False:
    line = input_file.readline()
    parts = line.split()
    if len(parts) < 2:
        EOF = True
        next_cluster = -1
    else:
        next_cluster = parts[13]
        weight = int(parts[11])
        
    if next_cluster == cluster and EOF == False:
        #add to data for current cluster
        if parts[0] in sample_list:
            RAD_sample.append(parts[0])
            RAD_seqs.append(parts[8])        
            RAD_qual.append(parts[9])    
            RAD_weight.append(weight)
            sample_depth[sample_list.index(parts[0])] += weight
            
    
    else:
        #analyze and print results for previous tag/cluster

        #codes: 0, no data; 1, passed; 2, low depth; 3, bad ratio; 4 extra reads; 5, extra reads, possible 3rd allele
        result_codes = num_samples * ['.']

        allele1_temp=[]
        allele2_temp=[]
        
        #count samples with threshold depth
        samples_with_depth = 0
        for i in range(num_samples):
            if sample_depth[i] > 6:
                samples_with_depth += 1
        Loci_with_depth[samples_with_depth] += 1
        #Loci_recovered[num_samples-result_codes.count(0)] += 1

        if sum(sample_depth)<100*num_samples:
            Total_depth[int(round(sum(sample_depth)/100))]+=1
        else:
            Total_depth[num_samples]+=1
        
        dline = bigclusters.readline()
        depth_data = dline.split()
        #print(depth_data)
        if sum(sample_depth) >= BLAST_depth and int(depth_data[3]) < 6000 and cluster not in failed: #NEED to fix to allow any size cluster
            #process each base and identify possible variable sites
            consensus=[]
            var_sites=[]
            total_qual=[]
            for i in range(len(RAD_seqs[0])):
                total_qual.append(0)
                ACGT_count = [0,0,0,0,0,0]
                for j in range(num_samples):
                    RAD_sample_bases[j] = [0,0,0,0,0,0]
                    RAD_sample_qual[j] = [0,0,0,0,0,0]
                    weighted_sample_bases[j] = [0,0,0,0,0,0]
                    sample_weight[j] = 0
                    
                for j in range(len(RAD_seqs)):
                    total_qual[i] += (ord(RAD_qual[j][i])-33)*RAD_weight[j]
                    base = ACGT.index(RAD_seqs[j][i])
                    ACGT_count[base] += int(RAD_weight[j]) * (ord(RAD_qual[j][i])-33)
                    sample = sample_list.index(RAD_sample[j])
                    RAD_sample_bases[sample][base] += int(RAD_weight[j])
                    if ord(RAD_qual[j][i])-33 > RAD_sample_qual[sample][base]:
                        RAD_sample_qual[sample][base] = ord(RAD_qual[j][i])-33

                for j in range(num_samples):
                    weighted_sample_bases[j] = [x*y for x,y in zip(RAD_sample_bases[j],RAD_sample_qual[j])]
                    sample_weight[j] = sum(weighted_sample_bases[j])
                
                major_count,minor_count = heapq.nlargest(2,ACGT_count[0:5])
                major_base = ACGT_count.index(major_count)
                if major_count > minor_count:
                    minor_base = ACGT_count.index(minor_count)
                else:
                    templist = ACGT_count[0:5]
                    templist[templist.index(major_count)]=0
                    minor_base = templist.index(minor_count)
                consensus.append(ACGT[major_base])
                    
                for j in range(num_samples):
                    if sample_depth[j] > 0:
                        sample_major_weight,sample_minor_weight = heapq.nlargest(2,weighted_sample_bases[j])
                        sample_major_base = weighted_sample_bases[j].index(sample_major_weight)
                        
                        # 3 tests for putative polymorphisms
                        if sample_minor_weight == 0:
                            if sample_major_base != major_base and sample_major_weight > 29:
                                var_sites.append(i)
                                break
                        elif sample_major_base == major_base:
                            if (float(sample_minor_weight)/sample_weight[j]) > 0.29 and sample_minor_weight > 29:
                                var_sites.append(i)
                                break
                        elif (float(sample_major_weight)/sample_weight[j]) > 0.29 and sample_major_weight > 29:
                            var_sites.append(i)
                            break

                #reject variable sites in restriction site but allow indels due to "misalignment"
                #if i < len(rsite):
                #    if i in var_sites and ''.join(consensus) == rsite[0:i+1] and minor_base != 4:
                #        var_sites.pop()
                #        print("Popped var site at: ",i,"\tCluster: ",cluster,sample_depth)
                                  
            #trim locus at first low quality var_site
            avg_qual=[]
            for i in range(len(RAD_seqs[0])):
                avg_qual.append(total_qual[i]/sum(sample_depth))
            excluded=[]
            for i in range(len(var_sites)):
                if avg_qual[var_sites[i]]<25:
                    excluded.append(var_sites[i])
                    break
            if len(excluded)>0:
                var_sites = var_sites[:var_sites.index(excluded[0])]
                consensus = consensus[0:min(excluded)]
                          
               
            ##***********************##
            ##****FIND HAPLOTYPES****##
            ##***********************##
            if len(var_sites)> 0:
                newhaps = []
                haplotypes = []
                haplotype_index = []
                haplotype_quals = []
                
                #identify and number haplotypes across all samples, save corresponding quals
                for i in range(len(RAD_seqs)):
                    haplotype = [RAD_seqs[i][j] for j in var_sites]
                    haplotype_qual = [RAD_qual[i][j] for j in var_sites]
                    haplotype_quals.append(haplotype_qual)
                    if haplotype not in haplotypes:
                        haplotypes.append(haplotype)
                        haplotype_index.append(len(haplotypes)-1)    
                    else:
                        haplotype_index.append(haplotypes.index(haplotype))

                #count haplotypes in each sample and save max quals for each
                RAD_sample_haplotypes = [[0]*len(haplotypes) for i in range(num_samples)]
                RAD_sample_quals = [[0]*len(haplotypes) for i in range(num_samples)]
                                                
                for i in range(len(RAD_seqs)):
                    sample = sample_list.index(RAD_sample[i])
                    if 'N' not in haplotypes[haplotype_index[i]]:
                        RAD_sample_haplotypes[sample][haplotype_index[i]] += RAD_weight[i]
                        if RAD_sample_quals[sample][haplotype_index[i]] == 0:
                            RAD_sample_quals[sample][haplotype_index[i]] = haplotype_quals[i]
                        else:
                            for j in range(len(haplotype_quals[i])):
                                if ord(RAD_sample_quals[sample][haplotype_index[i]][j]) < ord(haplotype_quals[i][j]):
                                    RAD_sample_quals[sample][haplotype_index[i]][j] = haplotype_quals[i][j]
                
                           
                #evaluate and print haplotypes for each sample
                for i in range(num_samples):
                    if sum(RAD_sample_haplotypes[i]) > 0:
                        N_haps = sum(x > 0 for x in RAD_sample_haplotypes[i])
                        if N_haps == 1:
                            major_hap_count = max(RAD_sample_haplotypes[i])
                            major_hap = RAD_sample_haplotypes[i].index(major_hap_count)
                            outseq = ''.join(haplotypes[major_hap])
                            if outseq not in newhaps:
                                newhaps.append(outseq)
                            major_quals = RAD_sample_quals[i][major_hap]
                            allele1_temp.append([sample_list[i],outseq,major_hap,major_hap_count,sample_depth[i],major_quals])
                            if sample_depth[i] > 4:
                                allele2_temp.append([sample_list[i],outseq,'.','.','.','.'])
                                result_codes[i] = 1
                                hom_test1 += 1
                            else:
                                allele2_temp.append([sample_list[i],'.','.','.','.','.'])
                                result_codes[i] = 2
                                fail_test1 += 1 #LOW DEPTH
                                
                        else:
                            weighted_haps = [0]*len(haplotypes)
                            weights = []
                            for j in range(len(haplotypes)):
                                if RAD_sample_haplotypes[i][j] > 0:
                                    weights.append(sum(ord(x)-33 for x in RAD_sample_quals[i][j]))
                                else:
                                    weights.append(0)
                            for j in range(len(haplotypes)):        
                                weighted_haps[j] = RAD_sample_haplotypes[i][j] * float(weights[j])/max(weights)
                            
                            top_weights = heapq.nlargest(N_haps,weighted_haps)
                            major_hap = weighted_haps.index(top_weights[0])
                            if top_weights[1] == top_weights[0]:
                                minor_hap = weighted_haps.index(top_weights[1],weighted_haps.index(top_weights[0])+1)
                            else:
                                minor_hap = weighted_haps.index(top_weights[1])
                            
                            major_hap_count = RAD_sample_haplotypes[i][major_hap]
                            major_quals = RAD_sample_quals[i][major_hap]
                            minor_hap_count = RAD_sample_haplotypes[i][minor_hap]
                            minor_quals = RAD_sample_quals[i][minor_hap]
                            count_ratio = (float(minor_hap_count)/(minor_hap_count+major_hap_count))
                            
                            weighted_counts = []
                            for k in range(len(RAD_sample_quals[i])):
                                min_qual_ratio = 1
                                if RAD_sample_quals[i][k] != 0:
                                    for j in range(len(haplotypes[major_hap])):
                                        if haplotypes[major_hap][j] != haplotypes[k][j]:
                                            qual_ratio = (float(ord(RAD_sample_quals[i][k][j])-33))/(ord(major_quals[j])-33)
                                            if qual_ratio < min_qual_ratio:
                                                min_qual_ratio = qual_ratio
                                    weighted_counts.append(min_qual_ratio*RAD_sample_haplotypes[i][k])
                                                        
                            weighted_counts.sort(reverse=True)
                            weighted_ratio1 = float(sum(weighted_counts[1:2]))/(sum(weighted_counts[1:2])+major_hap_count)
                            weighted_ratio2 = float(sum(weighted_counts[2:3]))/(sum(weighted_counts[2:3])+minor_hap_count+major_hap_count)

                            outseq = ''.join(haplotypes[major_hap])
                            if outseq not in newhaps:
                                newhaps.append(outseq)
                            allele1_temp.append([sample_list[i],outseq,major_hap,major_hap_count,sample_depth[i],major_quals])

                            if N_haps == 2:
                                if weighted_ratio1 < 0.07:
                                    if major_hap_count > 4:
                                        allele2_temp.append([sample_list[i],outseq,'.','.','.','.'])
                                        result_codes[i] = 1
                                        hom_test2 += 1
                                    else:
                                        allele2_temp.append([sample_list[i],'.','.','.','.','.'])
                                        result_codes[i] = 2
                                    
                                else:
                                    outseq = ''.join(haplotypes[minor_hap])
                                    if outseq not in newhaps:
                                        newhaps.append(outseq)
                                    allele2_temp.append([sample_list[i],outseq,minor_hap,minor_hap_count,'.',minor_quals])
                                                                        
                                    result_codes[i] = binom_test(minor_hap_count,major_hap_count)                                        
    
                            else: # N_haps = 3 or more 
                                if weighted_ratio1 < 0.07:
                                    allele2_temp.append([sample_list[i],outseq,'.','.','.','.'])
                                    result_codes[i] = 1
                                    hom_test3 += 1

                                else:
                                    outseq = ''.join(haplotypes[minor_hap])
                                    allele2_temp.append([sample_list[i],outseq,minor_hap,minor_hap_count,'.',minor_quals])
                                    if outseq not in newhaps:
                                        newhaps.append(outseq)
                                    if weighted_ratio2 < 0.07 or (weighted_counts[2]<=1 and minor_hap_count>=2):
                                        result_codes[i] = binom_test(minor_hap_count,major_hap_count)
                                    elif weighted_counts[2] > 9 or weighted_counts[2]/sum(weighted_counts) > 0.2:
                                        result_codes[i] = 5 #possible "3rd" allele
                                    else:
                                        result_codes[i] = 4
                    else:
                        #there are no data for sample - print blank lines
                        allele1_temp.append([sample_list[i],'.','.','.','.','.'])
                        allele2_temp.append([sample_list[i],'.','.','.','.','.'])
                        result_codes[i] = 0
                        no_data += 1

                #add result codes to allele_temp lists
                for i in range(num_samples):
                    allele1_temp[i].insert(6,result_codes[i])
                    allele2_temp[i].insert(6,result_codes[i])
                    
                #do pop check?
                if result_codes.count(3)+result_codes.count(4)+result_codes.count(5) < num_samples * 0.80:
                    if result_codes.count(3)>0:
                        check_list=[]
                        for i in allele1_temp:
                            if i[2] != "." and i[3] > 2:
                                if i[1] not in check_list:
                                    check_list.append(i[1])
                        for i in allele2_temp:
                            if i[6] == 1:
                                if i[1] not in check_list:
                                    check_list.append(i[1])
                        
                        for i in range(num_samples):
                            if result_codes[i]==3:
                                if allele2_temp[i][3] > 1:
                                    if allele2_temp[i][1] not in check_list:
                                        #print("BAD RATIO",allele1_temp[i][3],allele2_temp[i][3]) ## 144
                                        result_codes[i] = 3
                                    else:
                                        #print("OK as heterozygote",allele1_temp[i][3],allele2_temp[i][3]) ## 2,401
                                        result_codes[i] = 1
                                else:
                                    if allele1_temp[i][3] > 7:
                                        if allele2_temp[i][1] not in check_list:
                                            #print("OK as homozygote",allele1_temp[i][3],allele2_temp[i][3]) ## 199
                                            allele2_temp[i] = [sample_list[i],allele1_temp[i][1],'.','.','.','.',1]
                                            result_codes[i] = 1
                                        else:
                                            #print("possible heterozygote",allele1_temp[i][3],allele2_temp[i][3]) ## 466
                                            result_codes[i] = 3
                                    else:
                                        #print("possible low depth heterozygote",allele1_temp[i][3],allele2_temp[i][3]) ## 30,169
                                        result_codes[i] = 2
                                        if allele1_temp[i][3] == 1 and allele1_temp[i][1] not in check_list:
                                            if allele2_temp[i][1] in check_list:
                                                allele1_temp[i],allele2_temp[i] = allele2_temp[i],allele1_temp[i] 
                                                                                                                             
                                allele1_temp[i][6] = result_codes[i]
                                allele2_temp[i][6] = result_codes[i]
                                    
                #re-evaluate variable sites and fix output data
                variable=[]
                templist=[list(i) for i in zip(*newhaps)]
                for i in templist:
                    if max(i) != min(i):
                        variable.append(True)
                    else:
                        variable.append(False)

                while False in variable:
                    for i in range(num_samples):
                        if allele1_temp[i][1] != '.':
                            allele1_temp[i][1] = allele1_temp[i][1][:variable.index(False)]+allele1_temp[i][1][variable.index(False)+1:]
                            allele1_temp[i][5] = allele1_temp[i][5][:variable.index(False)]+allele1_temp[i][5][variable.index(False)+1:]
                        if allele2_temp[i][1] != '.':    
                            allele2_temp[i][1] = allele2_temp[i][1][:variable.index(False)]+allele2_temp[i][1][variable.index(False)+1:]
                        if allele2_temp[i][5] != '.':    
                            allele2_temp[i][5] = allele2_temp[i][5][:variable.index(False)]+allele2_temp[i][5][variable.index(False)+1:]
                    var_sites.pop(variable.index(False))
                    variable.pop(variable.index(False))
                                
                #renumber haplotypes
                if len(var_sites) > 0:
                    newhaps2=[]
                    for i in range(num_samples):
                        if allele1_temp[i][1] != '.':
                            if allele1_temp[i][1] not in newhaps2:
                                newhaps2.append(allele1_temp[i][1])
                            allele1_temp[i][2]=newhaps2.index(allele1_temp[i][1])
                            if allele2_temp[i][1] != '.':
                                if allele2_temp[i][1] not in newhaps2:
                                    newhaps2.append(allele2_temp[i][1])
                                allele2_temp[i][2]=newhaps2.index(allele2_temp[i][1])
                    
                    #reconstruct full seqs HERE!
                    fullhaps=[]
                    for i in range(len(newhaps2)):
                        hapseq = []
                        for m in range(len(consensus)):
                            if m in var_sites:
                                hapseq.append(newhaps2[i][var_sites.index(m)])
                            else:
                                hapseq.append(consensus[m])
                        hapstring = ''.join(hapseq)
                        fullhaps.append(hapstring)
                    
                    #print(fullhaps)
                    for i in range(num_samples):
                        if allele1_temp[i][5] != '.':
                            major_quals = allele1_temp[i][5]
                            major_quals_print =('%02d' % (ord(major_quals[0])-33))
                            for z in range(len(major_quals)-1):
                                major_quals_print = major_quals_print + ('.%02d' % (ord(major_quals[z+1])-33))
                            allele1_temp[i][5]=major_quals_print
                            #print(major_quals_print)
                        if allele2_temp[i][5] != '.':
                            minor_quals = allele2_temp[i][5]
                            minor_quals_print =('%02d' % (ord(minor_quals[0])-33))
                            for z in range(len(minor_quals)-1):
                                minor_quals_print = minor_quals_print + ('.%02d' % (ord(minor_quals[z+1])-33))
                            allele2_temp[i][5]=minor_quals_print

                        hap_index = allele1_temp[i][2]
                        if hap_index != '.':
                            allele1_temp[i].insert(1,fullhaps[hap_index])
                        else:
                            allele1_temp[i].insert(1,'.')
                        hap_index = allele2_temp[i][2]
                        if hap_index != '.':
                            allele2_temp[i].insert(1,fullhaps[hap_index])
                        else:
                            allele2_temp[i].insert(1,'.')

                    #strip alignment gaps
                    conseq=''.join(consensus)
                    if conseq.find('-')>-1:
                        conseq,allele1_temp,allele2_temp,var_sites = strip_alignment_gaps(conseq,allele1_temp,allele2_temp,var_sites)
                        consensus = list(conseq)
                                        
            #####################################            
            ### NEW CODE TO REFINE ALIGNMENTS ###
            #####################################            
            edit_code = 0
            if len(var_sites)> 0:
                ###check for indels###
                gaps = 0
                templist=[list(i) for i in zip(*newhaps2)]
                for i in range(len(templist)):
                    if '-' in templist[i]:
                        ###find unique gaps###
                        seqs_list1 = [x[1] for x in allele1_temp]+[y[1] for y in allele2_temp]
                        seqs_list2 = list(set(seqs_list1))
                        unique_gaps = []
                        for i in seqs_list2:
                            if i.find('-')>-1:
                                gap_pos = [m.start() for m in re.finditer('-',i)]
                                starts = [1]+[0] * (len(gap_pos)-1)
                                ends = [0] * (len(gap_pos)-1)+[1]
                                for j in range(len(gap_pos)-1):
                                    if gap_pos[j+1] > gap_pos[j]+1:
                                        starts[j+1] = 1
                                    if gap_pos[j] < gap_pos[j+1]-1:
                                        ends[j] = 1
                                start_pos = [k*l for k,l in zip(gap_pos,starts) if k*l > 0]
                                end_pos = [k*l for k,l in zip(gap_pos,ends) if k*l > 0]
                                if len(start_pos) < len(end_pos):
                                    start_pos.insert(0,0)
                                for j,k in zip(start_pos,end_pos):
                                    if [j,k] not in unique_gaps:
                                        unique_gaps.append([j,k])
                        unique_gaps.sort(reverse=True)
                        unique_gaps = sorted(unique_gaps, key=lambda x: x[1])
                        gaps = len(unique_gaps)
                        break
                        
                ###flag no-indel loci with >1 var sites near end###
                if gaps == 0:
                    endsites = [x for x in var_sites if x > len(conseq)-6]
                    if len(endsites) > 1:
                        edit_code = 1 #possible misalignment at end of sequence
                            
                else:
                    ###check for trailing gap and trim right end###
                    if unique_gaps[-1][1]==len(consensus)-1:
                        endsites=0
                        for i in var_sites:
                            if i > unique_gaps[-1][0]-6 and i < unique_gaps[-1][0]:
                                endsites+=1
                        if endsites > 1:
                            edit_code = 2 #ragged end with too many adjacent var sites
                        else:
                            ###clean ragged ends, renumber haps, adjust quals###
                            edit_code = 7 #trimmed ragged right end
                            var_sites = var_sites[:var_sites.index(unique_gaps[-1][0])]
                            consensus = consensus[:unique_gaps[-1][0]]
                            newhaps2=[]
                            for i in range(num_samples):
                                if allele1_temp[i][1] != '.':
                                    allele1_temp[i][1]=allele1_temp[i][1][:unique_gaps[-1][0]]
                                    if allele1_temp[i][1] not in newhaps2:
                                        newhaps2.append(allele1_temp[i][1])
                                        allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                    else:
                                        allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                    allele1_temp[i][2] = allele1_temp[i][2][:len(var_sites)]
                                    allele1_temp[i][6] = allele1_temp[i][6][:len(var_sites)*3-1]
                                if allele2_temp[i][1]!='.':
                                    allele2_temp[i][1]=allele2_temp[i][1][:unique_gaps[-1][0]]
                                    if allele2_temp[i][1] not in newhaps2:
                                        newhaps2.append(allele2_temp[i][1])
                                        allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                    else:
                                        allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                    allele2_temp[i][2] = allele2_temp[i][2][:len(var_sites)]
                                    allele2_temp[i][6] = allele2_temp[i][6][:len(var_sites)*3-1]
                                
                            for i in range(len(unique_gaps)-1,-1,-1):
                                if unique_gaps[i][0] >= len(newhaps2[0]):
                                    unique_gaps.pop(i)

                    gaps = len(unique_gaps)

                    ###check for leading gap and adjust accordingly###
                    if gaps > 0:
                        unique_gaps = sorted(unique_gaps, key=lambda x: x[0])
                        if unique_gaps[0][0] == 0:
                            if len(unique_gaps) > 1 and unique_gaps[1][0] <= unique_gaps[0][1] and unique_gaps[1][1] > unique_gaps[0][1]:
                                if edit_code == 0:
                                    edit_code = 3 #overlapping gaps at front end
                                else:
                                    edit_code += 30
                            else:
                                if edit_code == 0:
                                    edit_code = 8 #leading gaps trimmed
                                else:
                                    edit_code += 80
                                trim = unique_gaps[0][1]+1
                                consensus = consensus[trim:]                      
                                var_sites = var_sites[var_sites.index(unique_gaps[0][1])+1:]
                                newhaps2=[]
                                for i in range(num_samples):
                                    if allele1_temp[i][1]!='.':
                                        allele1_temp[i][1]=allele1_temp[i][1][unique_gaps[0][1]+1:]
                                        if allele1_temp[i][1] not in newhaps2:
                                            newhaps2.append(allele1_temp[i][1])
                                            allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                        else:
                                            allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                        allele1_temp[i][2] = allele1_temp[i][2][len(allele1_temp[i][2])-len(var_sites):]
                                        allele1_temp[i][6] = allele1_temp[i][6][(len(allele1_temp[i][6])-len(var_sites)*3)+1:]

                                    if allele2_temp[i][1]!='.':
                                        allele2_temp[i][1]=allele2_temp[i][1][unique_gaps[0][1]+1:]
                                        if allele2_temp[i][1] not in newhaps2:
                                            newhaps2.append(allele2_temp[i][1])
                                            allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                        else:
                                            allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                        allele2_temp[i][2] = allele2_temp[i][2][len(allele2_temp[i][2])-len(var_sites):]
                                        allele2_temp[i][6] = allele2_temp[i][6][(len(allele2_temp[i][6])-len(var_sites)*3)+1:]
                                
                                for i in range(len(unique_gaps)-1,-1,-1):
                                    if len(var_sites) > 0:
                                        if unique_gaps[i][0] < var_sites[0]:
                                            unique_gaps.pop(i)
                                    else:
                                        unique_gaps = []
                                for i in range(len(var_sites)):
                                    var_sites[i] -= trim
                                for i in range(len(unique_gaps)):
                                    unique_gaps[i][0] -= trim
                                    unique_gaps[i][1] -= trim                                                                      
                                
                        else:
                            ###test for slide right###
                            if unique_gaps[0][0] < len(rsite)+1 and unique_gaps[0][1]-unique_gaps[0][0] > len(rsite)-2:
                                if len(unique_gaps) > 1 and unique_gaps[1][0] <= unique_gaps[0][1] and unique_gaps[1][1] >= unique_gaps[0][1]:
                                    if edit_code == 0:
                                        edit_code = 4 #overlapping gaps for slide right
                                    else:
                                        edit_code += 40               
                                else:                               
                                    seqlist=[]
                                    for i in range(num_samples):
                                        if allele1_temp[i][1]!='.':
                                            seqlist.append(allele1_temp[i][1])
                                        if allele2_temp[i][1]!='.':
                                            seqlist.append(allele2_temp[i][1])
                                    bconseq = []
                                    seqlist2 = [list(i) for i in zip(*seqlist)]
                                    for i in seqlist2:
                                        bases = list(filter(('-').__ne__, list(set(i))))
                                        if len(bases)>0:
                                            bconseq.append(max(set(bases), key=bases.count))
                                        else:
                                            bconseq.append('-')
                                    gconseq = ''.join(bconseq)
                                    string1 = gconseq[:unique_gaps[0][0]]
                                    string2 = gconseq[unique_gaps[0][1]-len(string1)+1:unique_gaps[0][1]+1]
                                    diffs = sum(i != j for i,j in zip(string1,string2))
                                    #slide right
                                    if diffs == 1:
                                        if edit_code == 0:
                                            edit_code = 9 #SbfI site moved to right and trimmed
                                        else:
                                            edit_code += 90
                                        trim = unique_gaps[0][1]+1-unique_gaps[0][0]
                                        consensus = consensus[:unique_gaps[0][0]]+consensus[unique_gaps[0][1]+1:]
                                        testgap = '-'*(trim)
                                        var_sites[:] = [x - trim for x in var_sites]
                                        var_sites = [x for x in var_sites if x >= 0]
                                        newhaps2=[]
                                        for i in range(num_samples):
                                            if allele1_temp[i][1]!='.':
                                                if allele1_temp[i][1][unique_gaps[0][0]:unique_gaps[0][1]+1] == testgap:
                                                    allele1_temp[i][1]=allele1_temp[i][1][:unique_gaps[0][0]]+allele1_temp[i][1][unique_gaps[0][1]+1:]
                                                else:
                                                    allele1_temp[i][1]=allele1_temp[i][1][unique_gaps[0][1]+1-unique_gaps[0][0]:]
                                                if allele1_temp[i][1] not in newhaps2:
                                                    newhaps2.append(allele1_temp[i][1])
                                                    allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                                else:
                                                    allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                                allele1_temp[i][2] = allele1_temp[i][2][-len(var_sites):]                                                

                                            if allele2_temp[i][1]!='.':
                                                if allele2_temp[i][1][unique_gaps[0][0]:unique_gaps[0][1]+1] == testgap:
                                                    allele2_temp[i][1]=allele2_temp[i][1][:unique_gaps[0][0]]+allele2_temp[i][1][unique_gaps[0][1]+1:]
                                                else:
                                                    allele2_temp[i][1]=allele2_temp[i][1][unique_gaps[0][1]+1-unique_gaps[0][0]:]
                                                if allele2_temp[i][1] not in newhaps2:
                                                    newhaps2.append(allele2_temp[i][1])
                                                    allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                                else:
                                                    allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                                allele2_temp[i][2] = allele2_temp[i][2][-len(var_sites):]                                                

                                        for i in range(len(unique_gaps)-1,-1,-1):    
                                            if unique_gaps[i][0] < trim:
                                                unique_gaps.pop(i)
                                                
                                        #re-evaluate var_sites:
                                        var_sites = []
                                        templist=[list(i) for i in zip(*newhaps2)]
                                        for i in range(len(templist)):
                                            if len(list(set(templist[i]))) > 1:
                                                var_sites.append(i)
                                        for i in range(num_samples):
                                            if allele1_temp[i][1]!='.':
                                                allele1_temp[i][2] = ''.join(allele1_temp[i][1][j] for j in var_sites)
                                                allele1_temp[i][6] = allele1_temp[i][6][-(len(var_sites)*3-1):]
                                            if allele2_temp[i][1]!='.':
                                                allele2_temp[i][2] = ''.join(allele2_temp[i][1][j] for j in var_sites)
                                                allele2_temp[i][6] = allele2_temp[i][6][-(len(var_sites)*3-1):] 
                                        for i in range(len(unique_gaps)):
                                            unique_gaps[i][0] -= trim
                                            unique_gaps[i][1] -= trim                                                                      
                                                                        
                                    else:
                                        if edit_code == 0:
                                            edit_code = 5 #possible slide right, but too many mismatches
                                        else:
                                            edit_code += 50
                                                        
                            elif unique_gaps[0][0] < len(rsite):
                                if edit_code == 0:
                                    edit_code = 10 #small gap in rsite
                                else:
                                    edit_code += 100
                                
                    gaps = len(unique_gaps)
                    
                    ###check for split EcoRI and slide left###
                    if gaps > 0:
                        strayC=False
                        for i in range(num_samples):
                            if 'GAATT-' in allele1_temp[i][1]:
                                endstring = allele1_temp[i][1][allele1_temp[i][1].find("GAATT-")+6:]
                                endstring = ''.join(x for x in endstring if x != '-')
                                if endstring == 'C':
                                    SCstart = allele1_temp[i][1].find("GAATT-")+5
                                    SCend = allele1_temp[i][1][SCstart+1:].find("C")+SCstart
                                    last_gap = [SCstart,SCend]
                                    strayC = True
                                    break
                            if 'GAATT-' in allele2_temp[i][1]:
                                endstring = allele2_temp[i][1][allele2_temp[i][1].find("GAATT-")+6:]
                                endstring = ''.join(x for x in endstring if x != '-')
                                if endstring == 'C':
                                    SCstart = allele2_temp[i][1].find("GAATT-")+5
                                    SCend = allele2_temp[i][1][SCstart+1:].find("C")+SCstart
                                    last_gap = [SCstart,SCend]
                                    strayC = True
                                    break
                        if strayC == True:
                            OK_to_trim = False
                            if last_gap == unique_gaps[-1]:
                                OK_to_trim = True
                                if len(unique_gaps) > 1:
                                    for i in unique_gaps[:-1]:
                                        if i[0]<last_gap[0] and i[1]>last_gap[0]:
                                            OK_to_trim = False
                            if OK_to_trim == True:
                                if edit_code == 0:
                                    edit_code = 11 #stray C moved to left and trimmed
                                else:
                                    edit_code += 110
                                consensus = consensus[:last_gap[0]]+consensus[last_gap[1]+1:]
                                newhaps2=[]
                                testgap = '-'*(last_gap[1]+1-last_gap[0])
                                for i in range(num_samples):
                                    if allele1_temp[i][1]!='.':
                                        if allele1_temp[i][1][last_gap[0]:last_gap[1]+1] == testgap:
                                            allele1_temp[i][1]=allele1_temp[i][1][:last_gap[0]]+allele1_temp[i][1][last_gap[1]+1]
                                        else:
                                            allele1_temp[i][1]=allele1_temp[i][1][:last_gap[0]+1]
                                        if allele1_temp[i][1] not in newhaps2:
                                            newhaps2.append(allele1_temp[i][1])
                                            allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                        else:
                                            allele1_temp[i][3] = newhaps2.index(allele1_temp[i][1])
                                        

                                    if allele2_temp[i][1]!='.':
                                        if allele2_temp[i][1][last_gap[0]:last_gap[1]+1] == testgap:
                                            allele2_temp[i][1]=allele2_temp[i][1][:last_gap[0]]+allele2_temp[i][1][last_gap[1]+1]
                                        else:
                                            allele2_temp[i][1]=allele2_temp[i][1][:last_gap[0]+1]
                                        if allele2_temp[i][1] not in newhaps2:
                                            newhaps2.append(allele2_temp[i][1])
                                            allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                        else:
                                            allele2_temp[i][3] = newhaps2.index(allele2_temp[i][1])
                                        
                                unique_gaps.pop(-1)
                                gaps = len(unique_gaps)
                                        
                                #re-evaluate var_sites
                                var_sites = []
                                templist=[list(i) for i in zip(*newhaps2)]
                                for i in range(len(templist)):
                                    if len(list(set(templist[i]))) > 1:
                                        var_sites.append(i)
                                for i in range(num_samples):
                                    if allele1_temp[i][1]!='.':
                                        allele1_temp[i][2] = ''.join(allele1_temp[i][1][j] for j in var_sites)
                                        allele1_temp[i][6] = allele1_temp[i][6][:len(var_sites)*3-1]
                                    if allele2_temp[i][1]!='.':
                                        allele2_temp[i][2] = ''.join(allele2_temp[i][1][j] for j in var_sites)
                                        allele2_temp[i][6] = allele2_temp[i][6][:len(var_sites)*3-1] 
                            else:
                                if edit_code == 0:
                                    edit_code = 6 #strayC but not OK to trim
                                else:
                                    edit_code += 60
                                
                                       
            if len(var_sites) == 0:    
                #SNPs,gaps,infinite_sites = 0,0,True
                for i in range(num_samples):
                    if sample_depth[i] > 0:
                        result_codes[i] = 1
                    else:
                        result_codes[i] = 0

            results=[]
            for i in range(6):
                results.append(result_codes.count(i))
            
            if results[0] > (num_samples * 0.2):
                if len(var_sites) == 0:
                    outfile_name = base_filename + 'CL.out'
                elif sum(results[3:]) == 0:
                    outfile_name = base_filename + 'VL.out'
                elif sum(results[3:]) < 3:
                    outfile_name = base_filename + 'VLX2.out'
                else:
                    outfile_name = base_filename + 'VLXX.out'      
            else:
                if len(var_sites) == 0:
                    outfile_name = base_filename + 'C.out'
                elif sum(results[3:]) == 0:
                    outfile_name = base_filename + 'V.out'
                elif sum(results[3:]) < 3:
                    outfile_name = base_filename + 'VX2.out'
                else:
                    outfile_name = base_filename + 'VXX.out'

            temp_name = outfile_name[len(base_filename):]
            short_name = temp_name[:temp_name.find(".out")]
            ##1 new line of code version 1.09
            outfile_name = base_filename + 'out'
            
            for i in range(num_samples):
                vars()['results'+short_name][i][result_codes[i]]+=1
            

            vars()[outfile_name].write('Clstr: '+str(cluster)+'\tTot_Depth:\t'+str(sum(sample_depth)))

            #get correct BLAST results
            if sum(sample_depth) >= BLAST_depth:
                while cluster!=BLASTparts[0]:
                    BLASTline = BLASTinput.readline()
                    BLASTparts = BLASTline.strip('\n').split('\t')
            else:
                BLASTparts[2]=('NEED')
                BLASTparts[3]=('TO')
                BLASTparts[10]=('SEARCH')

            print("CLUSTER",cluster,BLASTparts[2],BLASTparts[3],BLASTparts[4],BLASTparts[5])
            vars()[outfile_name].write('\tBLAST_result:\t'+BLASTparts[2]+'\t'+BLASTparts[3]+'\t'+BLASTparts[4]+'\t'+BLASTparts[5])

            
            heterozy = 0
            if len(var_sites)>0:
                #check for homozygotes generated by alignment edits
                if edit_code > 6:
                    for i in range(num_samples):
                        if allele2_temp[i][4] != '.' and allele1_temp[i][3] == allele2_temp[i][3]:
                            allele1_temp[i][4] = allele1_temp[i][4] + allele2_temp[i][4]
                            #print(allele1_temp[i][4])
                            allele2_temp[i][4] = '.'
                            allele2_temp[i][6] = '.'
                            if allele1_temp[i][4] < 5:
                                allele2_temp[i][1] = '.'
                                allele2_temp[i][2] = '.'
                haps = len(newhaps2)
                infinite_sites = True
                conseq=''.join(consensus)
                if gaps > 0:
                    #reconfigure var sites and gaps
                    SNPs = 0
                    SNP_list=[]
                    newhaps2=[]
                    for i in range(num_samples):
                        if allele1_temp[i][1]!='.':
                            if allele1_temp[i][1] not in newhaps2:
                                newhaps2.append(allele1_temp[i][1])                        
                        if allele2_temp[i][1]!='.':
                            if allele2_temp[i][1] not in newhaps2:
                                newhaps2.append(allele2_temp[i][1])
                    templist=[list(i) for i in zip(*newhaps2)]
                    for i in range(len(templist)):
                        nucls = list(set([x for x in templist[i] if x !='-']))
                        if len(nucls) > 1:
                            SNP_list.append(i)
                            SNPs+=1
                    var_sites = SNP_list+unique_gaps
                    all_sites = SNP_list[:]
                    for i in unique_gaps:
                        for j in range(i[0],i[1]+1):
                            all_sites.append(j)
                    all_sites = list(set(all_sites))
                    all_sites.sort()
                    retain_quals = []
                    for i in all_sites:
                        if i in SNP_list:
                            retain_quals.append(True)
                        else:              
                            retain_quals.append(False)
                    
                    for i in range(num_samples):
                        if allele1_temp[i][1]!='.':
                            allele1_temp[i][2] = ''.join(allele1_temp[i][1][j] for j in SNP_list)
                            for j in unique_gaps:
                                if len(set(allele1_temp[i][1][j[0]:j[1]+1])) == 1 and allele1_temp[i][1][j[0]] == '-':
                                    allele1_temp[i][2] = allele1_temp[i][2]+'0'
                                else:
                                    allele1_temp[i][2] = allele1_temp[i][2]+'1'                                    
                            allele1_temp[i][6] = allele1_temp[i][6].split('.')
                            allele1_temp[i][6] = [x for x,y in zip(allele1_temp[i][6],retain_quals) if y == True]
                            allele1_temp[i][6] = '.'.join(allele1_temp[i][6])
                        if allele2_temp[i][1]!='.':
                            allele2_temp[i][2] = ''.join(allele2_temp[i][1][j] for j in SNP_list)
                            for j in unique_gaps:
                                if len(set(allele2_temp[i][1][j[0]:j[1]+1])) == 1 and allele2_temp[i][1][j[0]] == '-':
                                    allele2_temp[i][2] = allele2_temp[i][2]+'0'
                                else:
                                    allele2_temp[i][2] = allele2_temp[i][2]+'1'                                    
                            allele2_temp[i][6] = allele2_temp[i][6].split('.')
                            allele2_temp[i][6] = [x for x,y in zip(allele2_temp[i][6],retain_quals) if y == True]
                            allele2_temp[i][6] = '.'.join(allele2_temp[i][6])
                else:
                    SNPs = len(var_sites)

                if infinite_sites == True:
                    templist=[list(i[2]) for i in allele1_temp if i[2] != '.']+[list(i[2]) for i in allele2_temp if i[2] != '.']
                    siteslist = [list(i) for i in zip(*templist)]
                    for i in siteslist:
                        if len(set(i)) > 2:
                            infinite_sites = False
                            break
                if infinite_sites == True and len(var_sites) > 1:
                    for i in range(len(siteslist)):
                        for j in range(i+1,len(siteslist)):
                            twobp_haps = []
                            for k in range(len(siteslist[i])):
                                twobp = (siteslist[i][k]+siteslist[j][k])
                                if twobp not in twobp_haps:
                                    twobp_haps.append(twobp)
                            if len(twobp_haps) > 3:
                                infinite_sites = False
                                break
                #HW test
                genotypes = [[i[3],j[3]] for i,j in zip(allele1_temp,allele2_temp)]
                homozy = [1 for i in genotypes if i[0]==i[1] and i[1]!='.']
                heterozy = [1 for i in genotypes if i[0]!=i[1] and i[1]!='.']
                    
                #calculate depth flag
                counts = [[i[4],j[4]] for i,j in zip(allele1_temp,allele2_temp) if j[3] != '.']
                hom_depths = [x[0] for x in counts if x[1] == '.']
                het_depths = [sum(x) for x in counts if x[1] != '.']
                if sum(heterozy) > 0 and sum(homozy) > 0:
                    depth_ratio = median(het_depths)/median(hom_depths)
                    depth_ratio = '%.2f' % depth_ratio
                else:
                    depth_ratio = '.'
                ###############
                                
                alleles = [i[3]for i in allele1_temp if i[3] != '.']+[i[3]for i in allele2_temp if i[3] != '.']
                allele_freqs = dict((i,alleles.count(i)) for i in alleles)
                allele_freqs = [allele_freqs[i] for i in allele_freqs]
                
                exp_homozy = sum(i*i for i in allele_freqs)/(sum(allele_freqs)*sum(allele_freqs))
                observed = [sum(homozy),sum(heterozy)]
                expected = [exp_homozy*sum(observed),(1-exp_homozy)*sum(observed)]
                if sum(observed) > 0:
                    chisquare = sum([(((i-j)*(i-j))/j) for i,j in zip(observed,expected)])
                    if observed[0] > expected[0]:
                        direction = '\t-'
                    elif observed[0] < expected[0]:
                        direction = '\t+'
                    else:
                        direction ='\t.'
                    HW = "%.2f" % chisquare + direction
                else:
                    HW = '.\t.'
                heterozy=sum(heterozy)
                                
                             
                #write excluded sites to output file
                if len(excluded) > 0:
                    vars()[outfile_name].write('\tExcl: '+str(excluded[0])+' ('+str(avg_qual[excluded[0]])+')')
                    for i in range(len(excluded)-1):
                        vars()[outfile_name].write(', '+str(excluded[i+1])+' ('+str(avg_qual[excluded[i+1]])+')')
                    vars()[outfile_name].write('\tTrunc: '+str(min(excluded))+'\n')
                else:
                    vars()[outfile_name].write('\n')

                #write variable sites to output file
                vars()[outfile_name].write('SNPs: ')
                if SNPs > 0:
                    vars()[outfile_name].write(str(var_sites[0]+1))
                    for i in range(SNPs-1):
                        vars()[outfile_name].write(', '+str(var_sites[i+1]+1))
                else:
                    vars()[outfile_name].write('none')
                vars()[outfile_name].write('\tIndels: ')
                if gaps > 0:
                    if unique_gaps[0][1] > unique_gaps[0][0]:
                        vars()[outfile_name].write(str(unique_gaps[0][0])+'-'+str(unique_gaps[0][1]))
                    else:
                        vars()[outfile_name].write(str(unique_gaps[0][0]))
                    for i in range(gaps-1):
                        if unique_gaps[i+1][1] > unique_gaps[i+1][0]:
                            vars()[outfile_name].write(', '+str(unique_gaps[i+1][0])+'-'+str(unique_gaps[i+1][1]))
                        else:
                            vars()[outfile_name].write(', '+str(unique_gaps[i+1][0]))
                else:
                    vars()[outfile_name].write('none')
                vars()[outfile_name].write('\n')

                #write data for each sample
                for i in range(num_samples):
                    for j in range(len(allele1_temp[i])-1):
                        vars()[outfile_name].write(str(allele1_temp[i][j])+'\t')
                    vars()[outfile_name].write(str(allele1_temp[i][-1])+'\n')
                    for j in range(len(allele2_temp[i])-1):
                        vars()[outfile_name].write(str(allele2_temp[i][j])+'\t')
                    vars()[outfile_name].write(str(allele2_temp[i][-1])+'\n')
                    
            else:
                haps,SNPs,gaps,infinite_sites,HW,depth_ratio = 1,0,0,True,'.\t.','.'
                vars()[outfile_name].write('\nSNPs: none\tIndels: none\n')
                conseq=''.join(x for x in consensus if x != '-')
                consensus = conseq
                for m in range(num_samples):
                    if sample_depth[m] > 4:
                        vars()[outfile_name].write(sample_list[m]+'\t'+conseq+'\t.\t0\t'+ str(sample_depth[m])+'\t'+str(sample_depth[m])+'\t.\t1\n')
                        vars()[outfile_name].write(sample_list[m]+'\t'+conseq+'\t.\t0\t.\t.\t.\t1\n')
                    elif sample_depth[m] > 0:
                        vars()[outfile_name].write(sample_list[m]+'\t'+conseq+'\t.\t0\t'+ str(sample_depth[m])+'\t'+str(sample_depth[m])+'\t.\t2\n')
                        vars()[outfile_name].write(sample_list[m]+'\t.\t.\t.\t.\t.\t.\t2\n')
                    else:
                        vars()[outfile_name].write(sample_list[m]+'\t.\t.\t.\t.\t.\t.\t0\n')
                        vars()[outfile_name].write(sample_list[m]+'\t.\t.\t.\t.\t.\t.\t0\n')            

            
            result = '\t'.join(str(i) for i in results)
            sumfile.write(str(cluster)+'\t'+BLASTparts[2]+'\t'+BLASTparts[3]+'\t'+BLASTparts[4]+'\t'+BLASTparts[5]
                +'\t'+short_name+'\t'+str(len(consensus))+'\t'+str(len(var_sites))+'\t'+str(SNPs)+'\t'+str(gaps)+'\t'
                +str(haps)+'\t'+str(edit_code)+'\t'+str(infinite_sites)+'\t'+HW+'\t'+depth_ratio+'\t'+str(heterozy)+'\t'+result+'\t'+conseq+'\t')
            depths = '\t'.join(str(i) for i in sample_depth)
            sumfile.write(depths+'\t')
            codes = '\t'.join(str(i) for i in result_codes)
            sumfile.write(codes+'\n')
        
        if EOF == False:
            sample_depth=[]
            for i in range(num_samples):
                sample_depth.append(0)
            cluster=next_cluster
            if parts[0] in sample_list:
                RAD_seqs = [parts[8]]
                RAD_qual = [parts[9]]
                RAD_sample = [parts[0]]
                RAD_weight = [weight]
                sample_depth[sample_list.index(parts[0])] += weight
            else:
                RAD_seqs = []
                RAD_qual = []
                RAD_sample = []
                RAD_weight = []


for i in range(num_samples):
    print(sample_list[i],resultsC[i],resultsCL[i],resultsV[i],resultsVX2[i],resultsVL[i],resultsVLX2[i],resultsVXX[i],resultsVLXX[i])

print('Loci_recovered: ',Loci_recovered)
print('Loci_with_depth: ',Loci_with_depth)
print('Constant_loci: ',Constant_loci)
print('Passed_variable_loci: ',Passed_variable_loci)
print('Imperfect_variable_loci: ',Imperfect_variable_loci)
print('Total_depth: ',Total_depth)
vars()[base_filename + 'out'].close()
print(hom_test1,hom_test2,hom_test3,het_test1,het_test2,fail_test1,fail_test2,fail_test3,fail_test4,fail_test5,no_data)
sumfile.close()
