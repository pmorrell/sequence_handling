#!/usr/bin/env python

#   This code is based on code originally written by Tom Kono
#   Used to be named Filter_VCF.py

#   A script to apply various arbitrary filters to a VCF
#       - minimum # reads
#       - maximum # reads
#       - read balance in heterozygotes
#       - GQ
#       - DP

#   Some versions of GATK may produce incompatible output for this script
#   This script writes the filtered VCF lines to standard output

import sys

#   If variants have below a PHRED-scaled quality of 40,
#   The quality score in the 5th column, we did the quantile for those values and found that 40 is correct score.
#quality_cutoff = 40
quality_cutoff = sys.argv[2]

#Filter the 90% of the number of heterozygotes for genotype call
het_cutoff = 161
#   We exclude 20% missing data
missing_cutoff = 36

#   If we have low genotype confidence, then we also want to exclude the SNP
#   The genotype qualities (GQ) are also stored as a PHRED-scaled probability. I calculate the quantile of the GQ and 10% of the reads, GQ>=9, So we pick 9 as thrshold 
gt_cutoff = 9
#n_gt_cutoff = 1
#   Our coverage cutoff is 5 reads per sample
per_sample_coverage_cutoff = 5
#all of the samples should be have high coverage
#n_low_coverage_cutoff = 0
#   The number of samples
#nsam = 1

#   Read the file in line-by-line
with open(sys.argv[1]) as f:
    for line in f:
        #   Skip the header lines - write them out without modification
        if line.startswith('#'):
            sys.stdout.write(line)
        else:
            tmp = line.strip().split('\t')
            #   If indel or more than 2 alleles, skip the site - don't write it
            if len(tmp[3]) != 1 or len(tmp[4]) != 1:
                continue
            format = tmp[8]
            sample_information = tmp[9:]
            #   Start counting up # of hets, # low quality genotypes, and # low coverage sites
            nhet = 0
            low_qual_gt = 0
            low_coverage = 0
            missing_data = 0
            for s in sample_information:
                #   For the GATK Variant Recalibrator, the per-sample information has GT:AD:DP:GQ:PGT:PID:PL or GT:AD:DP:GQ:PL
                info = s.split(':')
                gt = info[0]
                #   We have to check for missing data first, because if it is missing, then the other fields are not filled in
                if '.' in gt:
                    missing_data += 1
                else:
                    dp = info[2]
                    gq = info[3]
                    #   Split on / (we have unphased genotypes) and count how many hets we have
                    if len(set(gt.split('/'))) > 1:
                        nhet += 1
                    if dp == '.' or int(dp) < per_sample_coverage_cutoff or gq == '.' or int(gq) < gt_cutoff:
                        #low_coverage += 1
                        missing_data += 1
                    #if gq == '.' or int(gq) < gt_cutoff:
                    #    low_qual_gt += 1
                #print (s,missing_data)        
            #   The quality score is the sixth element
            #   If the quality score is missing or the site quality score is too low or any of the counts are above the cutoff, don't print the site
            if tmp[5] == '.' or float(tmp[5]) < quality_cutoff or nhet > het_cutoff  or missing_data > missing_cutoff:
            #if tmp[5] == '.' or float(tmp[5]) < quality_cutoff or nhet > het_cutoff or low_qual_gt > n_gt_cutoff or low_coverage > n_low_coverage_cutoff or missing_data > missing_cutoff:
               continue
            else:
                sys.stdout.write(line)