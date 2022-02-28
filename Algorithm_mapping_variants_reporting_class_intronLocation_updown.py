# algorithm for mapping identified variants onto the genomic location and reporting the hit class 

# python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf

import sys
from sys import argv
import csv

exon_file = None
snp_file = None
exon_buffer = None
intron_option = None

# increase field limit
limit = sys.maxsize
while True:
    # decrease limit by factor 10 as long as OverflowError occurs
    try:
        csv.field_size_limit(limit)
        break
    except OverflowError:
        limit = int(limit/10)

        
if(len(argv) == 3):
    exon_file = argv[1]
    snp_file = argv[2]
    print('The program assumes that you do NOT want to report how far the variant falls in the exon boundary.')
elif(len(argv) == 4):
    exon_file = argv[1]
    snp_file = argv[2]
    exon_buffer = int(argv[3])
    print('The program assumes that you DO want to report how far the variant falls in the exon boundary.')
    print('Only variants flanking their nearby exon within <=' + str(exon_buffer) + ' bp are reported')
    intron_option = True
else:
    print('The input commands do not meet the requirement. Please see the README file and try again.')
    quit()
    
output_file = snp_file + '.append'

class FeatureType:
    def __init__(self, name, ext):
        self.name = name
        self.ext = ext
        self.chrom = {} # start position -> stop position
        self.gene = {} # chrom_start -> stop position
        self.starts = {} # chrom -> start positions
        
        for line in list(csv.reader(open(exon_file + '.' + ext + '_link_shrink'), delimiter='\t')):
            self.chrom[line[0]] = int(line[1])
        for line in list(csv.reader(open(exon_file + '.' + ext + '_gene'), delimiter='\t')):
            self.gene[line[0]] = line[1]
        
        # store start positions for binary search
        for line in list(csv.reader(open(exon_file + '.' + ext), delimiter='\t')):
            chrom = line[0]
            arr = []
            for start in line[1:]:
                if start != '':
                    arr.append(int(start))
            if chrom not in self.starts:
                self.starts[chrom] = arr
            else:
                self.starts.extend(arr)
        
                        
def binary_search(starts, snp_start):
    lo = 0
    hi = len(starts) - 1
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if snp_start >= starts[mid]:
            lo = mid
        else:
            hi = mid - 1
    return lo

feature_types = [FeatureType('CDSHIT', 'cds'), FeatureType('INTRON', 'intron'), FeatureType('UTR5', 'utr5'), FeatureType('UTR3', 'utr3'), FeatureType('UPSTREAM', 'upstream'), FeatureType('DOWNSTREAM', 'downstream')]

# CDS - start is 0-indexed, end is 1-indexed
# UTR5 - start is 0-indexed, end is 0-indexed    
# UPSTREAM - start is 0-indexed, end is 0-indexed
# UTR3 - start is 1-indexed, end is 1-indexed
# DOWNSTREAM - start is 1-indexed, end is 1-indexed
# INTRON - start is 1-indexed, add buffer because end is 0-indexed
            
with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    
    # map vcf variant back onto genome location and report hit class
    for line in list(csv.reader(open(snp_file), delimiter='\t')):
        if line[0][0] == '#':
            continue
        snp_chr = line[0]
        if 'chr' not in snp_chr:
            snp_chr = 'chr' + snp_chr
        if snp_chr == 'chrMT':
            snp_chr = 'chrM'
        snp_start = int(line[1])
        
        for typ in feature_types:
            starts = typ.starts[snp_chr]
            mid = binary_search(starts, snp_start-1)
            output = line.copy()
            
            if typ.name == 'INTRON' and intron_option:
                if snp_start >= starts[mid] and snp_start < starts[mid] + exon_buffer or snp_start <= typ.chrom[snp_chr + '_' + str(starts[mid])] + 1 and snp_start > typ.chrom[snp_chr + '_' + str(starts[mid])] + 1 - exon_buffer:
                    if snp_start - starts[mid] < exon_buffer:
                        output.append('INTRONHIT.' + str(snp_start - starts[mid] + 1) + '.right')
                        output.append(typ.gene[snp_chr + '_' + str(starts[mid])])
                        writer.writerow(output)
                    elif feature_types[1].chrom[snp_chr + '_' + str(starts[mid])] + 1 - snp_start < exon_buffer:
                        output.append('INTRONHIT.' + str(typ.chrom[snp_chr + '_' + str(starts[mid])] - snp_start + 2))
                        output.append(typ.gene[snp_chr + '_' + str(starts[mid])])
                        writer.writerow(output)
            elif snp_start-1 >= starts[mid] and snp_start <= typ.chrom[snp_chr + '_' + str(starts[mid])]:
                output.append(typ.name)
                output.append(typ.gene[snp_chr + '_' + str(starts[mid])])
                writer.writerow(output)
        
        print("Done for SNP " + snp_chr + "_" + str(snp_start) + ".")