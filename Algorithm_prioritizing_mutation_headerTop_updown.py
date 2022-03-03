# algorithm for prioritizing mutation effects

# python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf
# python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf IntronExon_boundary_in_bp


import sys
from sys import argv
import csv

input_file = argv[1]
output_file = input_file + ".prioritized_out"

# increase field limit
limit = sys.maxsize
while True:
    # decrease limit by factor 10 as long as OverflowError occurs
    try:
        csv.field_size_limit(limit)
        break
    except OverflowError:
        limit = int(limit/10)
        

with open(input_file, mode='r') as line:
    line_array = line

data = list(csv.reader(open(input_file), delimiter='\t'))
columns = data[0]
data = data[1:]

score = []
output_row = {} # stores highest-scoring row number for each variant

for i in range(len(data)):
    line = data[i]
    score.append(-1)
    if line[12] == 'CDSHIT':
        # get rid of right reciprocal pairs
        if line[7] == 'SNP':
            for s in line[9].split(','):
                if s == 'NSM':
                    score[i] = 6
                if s == 'NSN':
                    score[i] = 5
                if s == 'SYN':
                    score[i] = 4
        else:
            score[i] = 3.5
    elif line[12] == 'UPSTREAMHIT':
        score[i] = 3
    elif line[12] == 'UTR3HIT':
        score[i] = 2
    elif line[12] == 'UTR5HIT':
        score[i] = 1
    elif line[12] == 'DOWNSTREAMHIT':
        score[i] = 0
    
    # use combination key: sample_chromosome_gene
    key = "_".join(line[0:3])
    
    # save row number with the highest score
    if key not in output_row or score[i] > score[output_row[key]]:
        output_row[key] = i

with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(columns)
    for i in output_row:
        writer.writerow(data[output_row[i]])