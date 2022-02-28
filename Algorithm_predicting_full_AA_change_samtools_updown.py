# algorithm for predicting amino acid changes

# python3 Algorithm_predicting_full_AA_change_samtools_updown.py VCF_input_file_in_tab_delimited_format.vcf.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt > VCF_input_file_in_tab_delimited_format.vcf.out.txt

import sys
from sys import argv
import csv
import re

snp_file = argv[1]
conversion_file = argv[2]
gene_outfile = argv[3]
cds_intron_file = argv[4]
output_file = snp_file + '.out.txt'

# increase field limit
limit = sys.maxsize
while True:
    # decrease limit by factor 10 as long as OverflowError occurs
    try:
        csv.field_size_limit(limit)
        break
    except OverflowError:
        limit = int(limit/10)
        
        
# codon -> AA conversion
aa_dict = {
    'ATG':'M',
    'TGG':'W',
    'TTT':'F',
    'TTC':'F',
    'TAT':'Y',
    'TAC':'Y',
    'TGT':'C',
    'TGC':'C',
    'CAT':'H',
    'CAC':'H',
    'CAA':'Q',
    'CAG':'Q',
    'AAT':'N',
    'AAC':'N',
    'AAA':'K',
    'AAG':'K',
    'GAT':'D',
    'GAC':'D',
    'GAA':'E',
    'GAG':'E',
    'ATT':'I',
    'ATC':'I',
    'ATA':'I',
    'CCT':'P',
    'CCC':'P',
    'CCA':'P',
    'CCG':'P',
    'ACT':'T',
    'ACC':'T',
    'ACA':'T',
    'ACG':'T',
    'GTT':'V',
    'GTC':'V',
    'GTA':'V',
    'GTG':'V',
    'GCT':'A',
    'GCC':'A',
    'GCA':'A',
    'GCG':'A',
    'GGT':'G',
    'GGC':'G',
    'GGA':'G',
    'GGG':'G',
    'TCT':'S',
    'TCC':'S',
    'TCA':'S',
    'TCG':'S',
    'AGT':'S',
    'AGC':'S',
    'TTA':'L',
    'TTG':'L',
    'CTT':'L',
    'CTC':'L',
    'CTA':'L',
    'CTG':'L',
    'CGT':'R',
    'CGC':'R',
    'CGA':'R',
    'CGG':'R',
    'AGA':'R',
    'AGG':'R',
    'TAG':'*',
    'TAA':'*',
    'TGA':'*',
}

# UCSC_ID -> gene symbol
gene_dict = {}
# UCSC_ID -> strand 
strand_dict = {}

for line in list(csv.reader(open(conversion_file), delimiter = '\t')):
    gene_dict[line[0]] = line[4]
for line in list(csv.reader(open(cds_intron_file), delimiter = '\t')):
    strand_dict[line[0]] = line[2]
       
def protein_translation(original_cds_line, cds_line, snp_loc, protein_flag, strand, case_count, cases_left, output):
    codon_change_string = ''
    truncate_original_cds_line = ''
    truncate_original_aa_line = ''
    checked_codon = None

    # translate CDS into AA sequence by moving 3 at a time
    for i in range(0, len(original_cds_line), 3):
        original_codon = original_cds_line[i:i+3]
        truncate_original_cds_line = truncate_original_cds_line + original_codon
        truncate_original_aa_line = truncate_original_aa_line + aa_dict[original_codon]
        
        # if this is the codon position
        if snp_loc - i >= 0 and snp_loc - i <= 2 and protein_flag:
            checked_codon = aa_dict[original_codon]
                                # first case of ALT
            if case_count == 1 or cases_left == case_count: 
                output.extend([strand, i//3 + 1, 'SNP', checked_codon + '(' + original_codon + ')'])
                print(*[strand, i//3 + 1, 'SNP', checked_codon + '(' + original_codon + ')'], end='', sep='\t')
            else:
                output.extend([checked_codon + '(' + original_codon + ')'])
                print(*[checked_codon + '(' + original_codon + ')'], end='', sep='\t')
    
    truncate_cds_line = ''
    truncate_aa_line = ''
    for i in range(0, len(cds_line), 3):
        codon = cds_line[i:i+3]
        truncate_cds_line = truncate_cds_line + codon
        truncate_aa_line = truncate_aa_line + aa_dict[codon]
        
        if snp_loc - i >= 0 and snp_loc - i <= 2 and protein_flag:
            snp_codon_mutant = codon
            snp_aa_mutant = aa_dict[codon]
            
            output[-1] = output[-1] + '->' + snp_aa_mutant + '(' + codon + ')'
            if case_count == 1:
                if checked_codon == snp_aa_mutant:
                    output.append('SYN')
                elif snp_aa_mutant == '*':
                    output.append('NSN')
                else:
                    output.append('NSM')
            else:
                output[-1] = output[-1] + ','
                if checked_codon == snp_aa_mutant:
                    codon_change_string = codon_change_string + 'SYN'
                elif snp_aa_mutant == '*':
                    codon_change_string = codon_change_string + 'NSN'
                else:
                    codon_change_string = codon_change_string + 'NSM'
                if(cases_left > 1): # more ALT cases left
                    codon_change_string = codon_change_string  + ','

    print(truncate_original_cds_line)
    print(truncate_cds_line)
    print(truncate_original_aa_line)
    print(truncate_aa_line)
    
    # newly added for printing full length of AA
    output.append(truncate_original_aa_line)
    output.append(truncate_aa_line)
    
    return codon_change_string        
        
def remove_intron(line):
    return re.sub(r'[acgtn]', '', line)

def rev_dna_comp(dna):
    return dna[::-1].translate(dna.maketrans('ACGTacgt', 'TGCAtgca'))    

cds_dict = {}
with open(gene_outfile, 'r') as in_file:
    for line in in_file:
        line = line.strip()
        if '>' in line:
            line_arr = line.split(' ')
            ucsc_id = line_arr[0][16:]
            if len(line_arr) < 2:
                continue
            arr_chr = line_arr[1].split(':')[0].split('=')
            arr = line_arr[1].split(':')[1].split('-')
            cds_start = int(arr[0])
            cds_end = int(arr[1])
            
            if ucsc_id not in cds_dict:
                cds_dict[ucsc_id] = {}
            cds_dict[ucsc_id][arr_chr[1]] = (cds_start, cds_end, line)
                                    
with open(output_file, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['Sample', 'Chromosome', 'Variant Position', 'Gene Symbol', 'UCSC ID', 'Strand', 'AA Position of Mutation (for CDSHIT)', 'Variant Type', 'Amino Acid Ref (Codon) -> AA SNP (Codon)', 'Variant Class', 'Ref AA chain',  'Alt AA chain',  'Hit Type', 'Known dbSNP', 'Ref nt', 'Alt nt', 'Quality', 'Depth', 'Allele Freq', 'Read Categories', 'Info'])

    # read SNP file
    for line in list(csv.reader(open(snp_file), delimiter = '\t')):
        snp_flag = True # initially assume line is an SNP
        protein_flag = False

        snp_chromosome = line[0]
        snp_loc = int(line[1])
        
        # assume ref has only one case
        ref_char = line[3].upper()
        snp_chars = line[4].upper().split(',')
        case_count = len(snp_chars)
        for s in snp_chars:
            if len(s) > 1: # not a single SNP
                snp_flag = False
                break

        # get depth and read category info
        info_array = line[7].split(';')
        depth_array = []
        alle_freq = []
        read_category = []

        hit_type = None
        ucsc_id = None

        if 'VDB' in line[7]: # samtools v0.1.18
            if len(line) == 14: # three samples run from samtools
                hit_type = line[12]
                ucsc_id = line[13]
                print('Three samples & old samtools...')
            elif len(line) == 13: # two samples run from samtools
                hit_type = line[11]
                ucsc_id = line[12]
                print('Two samples & old samtools...')
            elif len(line) == 12: # one sample run from samtools
                hit_type = line[10]
                ucsc_id = line[11]
                print('One sample & old samtools...')
        else:
#             hit_type = line[10]
#             ucsc_id = line[11]
#             print("One sample & old samtools...")
            print('One sample & new samtools...')
            hit_type = line[8]
            ucsc_id = line[9]

        ucsc_id_flag = False

        strand = None
        codon_change_string = ''
        output = []
        if hit_type == 'CDSHIT':
            print('\n============================')
            print(*line, gene_dict.get(ucsc_id, 'Missing'), sep='\t')
            print('The number of possible ALT cases for this line is ' + str(case_count))

            if len(ref_char) == 1 and snp_flag: # both ref and SNP calls are single SNPs
                protein_flag = True
            output = [snp_file.split('.')[0], line[0][3:], line[1], gene_dict.get(ucsc_id, 'Missing'), ucsc_id]
            
            # for each possible ALT case
            for case in range(0, case_count): 
                if case > 0:
                    print('\n--------------------------------------')
                cds_start = None
                cds_end = None
                line_flag = False 
                
                if ucsc_id in cds_dict:
                    print("found ucsc_id" + ucsc_id)
                    if snp_chromosome in cds_dict[ucsc_id]:
                        print("found snp_chrom" + snp_chromosome)
                        if snp_loc < cds_dict[ucsc_id][snp_chromosome][0] or snp_loc > cds_dict[ucsc_id][snp_chromosome][1]:
                            print('SNP location is outside of CDS region! No checking occurs!')
                            line2 = cds_dict[ucsc_id][snp_chromosome][2]
                            strand = cds_dict[ucsc_id][snp_chromosome][2].split(' ')[4][-1]
                            line_flag = True
                            ucsc_id_flag = True
                                    
                    # found CDS line
                    if line_flag:
                        cds_line = line2
                        line2_no_intron = remove_intron(line2)
                        before_snp_string = None

                        # replace original char with SNP
                        if strand == '+':
                            coordinate = snp_loc - cds_start
                            # count number of lowercase nucleotides (NON-CDS or intron nucleotide) before SNP location and adjust its corrdinate
                            before_snp_string = line2[:coordinate]
                            cds_line = cds_line[:coordinate] + snp_chars[case] + cds_line[coordinate + len(ref_char):]

                            # process new CDS string
                            cds_line_no_intron = remove_intron(cds_line)
                            if case_count == 1:
                                protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(re.findall(r'[a-z]', before_snp_string)), protein_flag, strand, case_count, case_count - case, output)
                            else:
                                codon_change_string = codon_change_string + protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(re.findall(r'[a-z]', before_snp_string)), protein_flag, strand, case_count, case_count - case, output)
                        else:
                            coordinate = cds_end - snp_loc
                            replaced_snp_char = rev_dna_comp(snp_chars[case]).upper()
                            before_snp_string = line2[:coordinate]

                            # process new CDS string
                            cds_line = cds_line[:coordinate - len(ref_char) + 1] + replaced_snp_char + cds_line[coordinate + 1:]
                            cds_line_no_intron = remove_intron(cds_line)
                            if case_count == 1:
                                protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(ref_char) + 1 - len(re.findall(r'[a-z]', before_snp_string)), protein_flag, strand, case_count, case_count - case, output)
                            else:
                                codon_change_string = codon_change_string + protein_translation(line2_no_intron, cds_line_no_intron, coordinate - len(ref_char) + 1 - len(re.findall(r'[a-z]', before_snp_string)), protein_flag, strand, case_count, case_count - case, output)
                        line_flag = False
                        break
            if len(ref_char) == 1 and snp_flag: # single SNP case
                depth_array = info_array[0].split('=')
                if 'VDB' in line[7]: # samtools v0.1.18
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                else: # samtools v0.1.12
                    alle_freq = info_array[1].split('=')
                    read_category = info_array[3].split('=')
                if len(codon_change_string) == 0: # no ALT cases
                    output.extend([hit_type, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
                else:
                    output.extend([codon_change_string, hit_type, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
            else: # indel case
                depth_array = info_array[0].split('=')
                if 'VDB' in line[7]: # samtools v0.1.18
                    alle_freq = info_array[3].split('=')
                    read_category = info_array[5].split('=')
                else: # samtools v0.1.12
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                output.extend([strand, '---', 'INDEL', '---', '---', '---', '---', hit_type, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
            if not ucsc_id_flag:
                print(ucsc_id + ' cannot be found!!')
        else: # INTRON, UTR, or UP-DOWNSTREAM
            output.extend([snp_file.split('.')[0], line[0][3:], line[1], gene_dict.get(ucsc_id, 'Missing'), ucsc_id])

            if len(ref_char) == 1 and snp_flag: # both ref and SNP calls are single SNPs
                depth_array = info_array[0].split('=')
                if 'VDB' in line[7]: # samtools v0.1.18
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                else: # samtools v0.1.12
                    alle_freq = info_array[1].split('=')
                    read_category = info_array[3].split('=')
                output.extend([strand_dict.get(ucsc_id, 'Missing'), '---', 'SNP', '---', '---', '---', '---', hit_type, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
            else:
                depth_array = info_array[1].split('=')
                if 'VDB' in line[7]: # samtools v0.1.18
                    alle_freq = info_array[3].split('=')
                    read_category = info_array[5].split('=')
                else: # samtools v0.1.12
                    alle_freq = info_array[2].split('=')
                    read_category = info_array[4].split('=')
                output.extend([strand_dict.get(ucsc_id, 'Missing'), '---', 'INDEL', '---', '---', '---', '---', hit_type, line[2], line[3], line[4], line[5], depth_array[1], alle_freq[1], read_category[1], line[7]])
        writer.writerow(output)

        
        
        
