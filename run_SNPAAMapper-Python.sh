#!/bin/bash

. $1

# ./run_SNPAAMapper-Python.sh config.txt 

# generate exon annotation file
START1=$(date +%s)
# pytho3n Algorithm_generating_annotation_exon.py $geneAnnotation
END1=$(date +%s)
DIFF1=$(( $END1 - $START1))
echo "It took $DIFF1 seconds to generate exon annotation file!"

# process exon annotation files and generate feature start and gene mapping files
START2=$(date +%s)
# python3 Algorithm_preprocessing_exon_annotation_RR.py "$geneAnnotation".exon
python3 Algorithm_preprocessing_exon_annotation_RR.py "$geneAnnotation".exon
END2=$(date +%s)
DIFF2=$(( $END2 - $START2))
echo "It took $DIFF2 seconds to generate feature start and gene mapping files!"

# classify variants by regions (CDS, Upstream, Downstream, Intron, UTRs...)
START3=$(date +%s)
# python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf
# OR
# python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf IntronExon_boundary_in_bp
python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py "$geneAnnotation".exon $vcfFile $intronBoundary
END3=$(date +%s)
DIFF3=$(( $END3 - $START3))
echo "It took $DIFF3 seconds to classify variants by regions!"


# predict amino acid change type
START4=$(date +%s)
# python3 Algorithm_predicting_full_AA_change_samtools_updown.py VCF_input_file_in_tab_delimited_format.vcf.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt > VCF_input_file_in_tab_delimited_format.txt.out.txt
python3 Algorithm_predicting_full_AA_change_samtools_updown.py "$vcfFile".append $conversionFile $sequenceFile $geneAnnotation > "$vcfFile".out.txt
END4=$(date +%s)
DIFF4=$(( $END4 - $START4))
echo "It took $DIFF4 seconds to predict amino acid change type!"

# prioritize mutation effects
START5=$(date +%s)
# python3 Algorithm_prioritizing_mutation_headerTop_updown.py VCF_input_file_in_tab_delimited_format.vcf.append.out.txt
python3 Algorithm_prioritizing_mutation_headerTop_updown.py "$vcfFile".append.out.txt
END5=$(date +%s)
DIFF5=$(( $END5 - $START5))
echo "It took $DIFF5 seconds to prioritize mutation effects!"