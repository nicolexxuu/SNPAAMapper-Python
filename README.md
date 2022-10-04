# SNPAAMapper-Python

Ma, K., N. Xu, A. He, and Y. Bai, “SNPAAMapper-Python: A highly efficient genome-wide SNP variant analysis pipeline for next-generation sequencing data.” Poster Presentation at the 18th Annual Conference for the Mid-South Computational Biology and Bioinformatics Society (MCBIOS 2022)

SNPAAMapper is a downstream variant annotation program that can effectively classify variants by region (e.g. exon, intron, etc.), predict amino acid change type (e.g. synonymous, non-synonymous mutation, etc.), and prioritize mutation effects (e.g. CDS versus 5'UTR, etc.).

## Features

- The pipeline accepts a VCF input file in tab-delimited format and processes the vcf input file containing all cases (G5, lowFreq, and novel)
- The variant mapping step allows users to select whether they want to report the base pair distance between each identified intron variant and its nearby exon
- Compatibility with VCF files called by different SAMTools versions (0.1.18 and older) and/or generated using SAMTools with two or three samples
-  The spreadsheet result file contains full protein sequences for both reference and alternative alleles, which makes it easier for downstream protein structure/function analysis tools to use

## Requirements

- python 3.x
- sys
- csv
- re
- shutil
- [Git LFS](https://git-lfs.github.com/)

## Instructions

If you haven't yet, initialize Git LFS by running 

```sh
git lfs install
```

Clone this repo as follows

```sh
git clone https://github.com/nicolexxuu/SNPAAMapper-Python
cd ./SNPAAMapper-Python
```

and download [hg19_CDSIntronWithSign.txt.out](https://drive.google.com/file/d/1n1qI-DvhyZhqVWl9YuyGGW91b0IZLFfe/view?usp=sharing) to your local repository.

Next, type

```sh
./run_SNPAAMapper-Python.sh config.txt
```

OR run the following steps in sequential order (Note: the first two steps were compiled for the human hg19 genome and output files have already been generated):

<!-- 1. Generate annotation file:

    ```sh
    python3 Algorithm_generating_annotation_exon.py ChrAll_knownGene.txt
    ```
    -->
2. Process exon annotation files and generate feature start and gene mapping files:

    ```sh
    python3 Algorithm_preprocessing_exon_annotation_RR.py ChrAll_knownGene.txt.exon
    ```
    
3. Classify variants by regions (CDS, Upstream, Downstream Intron, UTRs...)

    ```sh
    python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf
    ```
    
    OR
    
    ```sh
    python3 Algorithm_mapping_variants_reporting_class_intronLocation_updown.py ChrAll_knownGene.txt.exon VCF_input_file_in_tab_delimited_format.vcf IntronExon_boundary_in_bp
    ```
    
4. Predict amino acid change type

    ```sh
    python3 Algorithm_predicting_full_AA_change_samtools_updown.py VCF_input_file_in_tab_delimited_format.vcf.append kgXref.txt hg19_CDSIntronWithSign.txt.out ChrAll_knownGene.txt >VCF_input_file_in_tab_delimited_format.vcf.out.txt
    ```
    
5. Prioritize mutation effects

    ```sh
    python3 Algorithm_prioritizing_mutation_headerTop_updown.py VCF_input_file_in_tab_delimited_format.vcf.append.out.txt
    ```

***The final output file is \*.append.out.txt.prioritzed_out.***

## References
1. “The Human Genome Project.” *Genome.gov*, www.genome.gov/human-genome-project.
2. *Nature News*, Nature Publishing Group, www.nature.com/articles/d42473-021-00030-9.
3. Lewis, Tanya. “Human Genome Project Marks 10th Anniversary.” *LiveScience*, Purch, 14 Apr. 2013, www.livescience.com/28708-human-genome-project-anniversary.html. 
4. Barba, Marina, Czosnek, Henryk, Hadidi, Ahmed. “Historical Perspective, Development and Applications of next-Generation Sequencing in Plant Virology.” *Viruses*, MDPI, 6 Jan. 2014, www.ncbi.nlm.nih.gov/pmc/articles/PMC3917434/.
5. Bai, Yongsheng, and James Cavalcoli. “SNPAAMapper: An Efficient Genome-Wide SNP Variant Analysis Pipeline for next-Generation Sequencing Data.” Bioinformation, Biomedical Informatics, 16 Oct. 2013, www.ncbi.nlm.nih.gov/pmc/articles/PMC3819573/.
6. “UCSC Genome Browser Project History.” *Genome Browser History*, https://genome.ucsc.edu/goldenPath/history.html.
7. “The Perl Programming Language.” *TIOBE*, https://www.tiobe.c

## License

MIT
