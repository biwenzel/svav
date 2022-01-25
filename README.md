# SVAV
Strucure variants assembly and validation, a tool to evaluate and validate structure variants in diploid organisms using local assembly.
At this point it is only usable for deletions and insertions.

![vertical](https://user-images.githubusercontent.com/97948065/149920073-85aa78b3-88bd-4602-aa30-b45489d9925c.png)

#### USAGE
python svav.py -o output_dir SV_coordinates.csv -m chr haplotyped.bam > logfile.log

SV_coordinates.csv:
- the file has to be a csv, which includes the following columns:

tech,type,length,pos1,pos2,Chromosome

svim,INS,72,41263994,41263994,chr1
- a indicator column is not mandatory, but can be used

haplotyped.bam
- the haplotypes have to be tagged with the label HP with either 1 or 2

#### SCRIPTS
svav.py
- main script
- combines all steps of the workflow
- takes the coordinates of the variant and a (phased) BAM file
- the default reference is hg38
- different mapping options are possible


clustering_and_creating.py
- takes the coordinates and the BAM file
- extracting local reads
- tries to create two FASTA files with at least 10 reads per haplotype
- if this fails, it tries to cluster the reads
- if this fails, it's presumed to be homozygous


calling_sv.py 
- takes the newly created BAM file it searches each read (consensus sequence) for structure variants
- creates a table of found variants

compare_sv.py
- takes two variant informations
- compares the variants based on their type, size, and position
- creates a table for validated and not validated variants


#### REQUIREMENTS:
Python packages:
biopython
numpy
scikit-learn
matplotlib
pysam
os
subprocess
sys
pandas
time
argparse

Tools in PATH:
wtdbg2
sed
minimap2
samtools
