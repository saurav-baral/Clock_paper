# here are details on how to run a GSEA

# at the end is a script to take in two files that have been pre processed

# 1) you need an annotation file that has a list of all genes tab then
# comma separated GO terms (this can be geneated from an EggNOGG annotation file using
# cut -f1,6 (this is shown below)
# 2) you need a list of the candidate genes from that list you want to assess
# for being enriched compared to the genome.
#
# NOTE: these two files must have the following headers added to them.
annotation=japonicus_aa.fasta.emapper.annotations
candidate_file=japponicus_run6above95
cut -f1,6 $annotation > $annotation.tsv
echo 'Parent  GO_term' | cat - "$annotation".tsv > "$annotation".header.tsv
echo 'geneid' | cat - "$candidate_file" > "$candidate_file".header.tsv
#
# the resulting *header.tsv files are then used in the script below
# 1st argument is the annotation file
# 2nd argument is the candidate set list
#
# note that the current setting for GO node size is 5, which
# is hard coded below for more robust results.
#
# example run
Rscript GSEA_run_script.R gifuensis_aa.fasta.emapper.annotations.header.tsv gifuensis_run6above95.header.tsv

# descriptions
# http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html

# note that annotations from EggNOGG need to be filtered at a given "bit score" or evalue.
# for example here I'm looking at the distribution of bitscores across all annotations
cut -f4 MM_90qfcby5.emapper.annotations.tsv | head
##
score
133.0
134.0
704.0
635.0
84.0

# so there are some bitscores that are very low, and thus these annotatios shoudl not be used as they are
# likely false positives, etc.

#############
# generating our GO input file

wc -l MM_90qfcby5.emapper.annotations.tsv
# 12653

# filtering
awk -F"\t" '{if ($4 > 100 ) print $0}' MM_90qfcby5.emapper.annotations.tsv | wc -l
# 12004

cut -f1,10 MM_90qfcby5.emapper.annotations.tsv | head
#query	GOs
ENSPXYG00000027403	-
ENSPXYG00000028407	-
ENSPXYG00000029176	-
ENSPXYG00000028850	GO:0001816,GO:0002250,GO:0002252,GO:0002253,GO:0002367,GO:0002369,GO:0002376,GO:0002429,GO:0002440,GO:0002443,GO:0002449,GO:0002456,GO:0002460,GO

###
# generating our GO dataset
awk -F"\t" '{if ($4 > 100 ) print $0}' MM_90qfcby5.emapper.annotations.tsv | cut -f1,10 - | grep -v '#' > MM_90qfcby5.GOfiltered.tsv
#
head MM_90qfcby5.GOfiltered.tsv
ENSPXYG00000027403	-
ENSPXYG00000028407	-
ENSPXYG00000029176	-
ENSPXYG00000028850	GO:0001816,GO:0002250,GO:0002252,GO:0002253,GO:0002367,GO:0002369,GO:0002376,GO:0002429,GO:0002440,GO:0002443,GO:0002449,GO:0002456,GO:0002460,GO

#
echo 'Parent  GO_term' | cat - MM_90qfcby5.GOfiltered.tsv > MM_90qfcby5.GOfiltered.header.tsv
head MM_90qfcby5.GOfiltered.header.tsv
Parent  GO_term
ENSPXYG00000027403	-
ENSPXYG00000028407	-
ENSPXYG00000029176	-
ENSPXYG00000028850	GO:0001816,GO:0002250,GO:0002252,GO:0002253,GO:0002367,GO:0002369,GO:0002376,GO:0002429,GO:0002440,GO:0002443,GO:0002449,GO:0002456,GO:0002460,GO:0002682,GO:0002684,GO:0002757,GO:0002764,GO:0002768,GO:0003674,GO:0005085,GO:0005088,GO:0005488,GO:0005515,GO:0005543,GO:000557

#############
# candidate gene list
head candidate_genes_table.header.tsv
# this file needs to have as the first row "geneid"
# if needed, use echo 'geneid' | cat - "$candidate_file" > "$candidate_file".header.tsv
geneid
ENSPXYG00000000083
ENSPXYG00000000141
ENSPXYG00000000183
ENSPXYG00000000268
ENSPXYG00000000552
ENSPXYG00000000929



# now run R script for these two files. This will run both the BP and MF GO term analysis
# it will also generate a file for REVIGO, as well as a TSV file of the raw output, and an xlsx version
# in the table, it will have results from several different topGO analysis. Read the

Rscript GSEA_script.R MM_90qfcby5.GOfiltered.header.tsv candidate_genes_table.header.tsv


########
# Example analysis
########
#
# GSEA on the internal branches
# what you need is the anntation file for the correct species
# cd /cerberus/projects/chrwhe/gobo/GSEA/gsea_gobo
# sudo cat ../internal_run6above95 > internal_run6above95
#
# mkdir internal_test
# cd internal_test
# # header addition for each and their annotation files
# ln -s ../gifuensis_aa.fasta.emapper.annotations .
# ln -s ../internal_run6above95 .
# annotation=gifuensis_aa.fasta.emapper.annotations
# candidate_file=internal_run6above95
# cut -f1,6 $annotation > $annotation.tsv
# echo 'Parent  GO_term' | cat - "$annotation".tsv > "$annotation".header.tsv
# echo 'geneid' | cat - "$candidate_file" > "$candidate_file".header.tsv
#
#
# Rscript /data/programs/scripts/GSEA_run_script.R gifuensis_aa.fasta.emapper.annotations.header.tsv internal_run6above95.header.tsv

# or
annotation=Colias_crocea_editTrinity.okay.aa.emapper.annotations
candidate_file=top300_abdomen_no134_103_color_results
cut -f1,6 $annotation > $annotation.tsv
echo 'Parent  GO_term' | cat - "$annotation".tsv > "$annotation".header.tsv
echo 'geneid' | cat - "$candidate_file" > "$candidate_file".header.tsv

Rscript /data/programs/scripts/GSEA_run_script.R $annotation.header.tsv $candidate_file.header.tsv


zip internal_run6above95.zip *.GSEA_*

scp chrwhe@duke.zoologi.su.se:/cerberus/projects/chrwhe/gobo/GSEA/gsea_gobo/internal_test/internal_run6above95.zip .
