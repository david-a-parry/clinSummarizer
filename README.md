# clinSummarizer
Wrappers and scripts for targetted sequencing analysis

##Getting started

First we need to create our reference files.

You will need a BED file of your capture regions. This will be used for 
calculating depth of coverage from your BAM files. If you don't want to specify 
the name of this bed file for each of your runs place it in the 'bed_files' 
subdirectory and call it 'capture_regions_grch37.bed'.

You will also need a tab delimited file with four columns indicating the names 
of your target gene names, the associated inheritance pattern, whether they are 
in the reportable or non-reportable category and the associated medical 
condition. For an example see 'genes/example_gene_inheritance_and_diseases.txt' 
in this directory. The default name for this file is 
'gene_inheritance_and_diseases.txt' and it should be placed in the 'genes' 
subdirectory.

Once you have this information you can run the 'updateBedFiles.pl' script to 
create BED files for the target genes specified in your 
'gene_inheritance_and_diseases.txt' file.

     ./updateBedFiles.pl

 Next you need to create the SQLITE gene database for your genes of interest. 

#TODO... 
