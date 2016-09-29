# clinSummarizer
Wrappers and scripts for targetted sequencing analysis

##Installation

Clone this repository as follows:

    git clone --recursive https://github.com/gantzgraf/clinSummarizer.git

These programs are written in Perl and require several non-core modules. 
You are likely to get an error when attempting to run these scripts for the
first time, telling you about missing modules. Please see
http://www.cpan.org/modules/INSTALL.html for instructions on how to install
these modules. Below is a list of these non-core modules that you are likely to
need install:

    Term::ProgressBar
    Bio::DB::HTS::Tabix
    Excel::Writer::XLSX 
    HTTP::Tiny
    LWP::Simple
    JSON

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
The dbCreator.pl program will read provided gene names and identify the 
transcript IDs, uniprot IDs, search NCBI's conserved domain database for domain
features and optionally add Evolutionary Trace (http://mammoth.bcm.tmc.edu/ETserver.html), 
ClinVar and HGMD information. In order to add ClinVar and HGMD features you will
need to provide VCF files annotated with Ensembl's Variant Effect Predictor
(VEP). ClinVar VCFs can be downloaded from the ClinVar FTP site at NCBI. HGMD
VCFs can be created from HGMD Professional's MART TSV output using the
hgmdMartToVcf.pl script from VcfHacks (provided in the lib/vcfhacks/
subdirectory). For Evolutionary Trace information you will need the prescored
files for Human RefSeq protein's available from
http://mammoth.bcm.tmc.edu/ETserver.html.

    ./dbCreator.pl -d genes/gene_database.db \
    -l genes/gene_inheritance_and_diseases.txt
    -c ~/ref/clinvar/vep.clinvar.vcf.gz \
    -m ~/ref/hgmd/vep.HGMD_variants_converted.vcf \
    -e ~/ref/HumanProteinsEA/data/ 
       
Occasionally the final CDD searches will fail. If the previous stages have 
completed successfully you can rerun just these steps, where protein IDs will be
read from the gene database as follows:
    
    ./dbCreator.pl -d genes/gene_database.db -n


###When you need to add more genes to your database

If at a future date you need to add more genes to your database you will need to
update the 'genes/gene_inheritance_and_diseases.txt' diseases file appropriately
and run './updateBedFiles.pl' again. You will also need to add information for 
these genes to your gene database file. The quickest way to do this is to run:

    ./dbCreator.pl -d genes/gene_database.db \
    -i newGene1[newGene2 ... newGeneN]  
    -c ~/ref/clinvar/vep.clinvar.vcf.gz \
    -m ~/ref/hgmd/vep.new_gene_HGMD_variants_converted.vcf \
    -e ~/ref/HumanProteinsEA/data/ 

Alternatively, you may rerun the command using your modified 
'genes/gene_inheritance_and_diseases.txt' file, but this will be slower due to 
dbCreator.pl looking up information for all the genes in the database.

##Running an analysis

The sampleSummarizer.pl script analyzes data from a multisample VCF and 
optionally depth of coverage data and FASTQC reports for given data. It produces 
per-sample and (optionally) an aggregated report of variant data. It may either 
be run as a stand-alone program on pre-existing data, or as the final step of 
the makeAlignmentCommands.pl wrapper program.

###From scratch using makeAlignmentCommands.pl

The makeAlignmentCommands.pl program is designed to work on a server running 
GridEngine job submission software. It creates and submits a series of qsub 
scripts to parallelize FASTQC, FASTQ trimming, alignment, depth of coverage 
analysis and variant calling. The starting point is a set of FASTQ files named 
in the format:

    [samplename]_S[0-9]_L[lanenumber]_R[readnumber]_[numbers].fastq(.gz)
 or
    [samplename]_[tenLetterBarcode]_L[lanenumber]_R[readnumber]_[numbers].fastq(.gz)

This wrapper also needs to know where to access various executable files. The 
required files are listed in the help documentation for the script. 

Currently this is hardcoded to work with the Eddie3 at the University of
Edinburgh. Modules loaded by the various qsub scripts are those already
available on this cluster. In order to run on a different cluster the 'module
load' commands of the various qsub scripts will require editing.

###From a pre-existing VCF

For a VCF to be compatible with the sampleSummarizer.pl it must be annotated
with Ensembl's VEP using the LoF and MaxEntScan plugins. In order to use
frequency filtering for variants present in dbSNP, EVS or ExAC the relevant VCF
files must be provided when running the program. Similarly, if you want to CADD
score your output you must provide a directory containing pre-scored variants.
The rest of the information required by the script will be read from your gene
database file and your 'gene_inheritance_and_diseases.txt' file.

If you want sheets with details of depth of coverage summaries for target gene
exons like those provided by the makeAlignmentCommands.pl program, you will need
to provide both the relevant BED files and a directory containing coverage data
for each sample. The sampleSummarizer.pl program expects to find the per-base
and sample_summary coverage files named in the format:

    'depth.[samplename][.-_].*bqsr(.sample_summary)*'.

An example command is given below:

    ./sampleSummarizer.pl -t  ../genes/clinsummarizer_database_Aug2016.db \
    -l genes/gene_inheritance_and_diseases.txt
    -f 0.01 \
    -s ~/ref/GRCh37/dbSNP147/All_20160601.vcf.gz \
    -x ~/ref/GRCh37/exac/ExAC.r0.3.sites.vep.vcf.gz \
    -e ~/ref/GRCh37/ESP6500SI-V2-SSA137.GRCh38-liftover.combined.vcf.gz \
    -z ~/ref/GRCh37/cadd/v1.3/ \
    -i vep.var.panel.filters.vcf \
    -o sample_summaries/ \
    -u sample_summaries/panel_summary.xlsx \
    [-q fastqc_output  -c depth_output/  -r reportable_exons.bed \
    -n not_reportable_exons.bed --progress ]


##AUTHOR

David A. Parry
University of Edinburgh

##COPYRIGHT AND LICENSE

Copyright 2016  David A. Parry

These programs are free software: you can redistribute them and/or modify
them under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version. These programs are distributed in the hope that they will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
License for more details. You should have received a copy of the GNU General
Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
