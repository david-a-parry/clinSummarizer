 #!/usr/bin/env perl 
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use POSIX qw/strftime/;
use Bio::DB::Sam;
use FindBin qw($RealBin);
use lib "$RealBin/lib/vcfhacks/lib";
use VcfReader;

my $fasta              = "$RealBin/ref_bundle/hs37d5.fa";
my $dbsnp              = "$RealBin/ref_bundle/dbSNP146/All_20151104.vcf.gz";
my $evs                = "$RealBin/ref_bundle/ESP6500SI-V2-SSA137.GRCh38-liftover.combined.vcf.gz";
my $exac               = "$RealBin/ref_bundle/ExAC.r0.3.sites.vep.vcf.gz";
my $email              = "";
my $gatk               = "$RealBin/exe_bundle/GATK/v3.5/GenomeAnalysisTK.jar";
my $vep_dir            = "$RealBin/exe_bundle/variant_effect_predictor";
my $depth_intervals    = "$RealBin/bed_files/capture_regions_grch37.bed";
my $gene_db            = "$RealBin/genes/gene_database.db";
my $reportable_cov     = "$RealBin/bed_files/reportable_coding_exons.bed";
my $gene_list          = "$RealBin/genes/gene_inheritance_and_diseases.txt";
my $not_reportable_cov = "$RealBin/bed_files/not_reportable_coding_exons.bed";
my $tmp_dir            = "$ENV{HOME}/scratch/tmp/";
my $freq               = 0.01;
my @interval_list      = 
(
    "$RealBin/bed_files/reportable_genes_grch37.bed",
    "$RealBin/bed_files/not_reportable_genes_grch37.bed",
);
my $script_dir = "sub_scripts";

my @bams = (); 
my $vcf;

my $threads = 8;
my $mem = 4;
my $vmem = 8;
my $runtime = 4;
my $date = strftime( "%d-%m-%y", localtime );

my $outdir;
my $help;

my %config;
GetOptions(
    \%config,
    "b|bams=s{,}"         => \@bams,
    "d|depth_intervals=s" => \$depth_intervals,
    "e|email=s"           => \$email,
    "evs=s"               => \$evs,
    "exac=s"              => \$exac,
    "f|freq=f"            => \$freq,
    "fasta=s"             => \$fasta,
    "g|gatk=s"            => \$gatk,
    "gene_database=s"     => \$gene_db,
    "gene_list=s"         => \$gene_list,
    "h|help",
    "l|list=s{,}"         => \@interval_list,
    "m|mem=i"             => \$mem,
    "n|print_scripts",
    "o|outdir=s"          => \$outdir,
    "print_scripts",
    "q|qsub",
    "reportable_bed"      => \$reportable_cov,
    "not_reportable_bed"  => \$not_reportable_cov,
    "r|runtime=i"         => \$runtime,
    "s|script_dir=s"      => \$script_dir,
    "t|threads=i"         => \$threads,
    "vmem=i"              => \$vmem,
    "v|vcf=s"               => \$vcf,
    "x|tmp_dir=s"         => \$tmp_dir,
) or die "Syntax error!\n";
usage() if $config{h};

usage("At least one BAM file must be supplied to the -b/--bams argument") 
 if not @bams;

usage("A VCF must be supplied to the -v/--vcf argument") if not $vcf;

$vmem .= "G";
$mem .= "G";
if (not $outdir){
    $outdir = "alignments_and_vcfs_$date";
    if (-d $outdir){
        print STDERR "WARNING: $outdir already exists!\n";
    }
}
if (not (-d $script_dir)){
    mkdir($script_dir) or die "Could not make script directory $script_dir: $!\n";
}

if (not (-d $outdir)){
    mkdir($outdir) or die "Could not make output directory $outdir: $!\n";
}

if (not (-d $tmpdir)){
    mkdir($tmpdir) or die "Could not make tmp directory $tmpdir: $!\n";
}

my $call_gatk = "java -Djava.io.tmpdir=$tmp_dir -Xmx$mem -jar $gatk -R $fasta";

my %vcf_samples = VcfReader::getSamples
(
    vcf => $vcf, 
    get_columns => 1
);

my @bam_samples = getBamSamples();
checkSamples();

my %scripts = ();
make_depth_scripts();
make_summarizer_script();

if ($config{n}){
    print STDERR "Created following scripts:\n" . join("\n", map { @{$scripts{$_}} } sort keys %scripts) ."\n" 
}
if ($config{q}){
    submit_scripts();
}

#################################################
sub make_summarizer_script{
    my $samplesummarizer = "perl $RealBin/sampleSummarizer.pl";
    if (-e "$RealBin/sampleSummarizerEddieBinary"){
       $samplesummarizer = "$RealBin/sampleSummarizerEddieBinary" ;
    }
    my $sscript = "$script_dir/summarizer-$date.sh";
    $sscript = check_exists_and_rename($sscript);
    push @{$scripts{post_depth}}, $sscript;
    open (my $SUBSCRIPT, ">$sscript") or die "Can't open $sscript for writing: $!\n";
    my $sample_list = join(" ", @bam_samples); 
    my $region_vcf = $vcf;
    my $int_cmd; 
    if (@interval_list){
        $region_vcf = fileparse($vcf);
        $region_vcf =~ s/\.vcf(\.gz)*$//;
        my $dir = "$outdir/vcf";
        if (not -d $dir){
            mkdir $dir or die "Could not make filtered VCF directory $dir: $!\n";
        }
        $region_vcf = "$dir/$region_vcf.regions.vcf";
        $int_cmd = "$RealBin/lib/vcfhacks/getVariantsByLocation.pl -i $vcf -b "
         . join (" ", @interval_list) . " -o $region_vcf" ;
    }
    if ($email){
        print $SUBSCRIPT <<EOT
#\$ -M $email
#\$ -m b
#\$ -m e
EOT
;
    }
    print $SUBSCRIPT <<EOT
#\$ -cwd
#\$ -V
#\$ -e $sscript.stderr
#\$ -o $sscript.stdout
#\$ -l h_rt=$runtime:00:00
#\$ -l h_vmem=$vmem
. /etc/profile.d/modules.sh
module load igmm/libs/htslib/1.3
module load igmm/apps/samtools/1.2

$int_cmd

$samplesummarizer -i $region_vcf -t $gene_db  -n $not_reportable_cov -r $reportable_cov -s $dbsnp  -e $evs -x $exac -z $RealBin/ref_bundle/cadd/v1.3/ -q $outdir/fastqc -c $outdir/depth -f $freq -o $outdir/sample_summaries/ -l $gene_list --samples $sample_list

EOT
;
    close $SUBSCRIPT;
}


#################################################
sub make_depth_scripts{
    foreach my $bam (@bams){
        my $f = fileparse($bam); 
        (my $stub = $f) =~ s/\.bam$//;
        my $gscript = "$script_dir/depth-$f.sh";
        $gscript = check_exists_and_rename($gscript);
        my $dir = "$outdir/depth";
        if (not -d $dir){
            mkdir $dir or die "Could not make depth directory $dir: $!\n";
        }
        push @{$scripts{depth}}, $gscript;
        open (my $SCRIPT, ">$gscript") or die "Can't open $gscript for writing: $!\n";
        my $depth = "$dir/depth.$stub";
        if ($email){
            print $SCRIPT <<EOT
#\$ -M $email
#\$ -m a
EOT
;
        }
        print $SCRIPT <<EOT
#\$ -cwd
#\$ -V
#\$ -e $gscript.stderr
#\$ -o $gscript.stdout
#\$ -l h_rt=$runtime:00:00
#\$ -l h_vmem=$vmem
. /etc/profile.d/modules.sh
module load igmm/apps/jdk/1.8.0_66

$call_gatk -T DepthOfCoverage -I $bam -o $depth -ct 4 -ct 9 -ct 19 -ct 29 -ct 49 -ct 99 -ct 199 -mbq 17 -mmq 20  -L $depth_intervals

EOT
;
    close $SCRIPT;
    }
}


##################################################################################
sub submit_scripts{
    my @wait_ids = ();
    foreach my $sc (@{$scripts{depth}}){
        my $cmd = "qsub $sc";
        print STDERR "Executing: $cmd\n";
        my $output = `$cmd`;
        check_exit($?);
        if ($output =~ /Your job (\d+) .* has been submitted/){
            push @wait_ids, $1;
        }else{
            die "Error parsing qsub output for $cmd!\nOutput was: $output";
        }
    }
    
    my $hold = join(",", @wait_ids);
    foreach my $sc (@{$scripts{post_depth}}){
        my $cmd = "qsub -hold_jid $hold $sc";
        print STDERR "Executing: $cmd\n";
        my $output = `$cmd`;
        check_exit($?);
        if ($output =~ /Your job (\d+) .* has been submitted/){
            push @wait_ids, $1;
        }else{
            die "Error parsing qsub output for $cmd!\nOutput was: $output";
        }
    }
}

##################################################################################
sub check_exit{
    my $e = shift;
    if($e == -1) {
        print "Failed to execute: $!\n";
    }elsif($e & 127) {
        printf "Child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
    }elsif($e != 0){
        printf "Child exited with value %d\n", $e >> 8;
    }
}



##################################################################################
sub checkSamples{
    if (not @bam_samples){
        die "ERROR: No samples found in BAM input!\n";
    }
    foreach my $s (@bam_samples){
        if (not exists $vcf_samples{$s}){
            die "ERROR: Sample '$s' does not exist in $vcf!\n";
        }
    }
}

##################################################################################
sub getBamSamples{
    my %samples = (); 
    foreach my $bam (@bams){
        map {$samples{$_} = undef}  getSamplesFromBam($bam);
    }
    return keys %samples; 
}

##################################################################################
sub getSamplesFromBam{
    my $f = shift;
    my @samples = ();
    my $bam  = Bio::DB::Bam->open($f);
    my $header = $bam->header->text;
    my @rgs = grep {/^\@RG/} split ("\n", $header);
    foreach my $rg (@rgs){
        if ($rg =~ /SM:(\S+)/){#get read group sample name
            push @samples, $1;
        }else{
            die "Could not parse following read group from $f:\n$rg\n";
        }
    }
    return @samples;
}

##################################################################################
sub check_exists_and_rename{
    my $f = shift;
    if (-e $f){
        warn "WARNING: $f already exists!\n";
        (my $newname = $f) =~ s/(\.\w+)$//;
        my $ext = $1;
        $newname = $f if not $newname;
        my $add = 1;
        if ($newname =~ /_(\d+)$/){
            $add = $1 + 1;
            $newname =~ s/_(\d+)$//;
        }
        $newname .= "_$add";
        $newname .= "$ext" if $ext;
        return check_exists_and_rename($newname);
    }
    return $f;
}


##################################################################################
sub usage{
    my $msg = shift; 
    print "ERROR: $msg\n" if $msg;
    print <<EOT
    
    USAGE: perl $0 -b [bam1 bam2 ...] -v [vcf] [options]
    
    From given BAM files and a VEP annotated VCF this script will create a series of qsub scripts to analyze a gene panel. 
    
    
    All samples from given BAM files must be present in the provided VCF. The VCF must be annotated with Ensembl's Variant Effect Predictor using the LoF and MaxEntScan plugins.     

    OPTIONS: 
    
    -b    --bams               [one or more input BAM files]
    -v    --vcf                [input multisample  VCF file, VEP annotated using 
                                LoF and MaxEntScan plugins]
    -o    --outdir             [optional location for output files 
                                - default will be a new subdirectory called alignments_and_vcfs_$date];
    -s    --script_dir         [directory to write sub scripts to - default = sub_scripts]
    -q    --qsub               [submit created scripts once finished - this automatically enables queueing of the post-alignment commands for convenience]
    -f    --freq               [allele frequency cutoff for reporting variants in final summary document]
    -d    --depth_intervals    [bed file for DepthOfCoverage 
                                - default = $RealBin/bed_files/capture_regions_grch37.bed]
    -n    --print_scripts      [print names of created scripts once finished]
    -g    --gatk               [location of GATK jar file 
                                - default = $RealBin/exe_bundle/GATK/v3.5/GenomeAnalysisTK.jar]
    -l    --list               [interval list(s) to use with GATK commands 
                                - default = $RealBin/bed_files/reportable_genes_grch37.bed and 
                                 $RealBin/bed_files/not_reportable_genes_grch37.bed]
    -e    --email              [address to email script messages to]
    -m    --mem                [amount of memory (in GB) for each GATK/Picard command - default = 4]
          --vmem               [amount of vmem for each script - default = 8]
    -t    --threads            [no. threads for each command - default = 8]
    -r    --runtime            [runtime in hours for each script - default = 4] 
    -x    --tmp_dir            [directory to use for temporary storage for commands 
                                - default =  $ENV{HOME}/scratch/tmp/]
          --dbsnp              [dbSNP file 
                                - default = $RealBin/ref_bundle/dbSNP146/All_20151104.vcf.gz]
          --evs                [EVS VCF file (all chromosomes in one file) 
                                - default = $RealBin/ref_bundle/ESP6500SI-V2-SSA137.GRCh38-liftover.combined.vcf.gz]
          --exac               [ExAC VCF file 
                                - default = $RealBin/ref_bundle/ExAC.r0.3.sites.vep.vcf.gz]
          --fasta              [genome reference fasta sequence 
                                - default = $RealBin/bed_files/PD_capture_capture_v4_1_Regions_b37.bed
          --gene_list          [TSV file of gene symbol, expected inheritance pattern and associated conditions for annotating sampleSummarizer output 
                                - default = $RealBin/genes/gene_inheritance_and_diseases.txt]
          --gene_database      [SQLITE database created using dbCreator.pl for all genes. This is used for annotating sampleSummarizer output 
                                - default = $RealBin/genes/gene_database.db]
          --reportable_cov     [Bed file of regions in reportable genes to analyze coverage in (e.g. coding exons) 
                                - default = "$RealBin/bed_files/reportable_coding_exons.bed]
          --not_reportable_cov [Bed file of regions in non-reportable genes to analyze coverage in 
                                - default = "$RealBin/bed_files/not_reportable_coding_exons.bed]
    -h    --help               [show this help message and exit]

EOT
;
    exit 1 if $msg;
    exit;
}
