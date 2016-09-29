#!/usr/bin/env perl 
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Data::Dumper;
use POSIX qw/strftime/;
use FindBin qw($RealBin);

my $fasta              = "$RealBin/ref_bundle/hs37d5.fa";
my $mills              = "$RealBin/ref_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf";
my $indels             = "$RealBin/ref_bundle/1000G_phase1.indels.b37.vcf";
my $dbsnp              = "$RealBin/ref_bundle/dbSNP146/All_20151104.vcf.gz";
my $evs                = "$RealBin/ref_bundle/ESP6500SI-V2-SSA137.GRCh38-liftover.combined.vcf.gz";
my $exac               = "$RealBin/ref_bundle/ExAC.r0.3.sites.vep.vcf.gz";
my $cadd               = "$RealBin/ref_bundle/cadd/v1.3/";
my $email              = "";
my $gatk               = "$RealBin/exe_bundle/GATK/v3.5/GenomeAnalysisTK.jar";
my $picard             = "$RealBin/exe_bundle/picard/picard-tools-2.1.1/picard.jar";
my $vep_dir            = "$RealBin/exe_bundle/variant_effect_predictor";
my $maxentscan         = "$RealBin/exe_bundle/maxentscan/fordownload/";
my $depth_intervals    = "$RealBin/bed_files/capture_regions_GRCh37.bed";
my $reportable_cov     = "$RealBin/bed_files/reportable_coding_exons_GRCh37.bed";
my $gene_list          = "$RealBin/genes/gene_inheritance_and_diseases.txt";
my $gene_db            = "$RealBin/genes/gene_database.db";
my $not_reportable_cov = "$RealBin/bed_files/not_reportable_coding_exons_GRCh37.bed";
my $tmp_dir            = "$ENV{HOME}/scratch/tmp/";
my $freq               = 0.01;
my @interval_list      = 
(
    "$RealBin/bed_files/reportable_genes_GRCh37.bed",
    "$RealBin/bed_files/not_reportable_genes_GRCh37.bed",
);
my $script_dir = "sub_scripts";

my $threads = 8;
my $mem = 4;
my $vmem = 8;
my $runtime = 4;
my $date = strftime( "%d-%m-%y", localtime );
my $vcf_name = 'panel';
my $outdir;
my $help;

my %config;
GetOptions(
    \%config,
    "c|cadd_dir=s"        => \$cadd,
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
    "i|indels=s"          => \$indels,
    "l|list=s{,}"         => \@interval_list,
    "maxenstscan=s"       => \$maxentscan,
    "m|mem=i"             => \$mem,
    "mills=s"             => \$mills,
    "n|print_scripts",
    "o|outdir=s"          => \$outdir,
    "print_scripts",
    "p|picard=s"          => \$picard,
    "q|qsub",
    "reportable_bed"      => \$reportable_cov,
    "not_reportable_bed"  => \$not_reportable_cov,
    "r|runtime=i"         => \$runtime,
    "s|script_dir=s"      => \$script_dir,
    "t|threads=i"         => \$threads,
    "v|vmem=i"            => \$vmem,
    "vcf_name=s"          => \$vcf_name,
    "vep_dir=s"           => \$vep_dir,
    "x|tmp_dir=s"         => \$tmp_dir,
) or die "Syntax error!\n";
usage() if $config{h};

usage("A file list is required!") if not @ARGV ;

$vmem .= "G";
$mem .= "G";
my $interval_string = "";
foreach my $int (@interval_list){
    $interval_string .= "-L \"$int\" ";
}

if (not $outdir){
    $outdir = "alignments_and_vcfs_$date";
    if (-d $outdir){
        print STDERR "WARNING: $outdir already exists!\n";
    }
}

if (not (-d $tmp_dir)){
    mkdir($tmp_dir) or die "Could not make tmp directory $tmp_dir: $!\n";
}

if (not (-d $script_dir)){
    mkdir($script_dir) or die "Could not make script directory $script_dir: $!\n";
}

if (not (-d $outdir)){
    mkdir($outdir) or die "Could not make output directory $outdir: $!\n";
}

my %fastqs = ();
foreach my $file (@ARGV){
    chomp($file);
    my ($sample, $snum, $lane, $read); 
    my ($f, $d) = fileparse($file);
    if ($f !~ /\.fastq(\.gz)$/){
        print STDERR "$d$f does not look like a FASTQ file - skipping.\n";
        next;
    }
=cut
    if ($d =~ /Sample_(\w+)([-_]\w+)*/){
        $sample = $1;
        $sample .= $2 if $2;
    }elsif ($f =~ /^([\w\-_]+)_L00/){
        $sample = $1;
        $sample .= $2 if $2;
=cut
    if ($f =~ /^([\w\-_]+)_(S\d+|[ACTG]{10})_L00/){
        $sample = $1;
        $snum = $2;
    }else{
        print STDERR "WARNING: Couldn't process sample name for $d$f - skipping.\n";
        next;
    }
    if ($f =~ /\S+_L(\d{1,3})_R(\d)_\d+/){
        $lane = $1;
        $read = $2;
        #add to hash
    }else{
        print STDERR "WARNING: Couldn't process lane and read details for $d$f - skipping.\n";
        next; 
    }
    if ($sample =~ /undetermined/i){
        print STDERR "INFO: Skipping undetermined file: $d$f\n";
        next;
    }
    $sample = check_sample_name($sample, $snum, $lane, $read, $d, $f); 
    $fastqs{$sample}->{$snum}->{$lane}->{$read} = "$d$f";
}
die "No fastqs to process!\n" if not keys %fastqs;
my $call_gatk = "java -Djava.io.tmpdir=$tmp_dir -Xmx$mem -jar $gatk -R $fasta";
my @bams = ();
my @gvcfs = ();
my %scripts = ();

make_fastqc_scripts();
make_alignment_scripts();
make_genotype_script();
make_depth_scripts();

if ($config{n}){
    print STDERR "Created following scripts:\n" . join("\n", map { @{$scripts{$_}} } sort keys %scripts) ."\n" 
}
if ($config{q}){
    submit_scripts();
}

check_duplicates();



#################################################
sub make_fastqc_scripts{
    my $dir = "$outdir/fastqc";
    if (not -d $dir){
        mkdir $dir or die "Could not make fastqc directory $dir: $!\n";
    }
    foreach my $s (sort keys %fastqs){
        foreach my $d (sort keys %{$fastqs{$s}}){
            foreach my $l (sort keys %{$fastqs{$s}->{$d}}){
                foreach my $r (sort  keys %{$fastqs{$s}->{$d}->{$l}}){
                    my $fq = $fastqs{$s}->{$d}->{$l}->{$r};
                    my $f = fileparse($fq);
                    my $script = "$script_dir/fastqc-$f.sh";
                    $script = check_exists_and_rename($script);
                    push @{$scripts{fastqc}}, $script;
                    open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
                    if ($email){
                        print $SCRIPT <<EOT
#\$ -M $email
#\$ -m b
#\$ -m e
EOT
;
                    }
                    print $SCRIPT <<EOT

#\$ -cwd
#\$ -V
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -l h_rt=$runtime:00:00
#\$ -l h_vmem=$vmem
. /etc/profile.d/modules.sh
module load igmm/apps/FastQC/0.11.4

fastqc -d $tmp_dir --extract -o $dir $fq

EOT
;
                    close $SCRIPT;
                }
            }
        }
    }
}
#################################################
sub make_alignment_scripts{
    my $dir = "$outdir/alignments";
    if (not -d $dir){
        mkdir $dir or die "Could not make alignments directory $dir: $!\n";
    }
    my $trim_dir = "$outdir/trimmed_fastq";
    if (not -d $trim_dir){
        mkdir $trim_dir or die "Could not make trimmed fastq directory $trim_dir: $!\n";
    }
    foreach my $s (sort keys %fastqs){
        foreach my $d (sort keys %{$fastqs{$s}}){
            foreach my $l (sort keys %{$fastqs{$s}->{$d}}){
                my @reads = ();
                my @trimmed = ();
                my @intermediate = ();
                foreach my $r (sort  keys %{$fastqs{$s}->{$d}->{$l}}){
                    push @reads, $fastqs{$s}->{$d}->{$l}->{$r};
                }
                if (scalar @reads != 2){
                    warn "WARNING: $s - $d - $l has " . scalar(@reads) . " reads\n";
                }
                my $bam = "$dir/$s-$d-$l.bam";
                push @{$fastqs{$s}->{bams}}, $bam;
                my $script = "$script_dir/align-$s-$d-$l.sh";
                $script = check_exists_and_rename($script);
                push @{$scripts{align}}, $script;
                if (-e $bam){
                    warn "WARNING: $bam BAM file already exists!\n";
                }
                open (my $SCRIPT, ">$script") or die "Can't open $script for writing: $!\n";
                my $rg = "'\@RG\\tID:$s-$d-$l\\tLB:$s-$d\\tSM:$s\\tPL:ILLUMINA'";
                my $readstring = "\"$reads[0]\"";
                foreach my $read (@reads){
                    my $rf = fileparse($read);
                    my $t = "$trim_dir/$rf";
                    my $p;
                    if ($rf =~ /L\d{1,3}_R(\d)/){
                        $p = $1;
                    }else{
                        die "Error parsing fastq '$read' - unable to determine read number for trim output.\n";
                    }
                    $t =~ s/\.fastq(\.gz)*$/_val_$p.fq.gz/ ;
                    push @trimmed, $t;
                }
                my $trimmedstring = "\"$trimmed[0]\"";
                my $trim_galore = "trim_galore";
                if (@reads > 1){
                    $readstring .= " \"$reads[1]\"";
                    $trimmedstring .= " \"$trimmed[1]\"";
                    $trim_galore .= " --paired";
                }
                if ($email){
                    print $SCRIPT <<EOT
#\$ -M $email
#\$ -m b
#\$ -m e
EOT
;
                }
                print $SCRIPT <<EOT

#\$ -cwd
#\$ -V
#\$ -e $script.stderr
#\$ -o $script.stdout
#\$ -pe sharedmem 8
#\$ -l h_rt=$runtime:00:00
#\$ -l h_vmem=$vmem
. /etc/profile.d/modules.sh
module load igmm/apps/bwa/0.7.12-r1039
module load igmm/apps/samtools/1.2
module load igmm/apps/R/3.2.2
module load igmm/apps/jdk/1.8.0_66
module load igmm/apps/TrimGalore/0.4.1

$trim_galore $readstring -o $trim_dir --path_to_cutadapt /exports/igmm/software/pkg/el7/apps/python/2.7.10/bin/cutadapt

bwa mem $fasta  -t $threads -M -R $rg $trimmedstring | samtools view -Sb - > "$bam"

EOT
;
            (my $out_bam = $bam) =~ s/\.bam$//;
            push @intermediate, $bam;#remove unsorted bam, keep sorted bam
            $out_bam .= "_sorted.bam";
            print $SCRIPT <<EOT
java -Djava.io.tmpdir=$tmp_dir -Xmx$mem -jar $picard SortSam SO=coordinate TMP_DIR=$tmp_dir CREATE_INDEX=TRUE I="$bam" O="$out_bam"

EOT
;

    #DEDUP
            my $in_bam = $out_bam;
            $out_bam =~ s/\.bam$//;
            my $metrics_file = "$out_bam" ."_rmdups.metrics";
            $out_bam .= "_rmdups.bam";
            push @intermediate, $out_bam;#remove dedupped
            print $SCRIPT <<EOT
java -Djava.io.tmpdir=$tmp_dir -Xmx$mem -jar $picard MarkDuplicates I="$in_bam" O="$out_bam" M="$metrics_file" CREATE_INDEX=TRUE TMP_DIR=$tmp_dir 

EOT
;


    #INDELREALIGN
            $in_bam = $out_bam;
            (my $stub = $in_bam) =~ s/\.bam$//;
            my $intervals = $stub . "_indelrealign.intervals";
            $out_bam = $stub . "_indelrealign.bam";
            push @intermediate, $out_bam;#remove realigned
            print $SCRIPT <<EOT
$call_gatk -T RealignerTargetCreator -known $indels -known $mills -I "$in_bam" -o "$intervals" $interval_string -nt $threads 


$call_gatk -T IndelRealigner -known $indels -known $mills -I "$in_bam" -targetIntervals "$intervals" -o "$out_bam" $interval_string 

EOT
;
    #BQSR           
            $in_bam = $out_bam;
            ($stub = $in_bam) =~ s/\.bam$//;
            my $grp = $stub . "_bqsr.grp";
            my $postrecal = $stub . "_bqsr_postrecal.grp";
            my $plots = $stub . "_bqsr_postrecal.plots.pdf";
            $out_bam = $stub . "_bqsr.bam";#keep this final bam
            my $gdir = "$outdir/gvcfs";
            if (not -d $gdir){
                mkdir $gdir or die "Could not make vcfs directory $gdir: $!\n";
            }
            my $f = fileparse($out_bam);
            $f =~ s/\.bam$//;
            my $gvcf = "$gdir/var.$f.g.vcf.gz";
            print $SCRIPT <<EOT
#get recalibration model
$call_gatk -T BaseRecalibrator -knownSites $dbsnp -knownSites $indels -knownSites $mills -I "$in_bam" -o "$grp" $interval_string -nct $threads 

#apply recalibration model
$call_gatk -T PrintReads -I "$in_bam" -o "$out_bam" -BQSR "$grp" $interval_string -nct $threads

#check recalibration model
$call_gatk -T BaseRecalibrator -knownSites $dbsnp -knownSites $indels -knownSites $mills -I "$in_bam" -BQSR "$grp" -o "$postrecal" $interval_string  -nct $threads 

#produce plots
$call_gatk -T AnalyzeCovariates -before "$grp" -after "$postrecal" -plots "$plots" $interval_string 

#make gvcf
$call_gatk -T HaplotypeCaller -I "$out_bam" --emitRefConfidence GVCF  -variant_index_type LINEAR -variant_index_parameter 128000 -o $gvcf $interval_string 

EOT
;           
            #remove intermediate bams
            @intermediate = map { "'$_'" } @intermediate;
            my $to_remove = join(" ", @intermediate); 
            print $SCRIPT <<EOT
#remove intermediate bams
rm $to_remove

EOT
;
                push @bams, $out_bam;
                push @gvcfs, $gvcf;
            }   
        }  
    }
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
        push @{$scripts{post_align}}, $gscript;
        open (my $SCRIPT, ">$gscript") or die "Can't open $gscript for writing: $!\n";
        my $depth = "$dir/depth.$stub";
        if ($email){
            print $SCRIPT <<EOT
#\$ -M $email
#\$ -m b
#\$ -m e
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


#################################################
sub make_genotype_script{
    my $gscript = "$script_dir/genotype-all-$date.sh";
    $gscript = check_exists_and_rename($gscript);
    my $dir = "$outdir/vcf";
    if (not -d $dir){
        mkdir $dir or die "Could not make vcfs directory $dir: $!\n";
    }

    my $samplesummarizer = "perl $RealBin/sampleSummarizer.pl";
    if (-e "$RealBin/sampleSummarizerEddieBinary"){
       $samplesummarizer = "$RealBin/sampleSummarizerEddieBinary" ;
    }

    push @{$scripts{post_align}}, $gscript;
    open (my $SCRIPT, ">$gscript") or die "Can't open $gscript for writing: $!\n";
    my $gs = "-V " . join(" -V ", @gvcfs);
    if ($email){
        print $SCRIPT <<EOT
#\$ -M $email
#\$ -m b
#\$ -m e
EOT
;
    }
    print $SCRIPT <<EOT
#\$ -cwd
#\$ -V
#\$ -e $gscript.stderr
#\$ -o $gscript.stdout
#\$ -pe sharedmem 8
#\$ -l h_rt=$runtime:00:00
#\$ -l h_vmem=$vmem
. /etc/profile.d/modules.sh
module load igmm/apps/R/3.2.2
module load igmm/apps/jdk/1.8.0_66
module load igmm/libs/htslib/1.3
module load igmm/apps/samtools/1.2

$call_gatk -T GenotypeGVCFs -stand_call_conf 30 -stand_emit_conf 4  -D $dbsnp $interval_string $gs -o $dir/var.$vcf_name-$date.raw.vcf.gz $interval_string

$call_gatk -T SelectVariants -selectType SNP -selectType MNP -V $dir/var.$vcf_name-$date.raw.vcf.gz -o $dir/var.snvs.$vcf_name-$date.raw.vcf.gz

$call_gatk -T VariantFiltration -V $dir/var.snvs.$vcf_name-$date.raw.vcf.gz -o $dir/var.snvs.$vcf_name-$date.filters.vcf.gz --filterExpression "QD < 2.0 " --filterName "QDfilterSNV" --filterExpression "FS > 60.0" --filterName "FSfilterSNV" --filterExpression "MQ < 40.0" --filterName "MQfilterSNV" --filterExpression "MQRankSum < -12.5" --filterName "MQRankSumFilterSNV" --filterExpression "ReadPosRankSum < -8.0" --filterName "ReadPosRankSumFilterSNV"
 
$call_gatk -T SelectVariants -selectType INDEL -selectType MIXED -V $dir/var.$vcf_name-$date.raw.vcf.gz -o $dir/var.indels.$vcf_name-$date.raw.vcf.gz

$call_gatk -T VariantFiltration -V $dir/var.indels.$vcf_name-$date.raw.vcf.gz -o $dir/var.indels.$vcf_name-$date.filters.vcf.gz --filterExpression "QD < 2.0" --filterName "QDindelFilter" --filterExpression "FS > 200.0" --filterName "FSindelFilter" --filterExpression "ReadPosRankSum < -20.0"  --filterName "ReadPosRankSumIndelFilter"

$call_gatk -T CombineVariants --genotypemergeoption UNSORTED -V $dir/var.snvs.$vcf_name-$date.filters.vcf.gz -V $dir/var.indels.$vcf_name-$date.filters.vcf.gz -o $dir/var.$vcf_name-$date.filters.vcf.gz

perl $vep_dir/variant_effect_predictor.pl  --dir $vep_dir/vep_cache --vcf --offline --everything --plugin LoF   --plugin MaxEntScan,$maxentscan --buffer_size 200 --check_alleles --gencode_basic --assembly GRCh37 -i $dir/var.$vcf_name-$date.filters.vcf.gz -o $dir/vep.var.$vcf_name-$date.filters.vcf --force --fork 8

bgzip -f $dir/vep.var.$vcf_name-$date.filters.vcf


EOT
;
    close $SCRIPT;

    my $sscript = "$script_dir/summarizer-$date.sh";
    $sscript = check_exists_and_rename($sscript);
    push @{$scripts{post_depth}}, $sscript;
    open (my $SUBSCRIPT, ">$sscript") or die "Can't open $sscript for writing: $!\n";
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

$samplesummarizer -i $dir/vep.var.$vcf_name-$date.filters.vcf.gz  -t $RealBin/genes/$gene_db   -n $not_reportable_cov -r $reportable_cov -s $dbsnp -e $evs -x $exac -z $cadd  -q $outdir/fastqc -c $outdir/depth -f $freq -o $outdir/sample_summaries/ -u $outdir/sample_summaries/summary.xlsx -l $gene_list

EOT
;
    close $SUBSCRIPT;
}

##################################################################################
sub submit_scripts{
    my @wait_ids = ();
    foreach my $sc (@{$scripts{fastqc}}){
        my $cmd = "qsub $sc";
        print STDERR "Executing: $cmd\n";
        my $output = `$cmd`;
        check_exit($?);
    }
    foreach my $sc (@{$scripts{align}}){
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
    @wait_ids = ();
    foreach my $sc (@{$scripts{post_align}}){
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
    $hold = join(",", @wait_ids);
    @wait_ids = ();
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
sub check_sample_name{
	my ($sample, $date, $lane, $read, $d, $f) = @_;
    if (not exists $fastqs{$sample}){
        return $sample;
    }
    my $rename = rename_sample($sample); 
    if ($fastqs{$sample}->{$date}){
        if ($fastqs{$sample}->{$date}->{$lane}){
            if ($fastqs{$sample}->{$date}->{$lane}->{$read}){
                warn "WARN: File $d/$f is a duplicate lane ($lane), read ($read) and sample ($sample). Trying sample name of $rename...\n";
                return check_sample_name($rename, $date, $lane, $read, $d, $f);
            }else{
                return $sample;#other read pair for sample - keep sample name
            }
        }
    }#sample exists and isn't a pair of an existing read - rename to something unique
    #warn "WARN: File $d/$f is a duplicate sample ($sample). Trying sample name of $rename...\n";
    #return check_sample_name($rename, $date, $lane, $read, $d, $f);
    return $sample;
}
##################################################################################
sub rename_sample{
    my $sample = shift;
	if ($sample =~ /\.(\d+)$/){
		my $dup = $1 + 1;
		$sample =~ s/\.(\d+)$//;
		$sample .= ".$dup";
	}else{
		$sample = "$sample.1";
	}
	return $sample;
}
##################################################################################
sub check_duplicates{
    my %bam_check = ();
    foreach my $bam (@bams){
        if (exists $bam_check{$bam}){
            die "Duplicate output bam: $bam!\n";
        }
    }
}
##################################################################################
sub usage{
    my $msg = shift; 
    print "ERROR: $msg\n" if $msg;
    print <<EOT
    
    USAGE: perl $0 /data/fastqs/*fastq.gz [options]
    
    From given FASTQs this script will create a series of qsub scripts to
    perform alignment and post-processing. 

    FASTQs must be provided as read-pairs and must be named in the following
    format: 
            
        [samplename]_S[0-9]_L[lanenumber]_R[readnumber]_[numbers].fastq(.gz)
    or
        [samplename]_[tenLetterBarcode]_L[lanenumber]_R[readnumber]_[numbers].fastq(.gz)
    
    e.g. 
        Sample1_S1_L001_R1_001.fastq
        Sample1_S1_L001_R2_001.fastq
        Sample2_ACGCGGACTA_L006_R1_002.fastq.gz
        Sample2_ACGCGGACTA_L006_R2_002.fastq.gz


    OPTIONS: 

    -o,--outdir DIR
        Optional location for output files. Default will be a new subdirectory 
        called alignments_and_vcfs_$date

    -s,--script_dir DIR
        Directory to write sub scripts to. Default = './sub_scripts'.
    
    -q,--qsub   
        Submit created scripts once finished - this automatically enables 
        queueing of the post-alignment commands for convenience.

    -f,--freq FLOAT
        Allele frequency cutoff for reporting variants in final summary document

    -d,--depth_intervals FILE
        Bed file for DepthOfCoverage. Default =
        $RealBin/bed_files/capture_regions_GRCh37.bed

    -n,--print_scripts
        Print names of created scripts once finished (if not using --qsub)

    -g,--gatk FILE
        Location of GATK jar file. Default =
        $RealBin/exe_bundle/GATK/v3.5/GenomeAnalysisTK.jar

    -p,--picard FILE
        location of picard jar file. Default =
        $RealBin/exe_bundle/picard/picard-tools-2.1.1/picard.jar

    -l,--list FILE(s)
        Interval list(s) to use with GATK commands. Default =
        $RealBin/bed_files/reportable_genes_GRCh37.bed and
        $RealBin/bed_files/not_reportable_genes_GRCh37.bed
    
    -e,--email STRING
        address to email script messages to

    -m,--mem INT
        amount of memory (in GB) for each GATK/Picard command. Default = 4

    -v,--vmem INT
        amount of vmem for each script. Default = 8

    -t,--threads INT
        no. threads for each command. Default = 8

    -r,--runtime INT
        runtime in hours for each script. Default = 4
 
    -x,--tmp_dir DIR
        directory to use for temporary storage for commands. Default =
        $ENV{HOME}/scratch/tmp/
    
    -i,--indels FILE
        1000 genomes phase1 indels. Default =
        $RealBin/ref_bundle/1000G_phase1.indels.b37.vcf

    -c,--cadd_dir DIR
        Directory containing pre-CADD-scored variants for scoring of output.
        Default = $RealBin/ref_bundle/cadd/v1.3/

    --dbsnp FILE
         dbSNP file. Default = $RealBin/ref_bundle/dbSNP146/All_20151104.vcf.gz

    --evs FILE
        EVS VCF file (all chromosomes in one file). Default =
        $RealBin/ref_bundle/ESP6500SI-V2-SSA137.GRCh38-liftover.combined.vcf.gz

    --exac  FILE         
        ExAC VCF file. Default =
        $RealBin/ref_bundle/ExAC.r0.3.sites.vep.vcf.gz]

    --fasta FILE     
        Genome reference fasta sequence. Default =
        $RealBin/bed_files/capture_regions_GRCh37.bed
    
    --gene_list FILE
        TSV file of gene symbol, expected inheritance pattern and associated
        conditions for annotating sampleSummarizer output. Default =
        $RealBin/genes/gene_inheritance_and_diseases.txt]

    --gene_database FILE     
        SQLITE database created using dbCreator.pl for all genes. This is used
        for annotating sampleSummarizer output. Default =
        $RealBin/genes/gene_database.db

    --mills FILE
        Mills and 1000 genomes indels for indelrealignment. Default =
        $RealBin/ref_bundle/Mills_and_1000G_gold_standard.indels.b37.vcf

    --reportable_cov FILE
        Bed file of regions in reportable genes to analyze coverage in (e.g.
        coding exons). Default =
        "$RealBin/bed_files/reportable_coding_exons_GRCh37.bed

    --not_reportable_cov FILE
        Bed file of regions in non-reportable genes to analyze coverage in -
        default = "$RealBin/bed_files/not_reportable_coding_exons_GRCh37.bed

    --vep_dir DIR       
        directory containing variant_effect_predictor.pl script and offline
        database named 'vep_cache'. Default =
        "$RealBin/exe_bundle/variant_effect_predictor

    --maxentscan DIR
        Directory containing maxentscan scoring programs. Default =
        "$RealBin/exe_bundle/maxentscan/fordownload/

    --vcf_name STRING
        Name for VCF outputs. Outputs will be named var.[vcf_name].raw.vcf etc.
        Default = 'panel'

    -h,--help
        Show this help message and exit


    AUTHOR:

    David A. Parry
    University of Edinburgh

EOT
;
    exit 1 if $msg;
    exit;
}
