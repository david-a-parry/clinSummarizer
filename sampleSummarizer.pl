#!/usr/bin/env perl 
#David A. Parry, April 2016

use strict;
use warnings;
use POSIX;
use DBI;
use Getopt::Long;
use Data::Dumper;
use Term::ProgressBar;
use POSIX qw/strftime/;
use List::Util qw(first sum max);
use Pod::Usage;
use File::Basename;
use FindBin;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use HTTP::Tiny;
use JSON;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Bio::DB::HTS::Tabix; 
use lib "$FindBin::Bin/dapPerlGenomicLib";
use VcfReader;
use SortGenomicCoordinates;

my %evs_acs = 
( 
    EA_AC => "EuropeanAmerican",
    AA_AC => "AfricanAmerican",
    TAC   => "Total",
) ;
my @exac_pops = qw / AFR AMR EAS FIN NFE SAS / ; 
my $progressbar;
my @allele_balance = ();
my @coverage_dirs  = ();
my @samples_to_analyze = (); 
my %opts = 
(
    b => \@allele_balance,
    f => 0,
    c => \@coverage_dirs,
    samples => \@samples_to_analyze,
    cadd_cutoff => 0,
);
GetOptions(
    \%opts,
    'a|AF=f',                     #AF cutoff within VCF
    'b|allele_balance=f{,}',      #min and optional max alt allele ratio per sample call
    'c|coverage=s{,}',            #directories with depth of GATK coverage data for each sample
    'cadd_cutoff=f',
    'd|depth=i',                  #optional min depth for sample call
    'e|evs=s',                    #evs VCF (concatanated)
    'f|allele_frequency=f',       #allele frequency cutoff for dbsnp, evs, exac
    'g|gq=f',                     #min GQ quality for calls
    'h|?|help',
    'i|input=s',                  #vcf input
    'l|gene_list=s',              #tsv of gene names, expected inheritence and condition
    'manual',
    'n|non_reportable_regions=s', #bed of regions in non-reportable genes for calculating coverage
    'o|output_prefix=s',          #output will be named prefixsamplename.xlsx
    'pathogenic',                 #only output variants if they're disease causing in clinvar/hgmd or LoF or affect CDD feature residue
    'p|progress',                 #show a progress bar?
    'q|fastqc_dir=s',             #directory with samples fastqc results
    'r|reportable_regions=s',     #bed of regions in reportable genes for calculating coverage
    'samples=s{,}',
    's|snp=s',                    #dbSNP VCF
    't|transcript_database=s',    #sqlite database of transcripts and protein info 
    'u|summary=s',                #summary file giving likely causative variants for samples
    'x|exac=s',                   #exac vcf
    'z|cadd_dir=s',               #directory containing tabix indexed CADD scores
) or pod2usage(-exitval => 2, -message => "Syntax error.\n"); 

pod2usage( -verbose => 1 ) if $opts{h};
pod2usage( -verbose => 2 ) if $opts{manual};
pod2usage( -exitval => 2, -message => "-i/--input is required" ) if (not $opts{i});
pod2usage( -exitval => 2, -message => "-t/--transcript_database is required" ) if (not $opts{t});

my $min_gq = defined $opts{g} ? 0 : $opts{g}; #default min GQ of 0
my %functional_classes = map {$_ => undef} 
qw /
    TFBS_ablation
    TFBS_amplification
    frameshift_variant
    inframe_deletion
    inframe_insertion
    initiator_codon_variant
    missense_variant
    protein_altering_variant
    regulatory_region_ablation
    regulatory_region_amplification
    splice_acceptor_variant
    splice_donor_variant
    stop_gained
    stop_lost
    transcript_ablation
    transcript_amplification
/;
my %lof_classes = map {$_ => undef} 
qw /
    frameshift_variant
    initiator_codon_variant
    splice_acceptor_variant
    splice_donor_variant
    stop_gained
    stop_lost
    transcript_ablation
/;


#open VCF, get samples and get VEP annotations
informUser("Checking input VCF\n");
my ($header, $first_var, $VCF)  = #this method means we can read from STDIN/pipe
 VcfReader::getHeaderAndFirstVariant($opts{i});
die "VCF header not OK for $opts{i}\n" 
 if not VcfReader::checkHeader(header => $header);

my %samples_to_columns = VcfReader::getSamples
(
    header => $header, 
    get_columns => 1
);
my %vep_fields = VcfReader::readVepHeader
(
    header => $header
);
my %info_fields = VcfReader::getInfoFields
(
    header => $header
);

if (@samples_to_analyze){
    foreach my $s (@samples_to_analyze){
        if (not exists $samples_to_columns{$s}){
            die "ERROR: User-specified sample '$s' does not exist in VCF!\n";
        }
    }
}else{
    @samples_to_analyze = sort 
    {
        $samples_to_columns{$a} <=> $samples_to_columns{$b} 
    } keys %samples_to_columns;
}

#store variants per sample in this hash for writing sample sheet
my %sample_vars = (); 

#setup our hash of headers
my %headers = getHeaders();

#set consequence ranks;
my %so_ranks = ();
setConsequenceRanks();

#check gene inheritance list
my %gene_conditions = ();
if ($opts{l}){
    checkGeneInheritanceFile();
}elsif($opts{u}){
    informUser
    (
        "ERROR: -u/--summary option requires a gene list with expected ".
        "inheritance and conditions for each gene in order to produce a ".
        "summary for likely pathogenic variants for each sample.\n"
    );
    $opts{u} = undef;
}

#check coverage directory contains required files if specified
my %sample_coverage = (); #hash of sample names to summary or depth files
my %regions = (); #hash to tabix iterators reportable/non-reportable regions
if (@coverage_dirs){
    checkCoverageDirs();
    #readCoverageBeds();
    #getBedIters();
}

#find CADD score files and get tabix iterators
my @cadd_iters = ();
if ($opts{z}){#cadd dir
    getCaddIters();
}

#get fastqc results files for each sample
my %fastqc = ();#hash of sample name to array of fastqc results files
if ($opts{q}){
    getFastqcSummaries();
}
#open and check transcript database
informUser("Checking transcript database\n");
my $dbh;
my $driver   = "SQLite";
my %search_handles = ();
my %transcript_ranks = ();
my %enst_xref = ();
readTranscriptDatabase();

#check dbSNP, evs and exac files
my %snp_search_args  = readDbSnpFile();
my %evs_search_args  = readEvsFile();
my %exac_search_args = readExacFile();

#create output directory if required
if ($opts{o}){
    my ($f, $d) = fileparse($opts{o});
    if ($d and not -d $d){
        mkdir($d) or die "Could not create output directory '$d' for output '$opts{o}: $!\n";
    }
}
my $summ_xl;#xlsx file for summaries of possible diagnoses    
my %summ_worksheets = (); 
my %genes_per_sample; #keep variants to check if compatible with inheritance
if ($opts{u}){
    setupSummaryFile(); 
}

    

#set up progress bar if user wants one
my $total_var;
my $next_update = 0;
my $n = 0;
if ($opts{p}){
    informUser("Calculating file length for progress bar...\n");
    $total_var = VcfReader::countVariants( $opts{i} );
    informUser("$opts{i} has $total_var variants.\n");
    $progressbar = Term::ProgressBar->new(
        { name => "Analyzing", 
          count => ($total_var), 
          ETA => "linear" 
        } 
    );
}

#get filehandle and start reading variants
informUser("Commencing variant analysis\n");
assessVariant($first_var);
$n++;
while (my $l = <$VCF>){
    next if $l =~ /^#/;
    assessVariant($l);
    $n++;
    if ($progressbar){
        $next_update = $progressbar->update($n) if $n >= $next_update;
    }
}
close $VCF;
if ($progressbar){
    $progressbar->update($total_var) if $total_var >= $next_update;
}

$progressbar = undef;

my %depth = ();
#$depth{bed_line}->{mean}->{sample_id} = mean coverage of bed region
#$depth{bed_line}->{pc20}->{sample_id} = percent bases covered 20X
#$depth{bed_line}->{mean|pc20}->{mean} = mean values for samples 
#$depth{bed_line}->{mean|pc20}->{median} = median values for samples 
getSampleDepths();#get depths first and calculate mean, median

writeSampleSummary();

writeGeneSummary();

informUser("Done.\n");

###########################################################
sub writeGeneSummary{
    return if not $summ_xl;
    my $row = 1;
    my $sample_row = 0;
    foreach my $s (@samples_to_analyze){
        my $gene_col = 0; 
        foreach my $gene (keys %{$genes_per_sample{$s}}){
            my $do_write = 0; 
            if ($gene_conditions{$gene}->{inheritance} =~ /^[AX]R$/){
                my $alleles = 0; 
                my %phased = (); 
                for (my $i = 0; $i < @{$genes_per_sample{$s}->{$gene}}; $i++){
                    my $h = $genes_per_sample{$s}->{$gene}->[$i];
                    if ($h->{count} == 2){#homozygous
                        $alleles = 2;
                        last;
                    }
                    my $pgt = getFormatField
                    (
                        var     => $h->{var}, 
                        field   => 'PGT',
                        sample  => $s,
                    ); 
                    if (not $pgt){
                        $alleles++;#if not phased,assume in trans with any other
                        last if $alleles > 1;
                        next;
                    }
                    my @pgt_al = split(/\|/, $pgt);
                    my $p1 = first { $pgt_al[$_] eq $h->{allele} } 0 .. $#pgt_al;
                    if (not defined $p1){#this allele not phased
                        $alleles++;
                        last if $alleles > 1;
                        next;
                    }
                    my $pid = getFormatField
                    (
                        var     => $h->{var}, 
                        field   => 'PID',
                        sample  => $s,
                    ); 
                    #get phased alleles and only count those in cis as one
                    $phased{"$pid:$p1"} = undef; 
                    
                }
                $alleles += keys %phased;#add phased alleles, 
                                         #counting those in cis only once
                if ($alleles > 1 ){#recessive and at least two alleles - write
                    $do_write = 1;
                }
            }else{
                $do_write = 1;#not recessive inheritence, write regardless
            }
            if ($do_write){
                my $col = 0; 
                foreach my $var (@{$genes_per_sample{$s}->{$gene}}){
                    my $col = 0;
                    foreach my $v (@{$var->{columns}}){
                        $summ_worksheets{Variants}->write($row, $col++, $v);
                    }
                    $row++;
                }
                if ($gene_col == 0){#haven't written sample name yet
                    $summ_worksheets{Samples}->write(++$sample_row, $gene_col++, $s);
                }
                $summ_worksheets{Samples}->write($sample_row, $gene_col++, $gene);
            }
        }
    }
    $summ_xl->close;
}

###########################################################
sub getFormatField{
    my %args = @_;
    my %format = VcfReader::getVariantFormatFields($args{var});
    return undef if not exists $format{$args{field}}; 
    my $value = VcfReader::getSampleGenotypeField
    (
        line              => $args{var},
        field             => $args{field},
        sample            => $args{sample}, 
        sample_to_columns => \%samples_to_columns,
    ); 
}
###########################################################
sub writeFastqcSummaries{
    my $sample = shift;
    my $ws = shift;
    my $hf = shift;
    if (exists $fastqc{$sample}){
        my %results = ();
        my %cats = ();
        informUser("Writing fastqc summary for $sample\n");
        foreach my $f (@{$fastqc{$sample}} ){
            readFastqcResults($f, \%results, \%cats);
        }
        my $col = 0;
        foreach my $h ("Catergory", sort keys %results ){
            $ws->write(0, $col++, $h, $hf);
        }
        my $row = 1;
        foreach my $c (sort keys %cats){
            $ws->write($row++, 0, $c);
        }
            
        $col = 1;#column of first fastq file
        foreach my $fq (sort keys %results){
            $row = 1;
            foreach my $c (sort keys %cats){
                if (exists $results{$fq}->{$c}){
                    $ws->write($row++, $col, $results{$fq}->{$c});
                }
            }
            $col++;
        }
    }
}

###########################################################
sub readFastqcResults{
    my $f = shift;
    my $res = shift;
    my $cat = shift;
    open (my $RES, $f) or informUser
    (
        "WARNING: Error opening FASTQC results file $f: $!\n"
    ) and return;
    my %stats = ();
    while (my $line  = <$RES>){
        chomp $line;
        next if not $line;
        my @split = split("\t", $line);
        my $fq = $split[2];
        my $stat = $split[1];
        $res->{$fq}->{$stat} =  $split[0];
        $cat->{$stat} = undef;
    }
}

###########################################################
sub getSampleDepths{
    return if not @coverage_dirs;
    return if not $opts{r} and not $opts{n};
    foreach my $s (@samples_to_analyze){
        if (exists $sample_coverage{$s}->{depth}){
            if ($opts{r}){
                getDepthPerRegion
                (
                    $sample_coverage{$s}->{depth}, 
                    $s,
                    $opts{r}, 
                );
            }
            if ($opts{n}){
                getDepthPerRegion
                (
                    $sample_coverage{$s}->{depth}, 
                    $s,
                    $opts{n}, 
                );
            }
        }
    }
    foreach my $region (keys %depth){
        foreach my $k (keys %{$depth{$region}}){
            my @values = ();
            foreach my $samp (keys  %{$depth{$region}->{$k}}){
                push @values, $depth{$region}->{$k}->{$samp};
            }
            my $mean = mean(@values);
            $depth{$region}->{$k}->{mean} = $mean;
            my $median = median(@values);
            $depth{$region}->{$k}->{median} = $median;
        }
    }
}

###########################################################
sub mean{
    my $total = sum(@_) || 0;
    return scalar @_ != 0 ? $total / @_ : 0;
}

###########################################################
sub median{
    return sum( ( sort { $a <=> $b } @_ )[ int( $#_/2 ), ceil( $#_/2 ) ] )/2;
}

###########################################################
sub writeDepthSummaries{
    my $sample = shift;
    my $worksheets = shift; 
    my $header_format = shift;
    if (exists $sample_coverage{$sample}->{summary}){
        informUser("Writing depth summary for $sample\n");
        writeSummaryDepth
        (
            $sample_coverage{$sample}->{summary}, 
            $worksheets->{Depth_Summary},
            $header_format,
        
        );
    }
    if (exists $sample_coverage{$sample}->{depth}){
        informUser("Writing depth per region for $sample\n");
        writeDepthPerRegion
        (
            $worksheets,
            $sample,
            $header_format,
        );
    }
}

###########################################################
sub writeSummaryDepth{
    my $f  = shift;
    my $ws = shift;
    my $hf = shift;
    open (my $SUM, $f) or informUser
    (
        "Could not open coverage summary file '$f': $!\n" 
    )  and return;
    my $row = 0; 
    while (my $line = <$SUM>){
        chomp $line;
        my @split = split("\t", $line);
        my $col = 0;
        foreach my $s (@split){
            if ($row == 0 or $col == 0){
                $ws->write($row, $col++, $s, $hf);
            }else{
                $ws->write($row, $col++, $s);
            }
        }
        $row++;
    }
    close $SUM;
}

###########################################################
sub writeDepthPerRegion{
    my ($worksheets, $sample, $hf) = @_;
    
    if ($opts{r}){
        writeDepthOutput
        (
            $sample,
            $opts{r}, 
            $worksheets->{Depth_ReportableExons},
            $hf,
        );
    }
    if ($opts{n}){
        writeDepthOutput
        (
            $sample,
            $opts{n}, 
            $worksheets->{Depth_NonReportableExons},
            $hf,
        );
    }
}

###########################################################
sub writeDepthOutput{
    my ($sample, $bed, $ws, $hf) = @_;
    open (my $BED, $bed) or informUser
    (
        "WARNING: Could not open $bed for reading: $!\n"
    ) and return;
    my $col = 0;
    foreach my $h (
    qw /
        Chrom 
        Start 
        End 
        MeanDepth 
        RunMeanDepth 
        RunMedianOfTheMeanDepth 
        %Covered_20X
        RunMean%Covered_20X
        RunMedian%Covered_20X
        OtherBedFields
    / ){
        $ws->write(0, $col++, $h, $hf);
    }
    my $row = 0; 
    while (my $line = <$BED>){
        $col = 0;
        $row++;
        chomp $line;
        my @bed = split("\t", $line); 
        my $length = $bed[2] - $bed[1];
        my @other = ();
        if (@bed > 3){
            @other = @bed[3..$#bed];
        }
        next if $length < 1;
        foreach my $v 
        (
            @bed[0..2], 
            $depth{$line}->{mean}->{$sample}, 
            $depth{$line}->{mean}->{mean}, 
            $depth{$line}->{mean}->{median}, 
            $depth{$line}->{pc20}->{$sample}, 
            $depth{$line}->{pc20}->{mean}, 
            $depth{$line}->{pc20}->{median}, 
            @other,
        ){
            $ws->write($row, $col++, $v, );
        }
    }
}

###########################################################
sub getDepthPerRegion{
    my $f = shift;
    my $s = shift;
    my $bed = shift;
    informUser("Calculating depth data for $s vs $bed\n");
    open (my $BED, $bed) or informUser
    (
        "WARNING: Could not open $bed for reading: $!\n"
    ) and return;
    open (my $DEPTH, $f) or informUser
    (
        "WARNING: Could not open coverage file '$f': $!\n" 
    )  and return;
    my @header = split("\t", <$DEPTH>);
    chomp @header;
    my $d_col = 0;
    $d_col++ until $header[$d_col] eq "Depth_for_$s" or $d_col > $#header;
    if ($d_col > $#header){
        informUser(
            "WARNING: Could not find sample depth column in $f for $s!\n"
        );
        return;
    }
    my $converted = "$f.ssconverted";
    open (my $DCONV, ">", $converted) or informUser
    (
        "WARNING: Could not write to converted coverage file '$converted': $!\n" 
    )  and return;
    print $DCONV join("\t", "#CHROM", "POS", @header[1..$#header]) . "\n";
    my %per_locus_depth = ();
    while (my $line = <$DEPTH>){
        chomp $line;
        my @split = split("\t", $line);
        my ($chrom, $pos) = split(":", $split[0]); 
        print $DCONV join("\t", $chrom, $pos, @split[1..$#split]) . "\n";
    }
    close $DEPTH;
    close $DCONV;
    my $ti = compressAndIndexBgzip($converted);
    return if not $ti;
    while (my $line = <$BED>){
        chomp $line;
        my @bed = split("\t", $line); 
        my $length = $bed[2] - $bed[1];
        my @other = ();
        next if $length < 1;
        my $iter = $ti->query("$bed[0]:" . ($bed[1] + 1) ."-$bed[2]");
        my @depths = ();
        while (my $hit =  $iter->next() ){ 
            push @depths, (split "\t", $hit)[$d_col];
        }
        my $mean = mean(@depths);
        my $bp_over_19 = grep { $_ > 19 } @depths;
        my $pc_covered_20 = ($bp_over_19/$length) * 100;
        $depth{$line}->{mean}->{$s} = $mean;
        $depth{$line}->{pc20}->{$s} = $pc_covered_20;
    }
    close $BED;
}

###########################################################
sub compressAndIndexBgzip{
    my $f = shift;
    my $er = `bgzip -f $f`; 
    if ($?){
        if ($? == -1){
            informUser
            (
                "ERROR: bgzip compression of depth file $f ".
                "failed - bgzip command failed to execute.\n"
            ) and return;
        }
        my $exit = $? >> 8;
            informUser
            (
                "ERROR: bgzip compression of depth file $f ".
                "failed code $exit: $er.\n"
            ) and return;
    }
    $er = `tabix -s 1 -b 2 -e 2 $f.gz`; 
    if ($?){
        if ($? == -1){
            informUser
            (
                "ERROR: tabix indexing of depth file $f.gz ".
                "failed - bgzip command failed to execute.\n"
            ) and return;
        }
        my $exit = $? >> 8;
            informUser
            (
                "ERROR: tabix indexing of depth file $f.gz ".
                "failed code $exit: $er.\n"
            ) and return;
    }
    return Bio::DB::HTS::Tabix->new(filename => "$f.gz");
}

###########################################################
sub writeSampleSummary{
    foreach my $s (@samples_to_analyze){
        informUser("Preparing output file for $s\n");
        my ($workbook, $worksheets, $h_format) = setupOutput($s);
        informUser("Writing variants for $s\n");
        foreach my $sheet (qw / Functional Other / ){
            my $row = 1;
            foreach my $var (@{$sample_vars{$s}->{$sheet}}){
                my $col = 0;
                foreach my $v (@$var){
                    $worksheets->{$sheet}->write($row, $col++, $v);
                }
                $row++;
            }
        }
        #read depth file per sample and output to Depth sheet if coverage dir supplied
        if (@coverage_dirs){
           writeDepthSummaries($s, $worksheets, $h_format);
        }
        #read fastqc summary files and output to Fastqc sheet if fastqc_dir supplied
        if (exists $opts{q}){
            writeFastqcSummaries($s, $worksheets->{Fastqc}, $h_format);
        }
        $workbook->close;
    }
}

###########################################################
sub setupSummaryFile{
    my $x = $opts{u} =~ /\.xlsx$/ ? $opts{u} : "$opts{u}.xlsx";
    $summ_xl = Excel::Writer::XLSX->new($x) or die "Error creating XLSX file $x: $!\n";
    $summ_worksheets{Samples} = $summ_xl->add_worksheet("Samples"); 
    $summ_worksheets{Variants} = $summ_xl->add_worksheet("Variants"); 
    my $header_formatting = $summ_xl->add_format(bold => 1);
    my $col = 0;
    foreach my $h (@{$headers{Variants}}){
        $summ_worksheets{Variants}->write(0, $col++, $h, $header_formatting);
    }
    $col = 0; 
    foreach my $h (qw /Sample Genes/){
        $summ_worksheets{Samples}->write(0, $col++, $h, $header_formatting);
    }
}

###########################################################
sub setupOutput{
    my @suffixes = (".vcf", ".vcf.gz", ".txt");
    my $s = shift;
    my $sample_out; 
    if (not $opts{o}){
        my ($out, $dir, $extension) = fileparse($opts{i}, @suffixes);
        $sample_out = "$dir/$out.$s.xlsx";
    }else{
        $sample_out = "$opts{o}$s.xlsx"; 
    }
    my $xl = Excel::Writer::XLSX->new($sample_out);
    my %worksheets = ();
    $worksheets{Functional} = $xl->add_worksheet("Functional"); 
    $worksheets{Other} = $xl->add_worksheet("Other"); 
    if (exists $sample_coverage{$s}->{summary}){
        $worksheets{Depth_Summary} = $xl->add_worksheet("Depth_Summary");
    }
    if (exists $sample_coverage{$s}->{depth}){
        if ($opts{r}){
            $worksheets{Depth_ReportableExons} = $xl->add_worksheet("Depth_ReportableExons");
        }
        if ($opts{n}){
            $worksheets{Depth_NonReportableExons} = $xl->add_worksheet("Depth_NonReportableExons");
        }
    }
    if (exists $fastqc{$s}){
        $worksheets{Fastqc} = $xl->add_worksheet("Fastqc");
    }
    my $header_formatting = $xl->add_format(bold => 1);
    foreach my $k (qw/ Functional Other /){
        my $col = 0;
        foreach my $h (@{$headers{Variants}}){
            $worksheets{$k}->write(0, $col++, $h, $header_formatting);
        }
    }
    return ($xl, \%worksheets, $header_formatting);
}



###########################################################
sub assessVariant{
    my $line = shift;
    my @split = split("\t", $line); 
    #simplify alleles
    my %min  = VcfReader::minimizeAlleles(\@split);

    #determine VEP transcript/consequence to report
    #get vep csqs
    my @vep_csq = getVepCsq(\@split);
    my $ref  = VcfReader::getVariantField(\@split, 'REF');
    my @alts = split(",", VcfReader::getVariantField(\@split, 'ALT'));
    my $qual = VcfReader::getVariantField(\@split, 'QUAL');
    my $filter = VcfReader::getVariantField(\@split, 'FILTER');
    my $an = VcfReader::getVariantInfoField(\@split, "AN");
    my @af = split(",",  VcfReader::getVariantInfoField(\@split, "AF") );
    my @ac = split(",",  VcfReader::getVariantInfoField(\@split, "AC") );
    my @vep_alts = VcfReader::altsToVepAllele
    (
        ref  => $ref,
        alts => \@alts,
    );
    my %alt_to_vep = ();
    @alt_to_vep{@alts} = @vep_alts;
    
    #process alleles one at a time
    foreach my $al (sort {$a<=>$b} keys %min){
        my @row = ();#array of values for output to spreadsheet row
        #skip missing alleles
        next if ($min{$al}->{ORIGINAL_ALT} eq '*');
        if ($opts{a}){#filter on AF field in VCF
            next if $af[$al-1] >= $opts{a};
        }
        #get hash of sample names to array of sample genotype fields
        my %sample_genos = getSamplesWithAllele(\@split, $min{$al});
        next if not keys %sample_genos;#no sample with valid genotype

        #get relevant consequence to report
        my $csq_to_report = getConsequenceToReport(\@vep_csq, \%alt_to_vep, $min{$al});
        next if not $csq_to_report;#i.e. no vep consequence overlaps a gene of interest
        my $most_damaging_csq = getMostDamagingConsequence($csq_to_report);
        
        #check dbSNP
        my ($dbsnp_freq, $dbsnp_id) = checkDbSnp($min{$al});
        next if $dbsnp_freq > $opts{f} and $opts{f};
        #check EVS
        my ($evs_freq, $evs_pop) = checkEvs($min{$al});
        next if $evs_freq > $opts{f} and $opts{f};
        #check ExAC
        my ($exac_freq, $exac_pop) = checkExac($min{$al});
        next if $exac_freq > $opts{f} and $opts{f};

        #add gene expected inheritance and condition
        foreach my $k (qw /inheritance condition reportable/){
            push @row, exists $gene_conditions{$csq_to_report->{symbol}}->{$k} ? $gene_conditions{$csq_to_report->{symbol}}->{$k} : "";
        }
        #output VEP annotations, CHROM, REF etc.
        foreach my $f (qw / symbol feature consequence hgvsc hgvsp exon intron / ){
            push @row, $csq_to_report->{$f};
        }
        push @row, $csq_to_report->{overlap_features};
        push @row, $csq_to_report->{cdd_feature_residues};
            
        foreach my $f (qw / lof lof_filter lof_info lof_flags polyphen sift maxentscan_diff maxentscan_alt maxentscan_ref/ ){
            push @row, $csq_to_report->{$f};
        }
        if (exists $search_handles{et} ){
            push @row, getEtScore($csq_to_report); 
        }else{
            push @row, '';
        }
        
        #get CADD scores
        my $cadd_score = getCaddScore($min{$al});
        push @row, $cadd_score;

        #get HGMD and ClinVar matches
        my ($hgmd, $hgmd_dm)  = addHgmdMatches($min{$al}, $csq_to_report, $most_damaging_csq);
        push @row, @$hgmd; 
        my ($clinvar, $clinvar_path) = addClinvarMatches($min{$al}, $csq_to_report, $most_damaging_csq);
        push @row, @$clinvar;
        
        #output frequency values
        push @row, $dbsnp_freq >= 0 ? $dbsnp_freq : "";
        push @row, $dbsnp_id ? $dbsnp_id : "";
        push @row, $evs_freq >= 0 ? $evs_freq : "";
        push @row, $evs_pop ? $evs_pop : "";
        push @row, $exac_freq >= 0 ? $exac_freq : "";
        push @row, $exac_pop ? $exac_pop : "";
            
        foreach my $f ( qw / CHROM ORIGINAL_POS ORIGINAL_REF ORIGINAL_ALT / ){
            push @row, $min{$al}->{$f};
        }
        #push @row, VcfReader::getVariantField(\@split, 'ID');
        push @row, $qual;
        push @row, $filter;
        push @row, $an;
        push @row, $af[$al-1];
        push @row, $ac[$al-1];
        #record these details against each sample with variant
        my $sheet = "Other";
        if (exists $functional_classes{$most_damaging_csq}){
            if ( $opts{pathogenic} ){
                if (( $clinvar_path or $hgmd_dm or 
                 exists $lof_classes{$most_damaging_csq} 
                 or $csq_to_report->{cdd_feature_residues} ) 
                ){
                    $sheet = "Functional";
                }
            }else{
                $sheet = "Functional";
            }
        }
        foreach my $s (keys %sample_genos){
            push @{$sample_vars{$s}->{$sheet}}, [@{$sample_genos{$s}}, @row];
            if ($opts{u} and $sheet eq "Functional"){#collect sample variants for inheritance analysis
                if ($cadd_score eq '.' or $cadd_score >= $opts{cadd_cutoff}){
                    my $h = 1;
                    my $alt = $min{$al}->{ORIGINAL_ALT};
                    $h = 2 if $sample_genos{$s}->[1]  =~  /^$alt[\/\|]$alt$/;
                    push @{$genes_per_sample{$s}->{$csq_to_report->{symbol}}}, 
                    {
                        count   => $h, 
                        columns => [@{$sample_genos{$s}}, @row],
                        var     => \@split,
                        allele  => $al,
                    };
                }
            }
        }
    }
}
###########################################################
sub getSamplesWithAllele{
    my $l = shift;
    my $var = shift;
    my %sample_vars = ();
    my %samp_gts = VcfReader::getSampleActualGenotypes
        (
              line => $l, 
              all => 1,
              sample_to_columns => \%samples_to_columns,
              minGQ => $min_gq,
        );
    my %samp_gqs = VcfReader::getSampleGenotypeField
        (
              line => $l, 
              field => "GQ",
              all => 1,
              sample_to_columns => \%samples_to_columns,
        );
    my %samp_pls = VcfReader::getSampleGenotypeField
        (
              line => $l, 
              field => "PL",
              all => 1,
              sample_to_columns => \%samples_to_columns,
        );

    my %samp_ads = (); 
    foreach my $s (@samples_to_analyze){
        my @ads = VcfReader::getSampleAlleleDepths 
         (
              line => $l,
              sample => $s,
              sample_to_columns => \%samples_to_columns,
        );
        $samp_ads{$s} = \@ads;
    }    
    foreach my $s (@samples_to_analyze){
        my @alts = split(/[\/\|]/, $samp_gts{$s});
        if (grep { $_ eq $var->{ORIGINAL_ALT} } @alts ){ 
            my @ads = @{$samp_ads{$s}}; 
            my $depth = sum(@ads);
            if ($opts{d}){
                next if $opts{d} > $depth;
            }
            if ($opts{g}){
                eval "$opts{g} <= $samp_gqs{$s}" or next;
            }
            if ($opts{pl}){
                my @pls = split(",", $samp_pls{$s});
                next if $pls[0] < $opts{pl};
            }
            my $ab = 0;
            if ( $depth > 0){
                $ab = $ads[$var->{ALT_INDEX}]/$depth;
            }
            if (@{$opts{b}}){
                next if $ab < $opts{b}->[0];
                if (scalar @{$opts{b}} > 1){
                    next if $ab > $opts{b}->[1];
                }
            }
            push @{$sample_vars{$s}},  $s, $samp_gts{$s}, join(",", @ads) , $ab, $samp_gqs{$s} ;
        }
    } 
    return %sample_vars; 
}
###########################################################
sub getVepCsq{
    my $l = shift;
    my @get_vep =  
      qw /
            Symbol
            Gene
            Consequence
            Allele
            feature
            canonical
            hgvsc
            hgvsp
            exon
            intron
            polyphen
            sift
            DOMAINS
            Amino_acids
            protein_position
            ensp
            LoF
            LoF_Filter
            LoF_info
            LoF_flags
            MaxEntScan_alt
            MaxEntScan_diff
            MaxEntScan_ref

      /;
    return VcfReader::getVepFields
    (
        line       => $l,
        vep_header => \%vep_fields,
        field      => \@get_vep,
    );
}

###########################################################
sub getConsequenceToReport{
    my $vep_csq = shift;
    my $alt_to_vep = shift;
    my $var = shift;
    my @csq_to_rank = ();
    foreach my $csq (@$vep_csq){ 
        #collect info for each transcript before selecting which to use in report
        #skip consequences for different alleles if variant site is multiallelic
        next if $csq->{allele} ne $alt_to_vep->{ $var->{ORIGINAL_ALT} };
        #skip if this gene symbol is not listed in our database
        if (keys %transcript_ranks){
            my $gene = uc($csq->{gene});
            next if not exists $transcript_ranks{$gene};
        }
        push @csq_to_rank , $csq;
    }
    return if not @csq_to_rank;#no matching gene symbols for allele
    return rankTranscriptsAndConsequences(\@csq_to_rank); 
}

###########################################################
sub rankConsequences{
    my @vars = @_;#array of hashes of VEP consequences
    return sort { 
        $so_ranks{getMostDamagingConsequence($a)} <=> 
        $so_ranks{getMostDamagingConsequence($b)} 
    } @vars;
    
}

###########################################################
sub getTranscriptsRanks{
    my $symbol = shift;
    my $transcript = shift; 
    return -1 if not exists $transcript_ranks{$symbol};
    return -1 if not exists $transcript_ranks{$symbol}->{$transcript};
    return $transcript_ranks{$symbol}->{$transcript};
}

###########################################################
sub rankTranscripts{
    my @vars = @_;#array of hashes of VEP consequences
    return sort { getTranscriptsRanks( $a->{gene}, $a->{feature} ) <=> getTranscriptsRanks( $b->{gene}, $b->{feature} ) } @vars;
}

###########################################################
sub rankTranscriptsAndConsequences{
    my $csq_array = shift;#ref to array of VEP consequences 
    @$csq_array = rankConsequences(@$csq_array); 
    my $most_damaging = getMostDamagingConsequence($csq_array->[0]) ;
    @$csq_array = rankTranscripts(@$csq_array); 
    if($opts{t}){
        #if not checking rules just use getCddAndUniprotOverlappingFeatures 
        #to annotate any overlapping features which will be added to $csq hash
        foreach my $csq (@$csq_array){
            getCddAndUniprotOverlappingFeatures($csq);
        }
    }
    return first { $_->{consequence} =~ $most_damaging } @$csq_array;
}

###########################################################
sub getMostDamagingConsequence{
    my $csq = shift;#hash ref to VEP consequences for single transcript/feature
    my @s_csq = split("&", $csq->{consequence} );
    @s_csq = sort { $so_ranks{$a} <=> $so_ranks{$b} } @s_csq;
    return $s_csq[0];
}
 
###########################################################
sub getCddAndUniprotOverlappingFeatures{
    #returns a ref to an array of hashes of overlapping features
    my @hits = (); 
    return \@hits if not $dbh;
    my $csq = shift;
    my $uniprot = $enst_xref{$csq->{feature}}->{uniprot}; 
    return \@hits if not $uniprot;
    return \@hits if not $csq->{ensp};#skip non-coding
    return \@hits if not $csq->{amino_acids};#skip any non-coding variant
    return \@hits if (split("/", $csq->{amino_acids}) < 2); 
    return \@hits if not $csq->{protein_position};
    my ($p_start, $p_end) = split("-", $csq->{protein_position}); 
    $p_end ||= $p_start;
    #get overlapping uniprot features
    $search_handles{uniprot}->execute($uniprot, $p_start, $p_end)
      or die "Error searching 'uniprot' table in '$opts{t}': " . 
      $search_handles{uniprot} -> errstr;
    while (my @row = $search_handles{uniprot}->fetchrow_array()) {
        my $hash = 
        { 
            type    => "uniprot", 
            start   => $row[2],
            end     => $row[3],
            feature => $row[4],
            note    => $row[5],
            grch37pos => $row[6],
            grch38pos => $row[7],
        }; 
        push @hits, $hash;
    }
    
    $search_handles{cdd}->execute($uniprot, $p_start, $p_end)
      or die "Error searching 'cdd' table in '$opts{t}': " . 
      $search_handles{cdd} -> errstr;
    while (my @row = $search_handles{cdd}->fetchrow_array()) {
        my $type = '';
        if ($row[2] eq 'Feature'){
            $type = 'cdd_feature';
        }elsif ($row[2] eq 'Hit'){
            $type = 'cdd_hit';
        }else{
            informUser("WARNING: Do not understand ResultType field '$row[2]' in cdd table of $opts{t} - ignoring.\n");
            next;
        }
        my $hash = 
        { 
            type      => $type, 
            feature   => $row[3],
            residues  => $row[4],#residues are only present in 'Feature' types, not 'Hit'
            start     => $row[5],
            end       => $row[6],
            grch37pos => $row[7],
            grch38pos => $row[8],
        }; 
        push @hits, $hash;
    }
    my @hit_summary = () ;
    my @res_summary = () ;
    foreach my $h (@hits){
        push @hit_summary, "$h->{type},$h->{feature},$h->{start}-$h->{end}";
        if ($h->{residues}){
            my @residues = split(/\,/, $h->{residues}); 
            foreach my $residue (@residues){
                if ($residue =~ /^([A-Z])(\d+)$/){
                    my ($res, $pos) = ($1, $2);
                    last if $pos > $p_end;
                    if ($pos >= $p_start and $pos <= $p_end){
                        push @res_summary, "$h->{feature},$residue";
                    }
                }
            }
        }
    }
    if (@hit_summary){
        $csq->{overlap_features} = join("/", @hit_summary); 
        $csq->{cdd_feature_residues} = join("/", @res_summary); 
    }
    return \@hits;
}


###########################################################
sub addHgmdMatches{
    my $var = shift;
    my $csq = shift;
    my $most_damaging_csq = shift;
    my @h_matches = getHgmdMatches($var);
    my $path = 0;
    my @results = ();
    if (@h_matches){
        $path++ if grep { $_->{variant_class} =~ /^D[MP]$/ } @h_matches;
        foreach my $f 
        ( qw /
                hgmd_id
                disease
                variant_class
                gene_symbol
                hgvs
            /
        ){#add HGMD results to row, multiple values per field separated by commas
            my $s = join(",", map { $_->{$f} } @h_matches );
            push @results, $s;
        }
    }else{
        push @results, map {"-"} (1..5); 
    }

    #get variants with same AA altered
    my @aa_matches = ();
    if (exists $search_handles{hgmd_aa} and 
        exists $search_handles{hgmd_id} and 
        ( $most_damaging_csq eq 'missense_variant' or 
          $most_damaging_csq eq 'protein_altering_variant' or 
          $most_damaging_csq =~  /^inframe_(inser|dele)tion$/ 
        )
    ){
        $search_handles{hgmd_aa}->execute
        (
            $csq->{feature},
            $csq->{protein_position},
        ) or die "Error searching 'HGMD_VEP' table in '$opts{t}': " . 
          $search_handles{hgmd_aa} -> errstr;
        while (my @db_row = $search_handles{hgmd_aa}->fetchrow_array()) {
            my ($hgmd_id, $feature, $consequence, $protein_position, $aa, $hgvsc, $hgvsp) = @db_row; 
            next if ($hgmd_id and grep {$_->{hgmd_id} eq $hgmd_id} @h_matches );
            my @s_csq = split("&", $consequence); 
            @s_csq = sort { $so_ranks{$a} <=> $so_ranks{$b} } @s_csq;
            if ($s_csq[0] ne 'missense_variant' and 
                $s_csq[0] ne 'protein_altering_variant' and 
                $s_csq[0] !~  /^inframe_(inser|dele)tion$/
            ){
                next;
            }
            my $desc = "same AA altered";
            if ($aa eq $csq->{amino_acids}){
                $desc = "same AA change";
            }
            $search_handles{hgmd_id}->execute($hgmd_id) 
              or die "Error searching 'HGMD' table in '$opts{t}': " . 
              $search_handles{hgmd_id} -> errstr;
            while (my ($vc, $disease) = $search_handles{hgmd_id}->fetchrow_array()) {
                  push @aa_matches, "$desc:HGMD_$hgmd_id:$vc:$disease:$hgvsc:$hgvsp";
                  $path++ if $vc =~ /^D[MP]$/;
            }
        }
    }
    push @results, join("\n", @aa_matches);
    return \@results, $path;
}
    
###########################################################
sub addClinvarMatches{
    my $var = shift;
    my $csq = shift;
    my $most_damaging_csq = shift;
    my @c_matches = getClinvarMatches($var);
    my @results = ();
    my $path = 0;
    if (@c_matches){
        $path++ if ( grep { $_->{clinical_significance} =~ /Pathogenic/ } @c_matches);
        foreach my $f 
        ( qw /
                measureset_id
                clinical_significance 
                all_traits 
                conflicted 
            /
        ){#add ClinVar results to row, multiple values per field separated by commas
            my $s = join(",", map { $_->{$f} } @c_matches );
            push @results, $s;
        }
    }else{
        push @results, map {"-"} (1..4); 
    }

    #get variants with same AA altered
    my @aa_matches = ();
    if ( exists $search_handles{clinvar_aa} and 
         exists $search_handles{clinvar_id} and  
         ($most_damaging_csq eq 'missense_variant' or 
          $most_damaging_csq eq 'protein_altering_variant' or 
          $most_damaging_csq =~  /^inframe_(inser|dele)tion$/ 
         )
        ){
        $search_handles{clinvar_aa}->execute
        (
            $csq->{feature},
            $csq->{protein_position},
        ) or die "Error searching 'ClinVar_VEP' table in '$opts{t}': " . 
          $search_handles{clinvar_aa} -> errstr;
        while (my @db_row = $search_handles{clinvar_aa}->fetchrow_array()) {
            my ($clinvar_id, $feature, $consequence, $protein_position, $aa, $hgvsc, $hgvsp) = @db_row; 
            next if (grep {$_->{measureset_id} eq $clinvar_id} @c_matches );
            my @s_csq = split("&", $consequence); 
            @s_csq = sort { $so_ranks{$a} <=> $so_ranks{$b} } @s_csq;
            if ($s_csq[0] ne 'missense_variant' and 
                $s_csq[0] ne 'protein_altering_variant' and 
                $s_csq[0] !~  /^inframe_(inser|dele)tion$/
            ){
                next;
            }
            my $desc = "same AA altered";
            if ($aa eq $csq->{amino_acids}){
                $desc = "same AA change";
            }
            $search_handles{clinvar_id}->execute($clinvar_id) 
              or die "Error searching 'ClinVar' table in '$opts{t}': " . 
              $search_handles{clinvar_id} -> errstr;
            while (my ($path, $clinsig, $conflicted, $disease) = 
                    $search_handles{clinvar_id}->fetchrow_array()
            ) {
                push @aa_matches,  "$desc:ClinVar_$clinvar_id:$clinsig:$disease:$hgvsc:$hgvsp";
                $path++ if $clinsig =~ /Pathogenic/; 
            }
        }
    }
    push @results, join("\n", @aa_matches);
    return \@results, $path;
}

###########################################################
sub getEtScore{
    #if more than one AA is altered return highest score
    my @scores = ();
    my $csq = shift;
    return '' if not $csq->{protein_position};
    my ($p_start, $p_end) = split("-", $csq->{protein_position}); 
    $p_end ||= $p_start;
    my @aa = split("/", $csq->{amino_acids});
    return '' if @aa != 2;
    s/-// for @aa; #turn deleted AAs to empty strings
    for (my $i = $p_start; $i <= $p_end; $i++){
        my $j = $p_start - $i;
        last if ($j > length($aa[1]));
        my $wt = substr($aa[0], 0, 1); 
        my $mut = substr($aa[1], 0, 1); 
        if ($wt ne $mut){
            $search_handles{et}->execute
            (
                $csq->{feature}, 
                $i, 
                $wt,
                $mut,
            ) or die "Error searching 'EvolutionaryTrace' table in '$opts{t}': " . 
              $search_handles{et} -> errstr;
            my $row = 0;
            while (my $score = $search_handles{et}->fetchrow_array()){
                push @scores, $score if $score;
                # multiple rows are possible because ensembl transcripts 
                # frequently map to more than one RefSeq protein
                #if ($row == 1){
                #    informUser("WARNING: Only expected one row from EvolutionaryTrace search!\n");
                #}
                $row++;
            }
        }
    }
    my $max = max(@scores); 
    return $max if $max;
    return '';
}
###########################################################
sub checkExac{
    my $var = shift;
    return -1 if not keys %exac_search_args;
    my @exac_hits = VcfReader::searchForPosition
    (
        %exac_search_args,
        chrom => $var->{CHROM},
        pos   => $var->{POS}
    );
    my @ids = ();  
    my $freq = 0;
    my $pop = '';
    foreach my $exac_line (@exac_hits) {
        #check whether the exac line(s) match our variant
        my @exac_split = split( "\t", $exac_line );
        if ( my $match = checkSnpMatches( $var, \@exac_split ) ) {
            foreach my $ex (@exac_pops){
                my $an = VcfReader::getVariantInfoField(\@exac_split, "AN_$ex");
                next if not $an;
                my $acs = VcfReader::getVariantInfoField(\@exac_split, "AC_$ex");
                next if not $acs;
                my $count = (split ",", $acs)[$match-1];
                my $this_freq = sprintf("%g", $count/$an);
                if ($this_freq > $freq){
                    $freq = $this_freq;
                    $pop = $ex;
                }
            }
        }
    }
    return ($freq, $pop);
}

###########################################################
sub checkEvs{
    my $var = shift;
    return -1 if not keys %evs_search_args;
    my @evs_hits = VcfReader::searchForPosition
    (
        %evs_search_args,
        chrom => $var->{CHROM},
        pos   => $var->{POS}
    );
    my @ids = ();  
    my $freq = 0;
    my $pop = '';
    foreach my $evs_line (@evs_hits) {
        #check whether the evs line(s) match our variant
        my @evs_split = split( "\t", $evs_line );
        if ( my $match = checkSnpMatches( $var, \@evs_split ) ) {
            foreach my $k (keys %evs_acs){
                my $counts = VcfReader::getVariantInfoField( \@evs_split, $k);
                my @af = split(",", $counts); 
                if (@af < $match){
                    die 
"Not enough allele frequencies found for allele counts for EVS line:\n$evs_line\n";
                }
                my $total = sum (@af);
                next if not $total;
                my $this_freq = sprintf("%g", $af[$match -1] / $total) ;
                if ($this_freq > $freq){
                    $freq = $this_freq;
                    $pop = $evs_acs{$k};
                }
            }
        }
    }
    return ($freq, $pop);
}

###########################################################
sub getCaddScore{
    return '.' if not @cadd_iters;
    my $var = shift;
    foreach my $ti (@cadd_iters){
        my $iter = $ti->query
        (
            "$var->{CHROM}:$var->{POS}-" . ($var->{POS} + 1)
        );
        while (my $result = $iter->next){
            chomp($result);
            my @res = split("\t", $result);
            next if $res[0] ne $var->{CHROM};
            my ($pos, $ref, $alt) = reduceRefAlt($res[1], $res[2], $res[3]);
            next if ($var->{POS} != $pos);
            next if ($var->{REF} ne $ref); #should we error here?
            next if ($var->{ALT} ne $alt); #diff alt allele
            return $res[5];
        } 
    }
    return '.';
}

###########################################################
sub reduceRefAlt{
    #reduce a single ref/alt pair to their simplest representation
    my ($pos, $ref, $alt) = @_;
    if (length($ref) > 1 and length($alt) > 1){
        #can only reduce if both REF and ALT are longer than 1
        my @r = split('', $ref);
        my @al = split('', $alt);
        while ($r[-1] eq $al[-1] and @r > 1 and @al > 1){
            #remove identical suffixes
            pop @r;
            pop @al;
        }
        while ($r[0] eq $al[0] and @r > 1 and @al > 1){
            #remove identical prefixes
            #increment position accordingly
            shift @r;
            shift @al;
            $pos++;
        }
        $ref = join('', @r);
        $alt = join('', @al);
    }
    return ($pos, $ref, $alt);
}

###########################################################
sub checkDbSnp{
    my $var = shift;
    return -1 if not keys %snp_search_args;
    my @snp_hits = VcfReader::searchForPosition
    (
        %snp_search_args,
        chrom => $var->{CHROM},
        pos   => $var->{POS},
    );
    my @ids = ();  
    my $freq = -1;
    foreach my $snp_line (@snp_hits) {
        #check whether the snp line(s) match our variant
        my @snp_split = split( "\t", $snp_line );
        if ( my $match = checkSnpMatches( $var, \@snp_split ) ) {
            my $id = VcfReader::getVariantField( \@snp_split, 'ID' );
            push @ids, $id if $id ne '.';
            my $c = VcfReader::getVariantInfoField( \@snp_split, 'CAF' );
            next if not $c;
            my @caf = split(",", $c); 
            next if $caf[$match] eq '.';
            $freq = $caf[$match] > $freq ? $caf[$match] : $freq;
        }
    }
    return ($freq, join(";", @ids) );
}

#################################################
sub checkSnpMatches {
    my ( $min_allele, $snp_line ) = @_;
    my %snp_min = VcfReader::minimizeAlleles($snp_line);
    foreach my $snp_allele ( keys %snp_min ) {
        (my $m_chr = $min_allele->{CHROM}) =~ s/^chr//;
        (my $s_chr = $snp_min{$snp_allele}->{CHROM}) =~ s/^chr//;
        next if $m_chr ne $s_chr;
        next if $min_allele->{POS} ne $snp_min{$snp_allele}->{POS};
        next if $min_allele->{REF} ne $snp_min{$snp_allele}->{REF};
        next if $min_allele->{ALT} ne $snp_min{$snp_allele}->{ALT};
        return $snp_allele;
    }
    return 0;
}

###########################################################
sub readEvsFile{
    return if not ($opts{e});
    my %info = checkDbFileAndGetInfo($opts{e});
    foreach my $field ( keys %evs_acs){
        if (not exists $info{$field}){
            die "ERROR: Could not find $field allele count field in $opts{e} EVS file!\n";
        }
    }
    return VcfReader::getSearchArguments($opts{e});
}

###########################################################
sub readExacFile{
    return if not ($opts{x});
    my %info = checkDbFileAndGetInfo($opts{x});
    foreach my $pop (@exac_pops){
        if (not exists $info{"AC_$pop"}){
            die "ERROR: Could not find AC_$pop allele count field in $opts{x} ExAC file!\n";
        }
        if (not exists $info{"AN_$pop"}){
            die "ERROR: Could not find AN_$pop allele number field in $opts{x} ExAC file!\n";
        }
    }
    return VcfReader::getSearchArguments($opts{x});
}

###########################################################
sub readDbSnpFile{
    return if not ($opts{s});
    my %info = checkDbFileAndGetInfo($opts{s});
    if (not exists $info{CAF}){
            die 
"ERROR - can't find CAF field in dbSNP file header for determining allele".
" frequency!\n";
    }
    return VcfReader::getSearchArguments($opts{s});
}

###########################################################
sub checkDbFileAndGetInfo{
    my $f = shift;
    my @head = VcfReader::getHeader($f);
    die "Header not ok for $f " 
        if not VcfReader::checkHeader( header => \@head );
    return VcfReader::getInfoFields( header => \@head);
}

###########################################################
sub getCaddIters{
    informUser("Checking CADD directory '$opts{z}'\n");
    opendir(my $CDIR, $opts{z}) or die "ERROR: Could not read directory '$opts{z}' - $!\n";
    my @c = map {"$opts{z}/$_"} grep { /\.(b)*gz$/ } readdir $CDIR;
    if (not @c){
        informUser
        (
            "WARNING: No .gz or .bgz files found in directory $opts{z}.\n"
        ); 
    }
    foreach my $cadd_file (@c){
        if (checkCaddFile($cadd_file)){
            push @cadd_iters, Bio::DB::HTS::Tabix->new(filename =>  $cadd_file);
        }
    }
    if (not @cadd_iters){
        informUser
        (
            "WARNING: No valid CADD files found in $opts{z}\n"
        );
    }else{
        informUser
        (
            scalar(@cadd_iters) . " valid CADD files found in $opts{z}\n"
        );
    }
}

###########################################################
sub checkCaddFile{
    my ($file) = @_;
    my $index = "$file.tbi";
    if (not -e $index ){
        informUser
        (
            "WARNING: Can't find tabix index for $file - please ensure $file ".
            "is bgzip compressed and indexed with tabix.\n"
        );
        return;
    }
    if ($file !~ /\.gz$/){
        informUser
        (
            "WARNING: CADD file $file does not have a '.gz' extension and is ".
            "presumably not bgzip compressed. Please ensure to use a ".
            "bgzip compressed CADD file.\n"
        );
        return;
    }
    my $FH = new IO::Uncompress::Gunzip $file or informUser 
    (
        "WARNING: IO::Uncompress::Gunzip failed while opening $file ".
        "for reading: \n$GunzipError"
    ) and return; 
    while (my $line = <$FH>){
        if ($line =~ /^##/){
            next;
        }elsif ($line =~ /^#/){
            if ($line !~ /^#Chrom\tPos\tRef\tAlt\t\w+Score\tPHRED/i){
                informUser
                (
                    "WARNING: Invalid header line for $file:\n\n$line\n\n".
                    "Expected to find a header as follows:\n\n".
                    "#Chrom\tPos\tRef\tAlt\tRawScore\tPHRED\n"
                ) and return;
            }else{
                return 1;
            }
        }else{
            last;
        }
    }
    informUser("WARNING: No header found for CADD file $file!\n");
    return;
}

###########################################################
sub checkCoverageDirs{
    foreach my $c (@coverage_dirs){
        informUser("Checking coverage directory '$c'\n");
        opendir(my $CDIR, $c) or die "ERROR: Could not read directory '$c' - $!\n";
        foreach my $f (readdir $CDIR){
            my $sample;
            my $is_summary = 0;
            if ($f =~ /^depth\.([\w\-_]+)[-_](S\d+|[ACTG]{10})[-_].*bqsr(\.sample_summary)*$/){
                $sample = $1;
                if ($3){
                    $is_summary = 1;
                }
            }elsif ($f =~ /^depth\.(\w+)[-_\.].*(bqsr|recal)(\.sample_summary)*$/){
                $sample = $1;
                if ($3){
                    $is_summary = 1;
                }
            }
            if ($sample){
                if ($is_summary){
                    if (exists $sample_coverage{$sample}->{summary}){
                        informUser("WARNING: Duplicate coverage summary files found ".
                                   "for $sample - will use the first encountered ".
                                   "($sample_coverage{$sample}->{summary}).\n");
                    }else{
                        $sample_coverage{$sample}->{summary} = "$c/$f";
                    }
                }else{
                    if (exists $sample_coverage{$sample}->{depth}){
                        informUser("WARNING: Duplicate depth files found ".
                                   "for $sample - will use the first encountered ".
                                   "($sample_coverage{$sample}->{depth}).\n");
                    }else{
                        $sample_coverage{$sample}->{depth}   = "$c/$f";
                    }
                }
            }
        }
    }
    foreach my $s (@samples_to_analyze){
        if (not exists $sample_coverage{$s}){
            informUser("WARNING: No coverage data found for sample $s\n");
        }else{
            if (not exists $sample_coverage{$s}->{summary}  ){
                informUser("WARNING: Missing summary depth file for sample $s\n");
            }
            if (not exists $sample_coverage{$s}->{depth} and ($opts{r} or $opts{n}) ){
                informUser("WARNING: Missing per-coordinate depth file for sample $s\n");
            }
        }
    }
}

###########################################################
sub getFastqcSummaries{
    informUser("Checking FASTQC directory '$opts{q}'\n");
    opendir(my $FDIR, $opts{q}) or die "ERROR: Could not read directory '$opts{q}' - $!\n";
    foreach my $f (readdir $FDIR){
        if (-d "$opts{q}/$f" and $f =~ /^([\w\-_]+)[-_](S\d+|[ACTG]{10})[-_].*fastqc$/){
            my $s = $1;
            opendir(my $QDIR, "$opts{q}/$f") or die "ERROR: Could not read directory '$opts{q}/$f' - $!\n";
            if (grep {$_ eq 'summary.txt'} readdir $QDIR){
                push @{$fastqc{$s}}, "$opts{q}/$f/summary.txt";
            }
        }
    }
}

###########################################################
sub checkGeneInheritanceFile{
    return if not $opts{l};
    open (my $LIST, $opts{l}) or die "Can't open gene list '$opts{l}' for reading: $!\n";
    while (my $line = <$LIST>){
        next if $line =~ /^#/;
        chomp $line;
        next if not $line;
        my @split = split("\t", $line); 
        $gene_conditions{$split[0]}->{inheritance} = $split[1];
        $gene_conditions{$split[0]}->{reportable} = $split[2];
        $gene_conditions{$split[0]}->{condition} = $split[3];
    }
}
###########################################################
sub setConsequenceRanks{
    my @so_terms = qw /
        transcript_ablation  
        splice_acceptor_variant
        splice_donor_variant
        stop_gained
        frameshift_variant
        stop_lost
        start_lost
        transcript_amplification
        inframe_insertion
        inframe_deletion
        missense_variant
        protein_altering_variant
        splice_region_variant
        incomplete_terminal_codon_variant
        stop_retained_variant
        synonymous_variant
        coding_sequence_variant
        mature_miRNA_variant
        5_prime_UTR_variant
        3_prime_UTR_variant
        non_coding_transcript_exon_variant
        intron_variant
        NMD_transcript_variant
        non_coding_transcript_variant
        upstream_gene_variant
        downstream_gene_variant
        TFBS_ablation
        TFBS_amplification
        TF_binding_site_variant
        regulatory_region_ablation
        regulatory_region_amplification
        feature_elongation
        regulatory_region_variant
        feature_truncation
        intergenic_variant 
    /;
    my $n = 0;
    %so_ranks = map { $_ => $n++ } @so_terms; 
}
###########################################################
sub readTranscriptDatabase{
    return if not $opts{t};
    $dbh = DBI->connect("DBI:$driver:$opts{t}", {RaiseError => 1})
      or die "Could not connect to sqlite database '$opts{t}': " . DBI->errstr . "\n";
    my %tables = map {$_ => undef} $dbh->tables;
    foreach my $t ( qw / transcripts uniprot cdd / ){#essential tables
        if (not exists $tables{"\"main\".\"$t\""}){
            die "ERROR: Could not find table '$t' in $opts{t} - did you use dbCreator.pl to create this database?\n";
        }
    }
    my %missing_tables = ();
    foreach my $t ( qw / HGMD HGMD_VEP ClinVar ClinVar_VEP / ){#non-essential tables
        if (not exists $tables{"\"main\".\"$t\""}){
            print STDERR "WARNING: $t table not found in $opts{t} - will skip variant annotations for this table.\n";
            $missing_tables{$t}++;
        }
    }
    my $q = "SELECT EnsemblGeneID, EnsemblTranscriptID, EnsemblProteinID, RefSeq_mRNA, RefSeq_peptide, CCDS,  Uniprot, TranscriptRank FROM transcripts";
    #we could do lazy loading for when we need to rank transcripts
    #but for current uses this isn't really worth the effort
    my $all = $dbh->selectall_arrayref($q);
    foreach my $tr (@$all){
        $transcript_ranks{$tr->[0]}->{$tr->[1]} = $tr->[7];
        my $i = 3;
        my %xrefs = map {$_ => $i++} qw / refseq_mrna refseq_peptide ccds uniprot / ; 
        foreach my $x (keys %xrefs){
            if (defined $tr->[$xrefs{$x}]){
                $enst_xref{$tr->[1]}->{$x} = $tr->[$xrefs{$x}];
            }
        }
    }
    %search_handles = 
    (
        cdd     =>  $dbh->prepare 
        (
            qq{ select * FROM cdd 
                WHERE UniprotId == ? 
                and End >= ? 
                and Start <= ? 
            } 
        ),

        uniprot =>  $dbh->prepare 
        (
            qq{ select * FROM uniprot 
                WHERE UniprotId == ? 
                and End >= ? 
                and Start <= ? 
            } 
        ),
     );
     if (not exists $missing_tables{HGMD}){
        $search_handles{hgmd_pos} = $dbh->prepare
        (
            qq{ select hgmd_id, disease, variant_class, gene_symbol, hgvs 
                FROM HGMD
                WHERE chrom == ? 
                and pos == ? 
                and ref == ?
                and alt == ?
            }
        );
        $search_handles{hgmd_id} =  $dbh->prepare
        (
            qq{ select variant_class, disease FROM HGMD
                WHERE hgmd_id == ? 
            }
        );
     }
    
     if (not exists $missing_tables{HGMD_VEP}){

        $search_handles{hgmd_aa} =  $dbh->prepare
        (
            qq{ select hgmd_id, feature, consequence, protein_position, amino_acids, hgvsc, hgvsp
                 FROM HGMD_VEP
                WHERE feature == ? 
                and protein_position == ? 
            }
        );
        
        $search_handles{hgmd_hgvs} =  $dbh->prepare
        (
            qq{ select hgvsc, hgvsp
                 FROM HGMD_VEP
                WHERE hgmd_id == ? 
            }
        );
      }

     if (not exists $missing_tables{ClinVar}){
        $search_handles{clinvar_pos} =  $dbh->prepare
        (
            qq{ select measureset_id, pathogenic, 
                clinical_significance, all_traits, conflicted 
                FROM ClinVar
                WHERE chrom == ? 
                and pos == ? 
                and ref == ?
                and alt == ?
            }
        );

        $search_handles{clinvar_id} =  $dbh->prepare
        (
            qq{ select pathogenic, clinical_significance, conflicted, all_traits 
                FROM ClinVar
                WHERE measureset_id == ? 
            }
        );
      }

     if (not exists $missing_tables{ClinVar_VEP}){
        $search_handles{clinvar_aa} =  $dbh->prepare
        (
            qq{ select measureset_id, feature, consequence, protein_position, amino_acids, hgvsc, hgvsp
                FROM ClinVar_VEP
                WHERE feature == ? 
                and protein_position == ? 
            }
        );


        $search_handles{clinvar_hgvs} =  $dbh->prepare
        (
            qq{ select hgvsc, hgvsp
                FROM ClinVar_VEP
                WHERE measureset_id == ? 
            }
        );
    }

    if (exists $tables{"\"main\".\"EvolutionaryTrace\""}){
        $search_handles{et} = $dbh->prepare
        (
            qq{ select score FROM EvolutionaryTrace
                WHERE EnsemblTranscriptID == ? 
                and Position == ? 
                and WildTypeResidue == ? 
                and MutantResidue == ? 
            }
        );
    }

}

###########################################################
sub getClinvarMatches{
    return if not  $search_handles{clinvar_pos} or not $search_handles{clinvar_hgvs};
    my $var = shift;
    my %clinvar = ();#keys are HGMD fields, values are array refs of values
    my @results = ();
    my @clinvar_fields = 
    ( qw /
            measureset_id
            pathogenic
            clinical_significance 
            all_traits 
            conflicted 
        /
    );
    $search_handles{clinvar_pos}->execute
    (
        $var->{CHROM}, 
        $var->{POS}, 
        $var->{REF}, 
        $var->{ALT}, 
    ) or die "Error searching 'ClinVar_VEP' table in '$opts{t}': " . 
      $search_handles{clinvar_pos} -> errstr;
    
    while (my @db_row = $search_handles{clinvar_pos}->fetchrow_array()) {
        my %clinvar = ();
        for (my $i = 0; $i < @db_row; $i++){
            $clinvar{$clinvar_fields[$i]} = $db_row[$i];
        }
        $search_handles{clinvar_hgvs}->execute($db_row[0]) 
          or die "Error searching 'ClinVar' table in '$opts{t}': " . 
          $search_handles{clinvar_hgvs} -> errstr;
        ($clinvar{hgvsc}, $clinvar{hgvsp}) = $search_handles{clinvar_hgvs}->fetchrow_array();
        push @results, \%clinvar;
    }
    return @results;
}

###########################################################
sub getHgmdMatches{
    return if not  $search_handles{HGMD} or not $search_handles{hgmd_hgvs};
    my $var = shift;
    my @results = ();
    my @hgmd_fields = 
    ( qw /
            hgmd_id
            disease
            variant_class
            gene_symbol
            hgvs
        /
    );
    $search_handles{hgmd_pos}->execute
    (
        $var->{CHROM}, 
        $var->{POS}, 
        $var->{REF}, 
        $var->{ALT}, 
    ) or die "Error searching 'HGMD' table in '$opts{t}': " . 
      $search_handles{hgmd_pos} -> errstr;
    while (my @db_row = $search_handles{hgmd_pos}->fetchrow_array()) {
        my %hgmd = ();
        for (my $i = 0; $i < @db_row; $i++){
            $hgmd{$hgmd_fields[$i]} = $db_row[$i];
        }
        $search_handles{hgmd_hgvs}->execute($db_row[0]) 
          or die "Error searching 'HGMD' table in '$opts{t}': " . 
          $search_handles{hgmd_hgvs} -> errstr;
        ($hgmd{hgvsc}, $hgmd{hgvsp}) = $search_handles{hgmd_hgvs}->fetchrow_array();
        push @results, \%hgmd;
    }
    return @results;
}

###########################################################
sub getHeaders{
    my %h = ();

    @{$h{Variants}} =  ( 
        qw/
            Sample
            GT
            AD
            AB
            GQ
            ExpectedInheritance
            AssociatedConditions
            Reportable\/Non-reportable
            Symbol
            Feature
            Consequence 
            HGVSc 
            HGVSp 
            Exon
            Intron
            Uniprot\/CDD_Features
            CDD_Feature_Residues
            LoF
            LoF_Filter
            LoF_info
            LoF_flags
            Polyphen
            SIFT
            MaxEntScan_diff
            MaxEntScan_alt
            MaxEntScan_ref
            EvolutionaryTrace
            CADD
            Hgmd_ID
            HGMD_Disease
            HGMD_variant_class
            HGMD_Symbol
            HGMD_HGVS
            HGMD_similar
            ClinVar_ID
            ClinVarSig
            ClinVarTrait
            ClinVarConflicted
            ClinVar_similar
            dbSNP_AF
            dbSNP_ID
            EVS_AF
            EVS_Pop
            ExAC_AF
            ExAC_Pop
            Chrom
            Pos
            Ref
            Alt
            Qual
            Filter
            AF
            AC
            AN
        /
    );

    return %h;
}
#################################################
sub informUser{
    my $msg = shift;
    my $time = strftime( "%H:%M:%S", localtime );
    if ($progressbar){
        $progressbar->message( "[INFO - $time] $msg" );
    }else{
        print STDERR "[INFO - $time] $msg";
    }
}
###########################################################

=head1 NAME

sampleSummarizer.pl - produce sample summary file for each sample in a VCF

=head1 SYNOPSIS

        sampleSummarizer.pl -i [vcf file] -t [transcript database] [options]
        sampleSummarizer.pl -h (display help message)
        sampleSummarizer.pl -m (display manual page)

=cut

=head1 ARGUMENTS

=over 8

=item B<-i    --input>

Input VCF file. Required.

=item B<-t    --transcript_database>

sqlite database of transcripts and protein info created using dbCreator.pl.
Required.

=item B<-o    --output_prefix>

Optional output prefix. Output files will be named "prefix[samplename].xlsx".
Default is to append sample name and .xlsx extention to input filename.

=item B<-c    coverage>

Directories with depth of GATK coverage data for each sample.

When reading the provided coverage directory, this program will get the sample
names from the VCF and look for corresponding files that match 

    'depth.[samplename][-_]S[0-9][anything].bqsr' or
'depth.[samplename][-_][tenLetterBarcode][anything]bqsr' 

and the corresponding files ending with '.sample_summary'. 

e.g.  depth.Sample1_S1_indelrealign.rmdups.bqsr
depth.Sample2_ACGCGGACTA_indelrealign.rmdups.bqsr

=item B<-q    --fastqc_dir>

Directory with samples fastqc results.

When looking for subdirectories of the provided fastqc directory, this program
will look for the summary.txt files within directories named as 

    [samplename]_S[0-9]_L[lanenumber]_R[readnumber]_[numbers].fastqc or
    [samplename]_[tenLetterBarcode]_L[lanenumber]_R[readnumber]_[numbers].fastqc

e.g.  Sample1_S1_L001_R1_001_fastqc Sample2_ACGCGGACTA_L001_R1_001_fastqc

=item B<-l    --gene_list>

.tsv file of gene names, expected inheritence and condition. The associated
inheritance pattern and conditions will be provided in output if this file is
provided.

=item B<-r    --reportable_regions>

Bed of regions in reportable genes for calculating coverage.

=item B<-n    --non_reportable_regions>

Bed file of regions in non-reportable genes for calculating coverage.

=item B<-s    --snp>

dbSNP VCF for annotating and filtering on allele frequency.

=item B<-x    --exac>

ExAC VCF for annotating and filtering on allele frequency.

=item B<-e    --evs>

Evs VCF for annotating and filtering on allele frequency. This must be a single
VCF made by combining the per-chromosome VCFs available for download.

=item B<-z    --cadd_dir>

Directory containing tabix indexed CADD scores for annotating CADD scores to
output.

=item B<-f    --allele_frequency>

Optional allele frequency cutoff for dbsnp, evs, exac. Variants with an allele
frequency above this value will be removed from the output. Valid values are
between 0.00 and 1.00.

=item B<-a    --AF>

Optional allele frequency cutoff for within VCF. Alleles with AF INFO fields 
higher than this value will be filtered.

=item B<-g    --gq>

Optional minimum genotype quality (GQ) for calls. Sample genotypes will only be
included in the output if greater than this value. Default = 0.

=item B<-d    --depth>

Optional minimum depth for sample calls. Genotype calls with a depth lower than
this will not be included in output.

=item B<-b    allele_balance>

Optional min and optional max alt allele ratio per sample call. If one value is
provided this will act as the minimum allele balance cutoff. If a second value
is provided this will be used as a maximum allele balance cutoff (if only
looking for heterozygous changes, for example). Valid values between 0.00 and
1.00.

=item B<--samples>

Only analyze the samples listed here. 

=item B<-u    --summary>

Summary XLSX file for outputting variants that are consistent with the 
proposed model of inheritance specified by the -l/--gene_list file for each 
sample.

=item B<--cadd_cutoff>

Only consider variants with a cadd score equal to or greater than this value
for outputting into the -u/--summary XLSX file.

=item B<--pathogenic>

Only output variants if they're disease causing in clinvar/hgmd, LoF or affect 
CDD feature residues.

=item B<-p    --progress>

Use this flag to show a progress bar while assessing variants.

=item B<-h    --help>

Display help message.

=item B<-m    --manual>

Show manual page.

=back 

=cut

=head1 DESCRIPTION

This program reads a VCF and produces a summary file for each sample within it,
listing variants and annotations and optionally creating summaries for depth and
fastqc reports. This requires a transcript database for the genes of interest
created using dbCreator.pl and various other database files as detailed in the
options section. 

It is generally recommended not to run this manually but to allow the commands
created by make_alignment_commands_mopd.pl to run this program as it will
automatically identify the required files. 

Separate sheets will be created for 'functional' variants (i.e. those generally
expected to effect the function of a transcript or its encoded protein) and for
'other' variants that do not fit in this class. If the necessary files are
provided a summary of the depth of coverage and coverage of individual regions
will be provided also.

=cut

=head1 AUTHOR

David A. Parry

=head1 COPYRIGHT AND LICENSE

Copyright 2016, David A. Parry

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

