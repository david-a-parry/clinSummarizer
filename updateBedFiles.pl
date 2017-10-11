#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($RealBin);
use Getopt::Long;

my %opts = (b => "GRCh37", g => "genes/gene_inheritance_and_diseases.txt");
GetOptions
(
    \%opts,
    'g|gene_list=s',
    'b|build=s',
    'u|ucsc',
    'h|?|help',
) or usage( "Syntax error." );
usage() if $opts{h};

sub usage{
    my $msg = shift;
    print STDERR "$msg\n" if $msg;
    print STDERR <<EOT

Usage: $0 [options]

Options:

    -b BUILD --build BUILD
        Genome build to use. Default = GRCh37
        Valid genome versions: hg19, hg38, GRCh37, GRCh38

    -g GENELIST --gene_list GENELIST
        File containing list of genes to retrieve coordinates for.
        Default = "genes/gene_inheritance_and_diseases.txt"

    -u    --ucsc
        Use UCSC genome browser for retrieving gene coordinates. 
        Useful if the Ensembl REST server is down.

Example: $0 
         $0 -b hg19 
         $0 -b hg19 -g gene_file.txt


EOT
    ;
    exit 2 if $msg;
    exit;
}


if ($opts{b} !~ /^(GRCh3[78]|hg(19|38))$/){
    die "Unrecognized build '$opts{b}' given\n";
}

my %genes = ();
my %non_coding = ();
getGenes();

getCodingExons();

if ($opts{u}){
    getGeneRegionsUcsc();
}else{
    getGeneRegions();
}

print STDERR "Done updating BED files.\n";

##################################################
sub getGeneRegions{
    foreach my $r (qw /reportable non-reportable/){
        next if not @{$genes{$r}};
        print STDERR "Retrieving genomic regions for $r genes.\n";
        my $cmd = "perl $RealBin/genomeUtils/coordinatesFromGenes.pl -j -g ". 
                  join(" ", @{$genes{$r}});
        if ($opts{b} eq 'hg19' or $opts{b} eq 'GRCh37'){
            $cmd .= ' -r ';
        }
        if ($opts{b} =~ /^hg/){
            $cmd .= "| perl -wne \'chomp; print \"chr\$_\\n\"\'";
        }
        my $output = executeCommand($cmd);
        my $bed = "bed_files/$r". "_genes_$opts{b}.bed";
        $bed =~ s/non-reportable/not_reportable/;
            #legacy reasons - we called the bed files 'not_reportable'
            #instead of 'non-reportable'
        open (my $OUT, ">", $bed) or die "Could not open $bed for writing: $!\n";
        print $OUT $output;
        close $OUT;
    }
}

##################################################
sub getGeneRegionsUcsc{
    foreach my $r (qw /reportable non-reportable/){
        next if not @{$genes{$r}};
        print STDERR "Retrieving gene coordinates for $r genes from UCSC.\n";
        my $g = join(" ", @{$genes{$r}}) ; 
        my $cmd = "perl $RealBin/getExonsFromUcsc/getExonsFromUcsc.pl -m -k -f 500 -w -g $g";
        if ($opts{b} eq 'GRCh38' or $opts{b} eq 'hg38'){
            $cmd .= " -b hg38";
        } 
        if ($opts{b} =~ /^GRC/){
            $cmd .= " | grep -vP '^chr\w+_' | sed s/chr// ";
        }
        my $output = executeCommand($cmd);
        my $bed = "bed_files/$r". "_genes_$opts{b}.bed";
        $bed =~ s/non-reportable/not_reportable/;
            #legacy reasons - we called the bed files 'not_reportable'
            #instead of 'non-reportable'
        open (my $OUT, ">", $bed) or die "Could not open $bed for writing: $!\n";
        print $OUT $output;
        close $OUT;
    }

}

##################################################
sub getCodingExons{
    foreach my $r (qw /reportable non-reportable/){
        next if not @{$genes{$r}};
        print STDERR "Retrieving coding exons for $r genes.\n";
        my $g = join(" ", @{$genes{$r}}) ; 
        if ($non_coding{$r}){
            $g = join(" ", grep{ !$non_coding{$r}->{$_} } @{$genes{$r}});
        }
        my $cmd = "perl $RealBin/getExonsFromUcsc/getExonsFromUcsc.pl -m -k -c -g $g";
        if ($opts{b} eq 'GRCh38' or $opts{b} eq 'hg38'){
            $cmd .= " -b hg38";
        } 
        $cmd .= " | sed s/chr// " if $opts{b} =~ /^GRC/;
        my $output = executeCommand($cmd);
        if ($non_coding{$r}){
            print STDERR "Appending exons of non-coding genes for $r genes.\n";
            my $ncg = join(" ", sort keys(%{$non_coding{$r}}));
            my $nc_cmd = "perl $RealBin/getExonsFromUcsc/getExonsFromUcsc.pl -m -k -g $ncg";
            if ($opts{b} eq 'GRCh38' or $opts{b} eq 'hg38'){
                $nc_cmd .= " -b hg38";
            } 
            my $nc_output = executeCommand($nc_cmd);
            $output .= $nc_output;
        }
        my $bed = "bed_files/$r". "_coding_exons_$opts{b}.bed";
        $bed =~ s/non-reportable/not_reportable/;
            #legacy reasons - we called the bed files 'not_reportable'
            #instead of 'non-reportable'
        open (my $OUT, ">", $bed) or die "Could not open $bed for writing: $!\n";
        print $OUT $output;
        close $OUT;
    }

}
##################################################
sub getGenes{
    open (my $GENES, "<", $opts{g}) or die "Could not open $opts{g}: $!\n";
    my $n = 0;
    $genes{reportable} = [];
    $genes{"non-reportable"} = [];
    while (my $line = <$GENES>){
        next if $line =~ /^#/;
        $n++;
        chomp $line;
        my ($s, undef, $r, undef, $nc) = split(/\t/, $line);
        next if not $s; #empty line?
        $s =~ s/^\s+//;#remove preceding whitespace
        $s =~ s/\W.+$//;#remove trailing crud
        if (lc($r) ne 'reportable' and lc($r) ne 'non-reportable'){
            die <<EOT 
Expected either 'reportable' or 'non-reportable' in column 3, line $n of 
$opts{g} but found '$r'. Full line was:
$line
Please edit this file and try again.
EOT
            ;
        }
        push @{$genes{lc$r}}, $s;
        if (defined $nc and lc($nc) eq 'non-coding'){
            $non_coding{lc$r}->{$s}++;
        }
    }
    die "No reportable genes identified from $opts{g}!\n" if not @{$genes{reportable}};
    print STDERR "Identified ". scalar(@{$genes{reportable}}) . " reportable genes from $opts{g}.\n";
    print STDERR "Identified ". scalar(@{$genes{"non-reportable"}}) . " non-reportable genes from $opts{g}.\n";
}

##################################################
sub executeCommand{
    my $cmd = shift;
    print STDERR "Executing: $cmd\n";
    my $output = `$cmd`;
    checkExit($?);
    if (not $output){
        die "ERROR: No output from command '$cmd'!\n";
    }
    return $output;
}

###################################################
sub checkExit{
    my $e = shift;
    my $msg = "";
    if($e == -1) {
        $msg = "Failed to execute: $!\n";
    }elsif($e & 127) {
        $msg = printf "Child died with signal %d, %s coredump\n",
        ($e & 127),  ($e & 128) ? 'with' : 'without';
    }elsif($e != 0){
        $msg = printf "Child exited with value %d\n", $e >> 8;
    }
    die $msg if $msg;
}
