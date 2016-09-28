#!/usr/bin/env perl
use strict;
use warnings;

my $build = shift;
$build = "GRCh37" if (not $build);
if ($build !~ /^(GRCh3[78]|hg(19|38))$/){
    die "Unrecognized build '$build' given\n";
}
my $gene_list = "genes/gene_inheritance_and_diseases.txt";

my %genes = getGenes();

getCodingExons();

getGeneRegions();

print STDERR "Done updating BED files.\n";

##################################################
sub getGeneRegions{
    foreach my $r (qw /reportable non-reportable/){
        next if not @{$genes{$r}};
        print STDERR "Retrieving genomic regions for $r genes.\n";
        my $cmd = "perl genomeUtils/coordinatesFromGenes.pl -j -g ". 
                  join(" ", @{$genes{$r}});
        if ($build eq 'hg19' or $build eq 'GRCh37'){
            $cmd .= ' -r ';
        }
        if ($build =~ /^hg/){
            $cmd .= "| perl -wne \'chomp; print \"chr\$_\\n\"\'";
        }
        my $output = executeCommand($cmd);
        my $bed = "bed_files/$r". "_genes_$build.bed";
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
        my $cmd = "perl getExonsFromUcsc/getExonsFromUcsc.pl -m -k -c -g ". 
                   join(" ", @{$genes{$r}}) ; 
        if ($build eq 'GRCh38' or $build eq 'hg38'){
            $cmd .= " -b hg38";
        } 
        $cmd .= " | sed s/chr// " if $build =~ /^GRC/;
        my $output = executeCommand($cmd);
        my $bed = "bed_files/$r". "_coding_exons.bed";
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
    open (my $GENES, "<", $gene_list) or die "Could not open $gene_list: $!\n";
    my %symbols;
    my $n = 0;
    $symbols{reportable} = [];
    $symbols{"non-reportable"} = [];
    while (my $line = <$GENES>){
        $n++;
        chomp $line;
        my ($s, undef, $r) = split(/\t/, $line);
        next if not $s; #empty line?
        $s =~ s/^\s+//;#remove preceding whitespace
        $s =~ s/\W.+$//;#remove trailing crud
        if (lc($r) ne 'reportable' and lc($r) ne 'non-reportable'){
            die <<EOT 
Expected either 'reportable' or 'non-reportable' in column 3, line $n of 
$gene_list but found '$r'. Full line was:
$line
Please edit this file and try again.
EOT
            ;
        }
        push @{$symbols{lc$r}}, $s;
    }
    die "No reportable genes identified from $gene_list!\n" if not @{$symbols{reportable}};
    print STDERR "Identified ". scalar(@{$symbols{reportable}}) . " reportable genes from $gene_list.\n";
    print STDERR "Identified ". scalar(@{$symbols{"non-reportable"}}) . " non-reportable genes from $gene_list.\n";
    return %symbols;
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
