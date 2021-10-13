#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long;

my ($inFile, $outDir, $help, $verbose); 

($inFile, $outDir) = ("data/sequenceList.tsv", "data/sequences");

&GetOptions( 
    "h|help"             => \$help,
    "v|verbose"          => \$verbose,
    "i|seqids=s"         => \$inFile,
    "o|outdir=s"         => \$outDir
    );

if( $help ) {
    &help();
    exit(1);
}

if (not -e $inFile or not -e $outDir){
    print "Either inFile [$inFile] or output directory [$outDir] does not exist!\n";
    &help();
    die(); 
}

print "Fetching sequences in [$inFile], writing these to [$outDir].\n" if (defined($verbose));

my $ids = readInFile($inFile); 
fetchSeqs($ids, $outDir);
reportStats($outDir);


exit(0); 

######################################################################

sub readInFile {
    
    my $inFile = shift;
    my @ids; 

    print "Opening [$inFile]\n" if (defined($verbose));
    open(IN,  "< $inFile") || die "FATAL: cannot open [$inFile]\n[$!]";
    while(my $in = <IN>){
	next if $in =~ m/^#/;
	chomp($in);
	if($in =~ /^(\S+)/  ){
	    push(@ids, $1);
	    print "Read sequence ID: [$1]\n" if (defined($verbose));
	}
    }
    return \@ids; 
}

######################################################################

sub fetchSeqs {
    my ($ids, $outDir) = @_;
    
    foreach my $id (@{ $ids } ){ 
	my $range;
	if($id =~ /(\S+):(\d+)\.\.(\d+)/){
	    $id = $1;
	    $range = "$2..$3";
	}

	my $outFileFasta = $outDir . "/" . $id . ".fasta";
	my $outFileEmbl  = $outDir . "/" . $id . ".embl";
	print "Fetching [$id], writing fasta to [$outFileFasta] & [$outFileEmbl].\n" if (defined($verbose));
	
	    #https://www.ebi.ac.uk/training/online/courses/ena-quick-tour/searching-and-visualising-data-ena/
	    #curl -X GET "https://www.ebi.ac.uk/ena/browser/api/fasta/BK006945.2?download=true&lineLimit=0&expanded=false&gzip=false&range=455933..457732" -H "accept: text/plain"
	if (not -e $outFileFasta){
	    my $cmd = "curl -o $outFileFasta -X GET \42https://www.ebi.ac.uk/ena/browser/api/fasta/$id?download=true&lineLimit=0&expanded=true&gzip=false\42 -H \42accept: text/plain\42";
	    $cmd . "&range=" . $range if (defined($cmd) && defined($range));
	    print "Running [$cmd].\n" if (defined($verbose));	    
	    system($cmd) and die("FATAL: failed to execute [$cmd]!\n[$!]");
	}
	else {
	    print "[$outFileFasta] exists, not refetching!\n" if (defined($verbose));
	}

	if (not -e $outFileEmbl){
	    my $cmd = "curl -o $outFileEmbl -X GET \42https://www.ebi.ac.uk/ena/browser/api/embl/$id?download=true&lineLimit=0&expanded=true&gzip=false\42 -H \42accept: text/plain\42";
	    $cmd . "&range=" . $range if (defined($cmd) && defined($range));
	    print "Running [$cmd].\n" if (defined($verbose));	    
	    system($cmd) and die("FATAL: failed to execute [$cmd]!\n[$!]");
	}
	else {
	    print "[$outFileEmbl] exists, not refetching!\n" if (defined($verbose));
	}

    }
    return 0;    
}
#######################################################################
sub seqStat {
        my $inFile = shift; 
        my %seqStat;  
        #Find number of sequences and total sequence length:
        my $runSeqstat = "esl-seqstat -ac --dna --informat fasta $inFile";
        open(SEQSTAT, "$runSeqstat |")  or die "FATAL: failed to open pipe: [$runSeqstat]\n[$!]";
        while (my $stat = <SEQSTAT>){
	    if ($stat =~ /^=\s+(\S+)/){
		$seqStat{'id'}=$1;
	    }
	    elsif ($stat =~ /^Total.*\s+(\d+)$/){
		$seqStat{'numRes'}=$1;
	    }
	    elsif ($stat =~ /^= (\S+)\s+(\d+)/){
		$seqStat{$1}=$2;
	    }
	    elsif ($stat =~ /^residue: (\S+)\s+\d+\s+(\S+)/){
		$seqStat{ 'freq' . $1}=$2;
	    }
        }
        close(SEQSTAT);
        return \%seqStat; 
}


######################################################################
sub reportStats {
    my $outDir = shift;
    my @fafiles  = glob( $outDir . '/*.fasta' );

    open(UT, "> $outDir/summary-stats.tsv" ) or die "FATAL: failed to open [$outDir/summary-stats.tsv] for writing!\n[$!]";
    printf UT "ENAID\tLength\tCG.content\tfilename\n";
    printf    "ENAID\tLength\tCG.content\tfilename\n" if (defined($verbose));
    foreach my $faFile (@fafiles){
	my %seqStat = %{ &seqStat($faFile) };
	printf UT $seqStat{'id'} . "\t" . $seqStat{'numRes'} . "\t%0.2f\t" . $faFile . "\n", $seqStat{ 'freqC' }+$seqStat{ 'freqG' };
	printf    $seqStat{'id'} . "\t" . $seqStat{'numRes'} . "\t%0.2f\t" . $faFile . "\n", $seqStat{ 'freqC' }+$seqStat{ 'freqG' } if (defined($verbose));
    }
    close(UT);

    my @emfiles  = glob( $outDir . '/*.embl' );
    system("embl2stats.pl head > $outDir/summary.txt") and die("FATAL: failed to run [embl2stats.pl head > summary.txt]!\n[$!]");
    foreach my $emFile (@emfiles){
	system("embl2stats.pl $emFile >> $outDir/summary.txt") and die("FATAL: failed to run [embl2stats.pl $emFile >> summary.txt]!\n[$!]");
    }
    system("R CMD BATCH bin/seqStats.R") and die("FATAL: failed to execute [R CMD BATCH bin/SeqStats.R]!\n[$!]");
}

######################################################################
sub help {
    print STDERR <<EOF;

fetchSeqs.pl: fetch sequences from the ENA database, report summary statistics for these.

Usage:   fetchSeqs.pl -i <seqfile> 
Options:       -h|--help                     Show this help.

               -i|--seqids <filename>        a TSV file, first column is ENA sequence IDs  
TODO:

EOF
}
