#!/usr/bin/perl

#compareNullSeqs.pl: Compare the randomized/null sequences to the original -- what features are preserved?
#           --nuc freqs, sequence similarity, coverage, #ORFs & length distributions, 


use warnings;
use strict;
use Getopt::Long;
use Time::HiRes 'gettimeofday', 'tv_interval';
my ($help, $verbose); 
my ($num,$nativeDir,$randomisedDir) = (10,"./","./"); 

&GetOptions(
    "n|num=s"            => \$num,
    "N|natdir=s"          => \$nativeDir,
    "r|randdir=s"         => \$randomisedDir,
    "h|help"             => \$help,
    "v|verbose"          => \$verbose
    );

if( $help ) {
    &help();
    exit(1);
}

#read sequence file names in:
my @nativeSeqFiles = glob("$nativeDir/*.fasta  $nativeDir/*.fa  $nativeDir/*.fna"); 
print "Native seqFiles: [@nativeSeqFiles]\n" if defined($verbose);

my @randomisedSeqFiles = glob("$randomisedDir/*.fasta  $randomisedDir/*.fa  $randomisedDir/*.fna"); 
print "Randomised seqFiles: [@randomisedSeqFiles]\n" if defined($verbose);
my %randomisedSeqFileHash=map{$_ => 0}@randomisedSeqFiles;

#What to return?
#-- Z-scores? (CASP)
#-- Empirical P-values?
#-- log-odds? 
#-- something else? 

#Paradigm is 1. native sequence vs N matched controls
#  

#Mono-stats:
#-- length
#-- mono-, di- & Tri- nuc frequencies: 
#-- Relative entropy of the sequences (\sum p*log2(p/q))

#For each native sequence loop
foreach my $seqFile (@nativeSeqFiles){
    my $natFileRoot = $seqFile; 
    if($seqFile =~ /$nativeDir\/(.*?)\.(fasta|fa|fna)$/){
	$natFileRoot = $1;
    }
    
    my @randomised=();

    #Find the matching null sequences: 
    foreach my $rf (@randomisedSeqFiles){
	if($rf =~ /\/$natFileRoot/){
	    push(@randomised, $rf);
	    print "[$rf]\n";
	}
    }
    
    
    print "seqFile: [$seqFile] natFileRoot: [$natFileRoot]\n" if defined($verbose);
    next if (scalar(@randomised) == 0); # no randomised sequences to compare too...

    
    my $nSeq = parseFastaSeq($seqFile);
    my %nativeFreqs = %{ &extract0123merFreqs($nSeq) };
    
    #For each corresponding null sequence loop
    for(my $fileI=0; $fileI<scalar(@randomised); $fileI++){
	
	my $randFileRoot = $randomised[$fileI]; 
	if($randomised[$fileI] =~ /$randomisedDir\/(.*?)\.(fasta|fa|fna)$/){
	    $randFileRoot = $1;
	}
	
	printf "$randFileRoot\n";

	#######################################################################
	my $nullStats = $randomised[$fileI] . ".nullstats";
	open(NUT, "> $nullStats") or die "FATAL: failed to open [$nullStats] for writing!\n[$!]\n";
	#Print null sequence stats loop...
	my ($rSeq, $rSeqCnt,$rName) = ("go", 1, "name"); 
	while(length($rSeq)){
	    
	    ($rName, $rSeq) = parseNthFastaSeq($randomised[$fileI], $rSeqCnt) if(-s $randomised[$fileI]);
	    #print "rSeqCnt:[$rSeqCnt] rSeq:[$rSeq]\n";
	    if (length($rSeq)){
		$rSeqCnt++;		
	    }
	    else{
		last;
	    }
	    
	    my %randomFreqs = %{ &extract0123merFreqs($rSeq) };

	    my $absDiffRat    = abs($nativeFreqs{0}{'len'} - $randomFreqs{0}{'len'})/$nativeFreqs{0}{'len'};
	    my $lenRelEntropy = 0.0;
	    $lenRelEntropy = $absDiffRat * log( $absDiffRat )/log(2) if($absDiffRat>0);
	    printf NUT "$rName\t0:len:%d/%d:%0.5f\t", $nativeFreqs{0}{'len'}, $randomFreqs{0}{'len'}, $lenRelEntropy;  
	    for (my $kmerLen=1; $kmerLen<=3; $kmerLen++){
		my $entropy=0;
		printf NUT "$kmerLen"; 
		foreach my $kmer (sort keys %{$nativeFreqs{$kmerLen}}) {
		    #printf "$kmerLen\t$kmer\t%0.4f\t%0.4f\n", $nativeFreqs{$kmerLen}{$kmer}, $randomFreqs{$kmerLen}{$kmer};
		    next if (not defined($nativeFreqs{$kmerLen}{$kmer}) or not defined($randomFreqs{$kmerLen}{$kmer}));
		    #printf ":$kmer:%0.5f/%0.5f", $nativeFreqs{$kmerLen}{$kmer}, $randomFreqs{$kmerLen}{$kmer};
		    next if (not isNumeric($nativeFreqs{$kmerLen}{$kmer}) or not isNumeric($randomFreqs{$kmerLen}{$kmer}) or
			     $nativeFreqs{$kmerLen}{$kmer} == 0 or $randomFreqs{$kmerLen}{$kmer} == 0);
		    $entropy+= $nativeFreqs{$kmerLen}{$kmer}*log( $nativeFreqs{$kmerLen}{$kmer} / $randomFreqs{$kmerLen}{$kmer} )/log(2);
		}
		printf NUT ":%0.5f\t", $entropy; 
	    }
	    printf NUT "\n";
	}
	close(NUT);
	#######################################################################

	#Sequence similarity -- based on 9-mers:  
	my $mashOut = $randomised[$fileI] . ".mashout";
	my $mashExe = "mash dist -w 1.0 -i -n -k 9 $seqFile $randomised[$fileI] > $mashOut 2> /dev/null";
	system( $mashExe ) and printf "WARNING: [$mashExe] failed to execute!\n\[$!]\n";
	######################################################################
	#ORF stats 
	
	
	#esl-translate -l 100 ./data/sequences-null/AB000109.1.phastsim-codon-001.fasta > blah && esl-seqstat -c blah
	
	######################################################################
	system("cat $nullStats  $mashOut ") if defined($verbose);
	#######################################################################
	
    }
    #last;
}


#Sequence similarity
#-- PID?
#-- longest identity
#-- bitscore

#Kmer based distances: 
#mash dist -w 1.0 -i -n -k 9  ./sequences/AM746676.1.fasta sequences-null/AM746676.1.rose-001.fasta
#reference-ID, query-ID, distance, p-value, shared-hashes
#ENA|AM746676|AM746676.1	rose10	0.000278824	7.37517e-12	995/1000

#Alignment
#minimap2 ./sequences/AM746676.1.fasta sequences-null/AM746676.1.rose-001.fasta
#1          2            3           4         5               6           7             8            9          10                 11           12
#queryName, queryLength, queryStart, queryEnd, relativeStrand, targetName, targetLength, targetStart, targetEnd, numResidueMatches, alnBlockLen, mappingQuality 
#rose10	13033779	9	13033772	+	ENA|AM746676|AM746676.1	13033779	9	13033772	12915052	13033763	60	tp:A:P	cm:i:2376604	s1:i:12915052	s2:i:8421	dv:f:0.0012
#rose8	13033779	9	13033772	+	ENA|AM746676|AM746676.1	13033779	9	13033772	12917437	13033763	60	tp:A:P	cm:i:2377725	s1:i:12917437	s2:i:8409	dv:f:0.0011

#ORF Module:
#esl-translate -W -l 100 ./data/sequences-null/AB000109.1.ushuffle-k10.fasta > blah && esl-seqstat blah

#Rfam & Pfam matches
#Find families for each sequence, search these against the controls -- sum bitscore differences? 








######################################################################
#parseFastaSeq: return first biosequence in a given fasta file: 
sub parseFastaSeq {
    my $inFile = shift;
    open(IN, "< $inFile") or die("FATAL: failed to open [$inFile]\n[$!]");
    my $name;
    my $seq='';
    while(my $in = <IN>){
	if($in =~ /^>(\S+)/ && not defined($name)){
	    $name=$1; 
	}
	elsif($in =~ /^>(\S+)/ && defined($name)){
	    last;
	}
	else {
	    chomp($in);
	    $seq .= $in; 
	}
    }
    close(IN);        
    
    return $seq; 
}

######################################################################
#parseNthFastaSeqs: return the Nth biosequence in a given fasta file: 
sub parseNthFastaSeq {
    my $inFile = shift;
    my $nthSeq = shift;

    die "FATAL: file:[$inFile] doesn't exist or is empty, or nthSeq[$nthSeq] undefined!\n"
	if(not -s $inFile or not defined($nthSeq) or not $nthSeq);
    open(IN, "< $inFile") or die("FATAL: failed to open [$inFile]\n[$!]");
    my ($name, $seq, $counter) =('','',1);
    while(my $in = <IN>){
	if($in =~ /^>\s*(\S+)/    && $nthSeq >  $counter){
	    $counter++;
	}
	elsif($in =~ /^>\s*(\S+)/ && $nthSeq == $counter){
	    $name=$1;
	    $counter++;
	    #print "parsing [$name], counter = [$counter], nthSeq = [$nthSeq]\n";
	}
	elsif($in =~ /^>\s*(\S+)/ && $nthSeq <  $counter){
	    $counter++;
	    last;
	}
	elsif ($nthSeq == ($counter-1)) {
	    chomp($in);
	    $seq .= $in; 
	}
    }
    close(IN);        
    
    return ($name,$seq); 
}


######################################################################

#extractNmerFreqs: return a pointer to a hash of 0- 1- 2- & 3-mer counts for an input sequence 
sub extract0123merFreqs {
    my $seq  = shift;
    #my $kmerLen = shift; 
    my $len = length($seq);
    
    my %kmerFreq;
    
    $kmerFreq{0}{'len'} = $len;

    for (my $kmerLen=1; $kmerLen<=3; $kmerLen++){
	
	for (my $i = 0; $i < $len+1-$kmerLen; $i++){
	    my $subseq = substr($seq, $i, $kmerLen);
	    if(defined($kmerFreq{$kmerLen}{$subseq})){
		$kmerFreq{$kmerLen}{$subseq}++;
	    }
	    else{
		$kmerFreq{$kmerLen}{$subseq}=1;
	    }
	}
	
	foreach my $kmer (keys %{$kmerFreq{$kmerLen}}) {
	    $kmerFreq{$kmerLen}{$kmer} = $kmerFreq{$kmerLen}{$kmer}/($len+1-$kmerLen);
	}
	
    }
    
    return \%kmerFreq; 
}

######################################################################
sub isNumeric {
    my $num = shift;
    if ($num=~/^-?\d+\.?\d*$/) { 
        return 1; 
    }
    else {
        return 0;
    }
}

