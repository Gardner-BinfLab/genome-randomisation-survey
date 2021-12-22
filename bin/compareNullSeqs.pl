#!/usr/bin/perl

#compareNullSeqs.pl: Compare the randomized/null sequences to the original -- what features are preserved?
#           --nuc freqs, sequence similarity, coverage, #ORFs & length distributions, 


use warnings;
use strict;
use Getopt::Long;
use Time::HiRes 'gettimeofday', 'tv_interval';
use Math::NumberCruncher;
my ($help, $verbose, $debug, $debugQuick); 
my ($nativeDir,$randomisedDir) = ("./","./"); 

&GetOptions(
    "N|natdir=s"         => \$nativeDir,
    "R|randdir=s"        => \$randomisedDir,
    "d|debug"            => \$debug,
    "dq|debugQuick"      => \$debugQuick,
    "h|help"             => \$help,
    "v|verbose"          => \$verbose
    );

if( $help ) {
    &help();
    exit(1);
}

#read sequence file names in:
my @nativeSeqFiles = glob("$nativeDir/*.fasta  $nativeDir/*.fa  $nativeDir/*.fna"); 
printf "Native seqFiles: [%d]\n", scalar(@nativeSeqFiles) if defined($verbose);

my @randomisedSeqFiles = glob("$randomisedDir/*.fasta  $randomisedDir/*.fa  $randomisedDir/*.fna"); 
printf "Randomised seqFiles: [%d]\n", scalar(@randomisedSeqFiles) if defined($verbose);
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

#ORF stats
# noSeqs, totRes, aveLen
# could use amino-acid frequencies too

#Sequence similarity 
# -k-mer similarity 


#For each native sequence...
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
	    #print "[$rf]\n";
	}
    }
    
    
    print "seqFile: [$seqFile] natFileRoot: [$natFileRoot]\n" if defined($verbose);
    next if (scalar(@randomised) == 0); # no randomised sequences to compare too...
    
    
    my $nSeq = parseFastaSeq($seqFile);
    my $nativeFreqs = &extract0123merFreqs($nSeq);
    
    #ORF stats:
    my $nOrfOut = $seqFile . ".orfout";
    unlink($nOrfOut) if (-s $nOrfOut);
    computeORFstats($seqFile, $nOrfOut); 
    
#    my $nOrfExe = "esl-translate -W --informat fasta -l 50 $seqFile > blah.translate.native && esl-seqstat -c blah.translate.native > $nOrfOut 2> /dev/null";
#    system( $nOrfExe ) and printf "WARNING: [$nOrfExe] failed to execute!\n\[$!]\n";
#    unlink( "blah.translate") if (not defined $debugQuick or not defined $debug);
    
    #For each corresponding null sequence loop file...
    for(my $fileI=0; $fileI<scalar(@randomised); $fileI++){
	
	my $randFileRoot = $randomised[$fileI]; 
	if($randomised[$fileI] =~ /$randomisedDir\/(.*?)\.(fasta|fa|fna)$/){
	    $randFileRoot = $1;
	}
	print "$randomised[$fileI]\n" if (defined $verbose);
	my $rOrfOutTemp = $randomised[$fileI] . ".orfout.temp";
	unlink($rOrfOutTemp) if (-s $rOrfOutTemp);
	
	#######################################################################
	my $nullSeqStats = $randomised[$fileI] . ".nullseqstats";
	open(NUT, "> $nullSeqStats") or die "FATAL: failed to open [$nullSeqStats] for writing!\n[$!]\n";
	#Print null sequence stats loop...
	my ($rSeq, $rSeqCnt,$rName) = ("go", 1, "name"); 
	#For each random/null sequence... 
	while(length($rSeq)){

	    if(-s $randomised[$fileI]){
		($rName, $rSeq) = parseNthFastaSeq($randomised[$fileI], $rSeqCnt);
	    }
	    else {
		last;
	    }
	    
	    #print "rName:[$rName], rSeq:[$rSeq]\n" if defined $verbose; 	    
	    #print "rSeqCnt:[$rSeqCnt] rSeq:[$rSeq]\n";
	    if (length($rSeq)){
		$rSeqCnt++;		
	    }
	    else{
		last;
	    }
	    
	    my $randomFreqs = &extract0123merFreqs($rSeq);
	    
	    #PRINT NUC STATS....
	    printNucFreqStats($nativeFreqs, $randomFreqs, $rName, \*NUT);
	    
	    #ORF stats: 
	    open(RUT, "> blah.rut");
	    print RUT ">$rName\n$rSeq\n";
	    close(RUT);
	    
	    computeORFstats("blah.rut", $rOrfOutTemp); 
	    
	    #FIX THIS!!! IF NO ORFs THERE IS AN ERROR!!! 
	    #my $rOrfExe = "esl-translate -W --informat fasta -l 100 blah.rut > blah.translate.random && esl-seqstat -c blah.translate.random >> $rOrfOutTemp 2> /dev/null";
	    #system( $rOrfExe ) and printf "WARNING: [$rOrfExe] failed to execute!\n\[$!]\n";
	    unlink("blah.rut") if (not defined $debugQuick or not defined $debug);
	    
	}
	close(NUT);
	#######################################################################
	
	######################################################################
	#Summarise ORF stats 
	my $rOrfOut = $randomised[$fileI] . ".orfout";
	printORFstats($nOrfOut, $rOrfOutTemp, $rOrfOut);
	#unlink($rOrfOutTemp) if (not defined $debugQuick or not defined $debug);
		
	#Sequence similarity -- based on 9-mers:  
	my $mashOut = $randomised[$fileI] . ".mashout";
	my $mashExe = "mash dist -w 1.0 -i -n -k 9 $randomised[$fileI] $seqFile  > $mashOut 2> /dev/null";
	system( $mashExe ) and printf "WARNING: [$mashExe] failed to execute!\n\[$!]\n";
	#output fields are [reference-ID, query-ID, distance, p-value, shared-hashes]
	
	######################################################################
	#Modify to compute Z-scores: 
	#system("cat $nullSeqStats  $mashOut $rOrfOut ") if defined($verbose);
	my $zFile = $randomised[$fileI] . ".zscores";
	computeZScores($nullSeqStats,  $mashOut, $rOrfOut, $zFile);
	
	printf "$randFileRoot\n" if defined($verbose);
	system("cat $zFile ") if defined($verbose);
	#######################################################################
	last if (defined $debugQuick                  ); #DEBUG MODE
    }
    last     if (defined $debugQuick or defined $debug); #DEBUG MODE
}

#Remove temp files:
#unlink("/tmp/blah.rut", "/tmp/blah");

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

exit(0);


######################################################################
#computeZScores: compute Z-scores for each set of statistic 
sub computeZScores {
    
    my ($nullSeqStats,  $mashOut, $rOrfOut, $outFile) = @_;

    ###################################
    #Sequence statistics (length, mono-, di- & tri-nucleotide frequencies 
    open(SSTAT, "< $nullSeqStats") or die("FATAL: failed to open [$nullSeqStats] for reading\n[$!]");
    my (@nullLen, @nullMonoEnt, @nullDiEnt, @nullTriEnt);
    my ($natLen, $cg);
    while( my $in = <SSTAT>){
	if($in =~ /0:len:(\d+)\/(\d+):\S+\s+1:(\S+)\s+2:(\S+)\s+3:(\S+)\s+cg:(\S+)/){
	    $natLen=$1; 
	    push(@nullLen,     $2); 
	    push(@nullMonoEnt, $3); 
	    push(@nullDiEnt,   $4); 
	    push(@nullTriEnt,  $5);
	    $cg              = $6;
	}
    }
    close(SSTAT); 
    
    my $header   = "fileName\tc+g\tlen\tlen.Z";    
    my $printStr .= $outFile . "\t$cg\t$natLen\t" . computeZscore(\@nullLen, $natLen);
    my @zscores; 
    
    #ADD WEIGHTED SUM OF Z-SCORES
    #ADD LENGTH & G+C!!!
    
    $header   .= "\tmonoNuc.Z";
    push(@zscores, computeZscore(\@nullMonoEnt, 0)); 
    $printStr .= "\t" . $zscores[0];
    $header   .= "\tdiNuc.Z";    
    push(@zscores, computeZscore(\@nullDiEnt,   0)); 
    $printStr .= "\t" . $zscores[1];

    $header   .= "\ttriNuc.Z";    
    push(@zscores, computeZscore(\@nullTriEnt,  0)); 
    $printStr .= "\t" . $zscores[2];
    
    ###################################
    #ORF stats: 
    open(ORFIN, "< $rOrfOut") or die("FATAL: failed to open [$rOrfOut] for reading\n[$!]");
    while( my $in = <ORFIN>){
	if($in =~ /orf.ZnoSeqs:(\S+)\s+orf.ZtotRes:(\S+)\s+orf.ZaveLen:(\S+)/){
	    $header   .= "\torf.noSeqs.Z";    
	    $printStr .= "\t" . $1;
	    push(@zscores,      $1); 
	    
	    $header   .= "\torf.totRes.Z";    
	    $printStr .= "\t" . $2;
	    push(@zscores,      $2); 

	    $header   .= "\torf.aveLen.Z";    
	    $printStr .= "\t" . $3;
	    push(@zscores,      $3); 	    
	}
    }
    close(ORFIN);
    ###################################
    #Mash -- sequence distance 
    open(MASHIN, "< $mashOut") or die("FATAL: failed to open [$mashOut] for reading\n[$!]");
    my (@mashDist);
    while( my $in = <MASHIN>){
	my @in = split(/\s+/, $in);  
	if(defined( $in[2] ) and isNumeric($in[2]) ){
	    #output fields are [reference-ID, query-ID, distance, p-value, shared-hashes]
	    $natLen=$1; 
	    push(@mashDist,     $in[2]);
	}
    }
    close(MASHIN); 
    
    $header   .= "\tmash.dist.Z";    
    push(@zscores, computeZscore(\@mashDist,  0)); 
    $printStr .= "\t" . $zscores[6];
    
#    $printStr .= "\t" . computeZscore(\@mashDist, 0);    
    
    $header   .= "\tweighted.Z";
    #                      Nucleotide Composition Scores          |  ORF Scores                             | Sequence Similarity Scores
    $printStr .= "\t" . ( ($zscores[0]+$zscores[1]+$zscores[2])/3 - ($zscores[3]+$zscores[4]+$zscores[5])/3 + $zscores[6]  ) ;
    
    
    open( OUTFILE, "> $outFile") or die("FATAL: failed to open [$outFile] for reading\n[$!]");
    print OUTFILE "$header\n$printStr\n"; 
    close(OUTFILE);
    
    return 0; 
}

######################################################################
#computeZscore: compute a Z-score from an array & x-value 
sub computeZscore {
    my $input = shift; #reference to an array of numbers
    my $x     = shift; #single value
    my $mean=Math::NumberCruncher::Mean(             $input);
    my $stdv=Math::NumberCruncher::StandardDeviation($input);
    if(defined($stdv) && $stdv > 0.00001){
	return sprintf( "%0.5f", ($mean-$x)/$stdv); 
    }
    else {
	return '0.00000';
    }
}


######################################################################
#computeORFstats: translate ORFs and compute statistics for these -- NB. appends output file:
sub computeORFstats {
    
    my ($inFile, $outFile) = @_;
    my $nOrfExe1 = "esl-translate -W --informat fasta -l 100 $inFile > blah.translate.native"; 
    system( $nOrfExe1 ) and printf "WARNING: [$nOrfExe1] failed to execute!\n\[$!]\n";
    
    if(-s "blah.translate.native" ){
	my $nOrfExe2 = "esl-seqstat -c blah.translate.native >> $outFile 2> /dev/null";
	system( $nOrfExe2 ) and printf "WARNING: [$nOrfExe2] failed to execute!\n\[$!]\n";
    }
    else {
	open(ORFOUT, ">> $outFile") or die("FATAL: failed to open [$outFile] for writing\n[$!]");
	print ORFOUT "Format:              FASTA
Alphabet type:       amino
Number of sequences: 0
Total # residues:    0
Smallest:            0
Largest:             0
Average length:      0

Residue composition:
";	
	my @aa = qw(A C D E F G H I K L M N P Q R S T V W Y);
	foreach my $aa (@aa){
	    print ORFOUT "residue: $aa        0  0.0000    0.0000\n";
	}	
	close(ORFOUT); 
    }
    
    unlink("blah.translate") if (not defined $debugQuick or not defined $debug);
    return 0;
    
    
}


######################################################################
#printORFstats: print the ave length and amino-acid frequency statistics:
sub printORFstats {
    my ($nOrfOut, $rOrfOutTemp, $rOrfOut) = @_;

    my (%nativeStats, %randomStats);
    %nativeStats = %{ parseEslSeqStat($nOrfOut    ) } if( -s $nOrfOut    );
    %randomStats = %{ parseEslSeqStat($rOrfOutTemp) } if( -s $rOrfOutTemp);
    
    open(ORFUT, "> $rOrfOut") or die("FATAL: failed to open [$rOrfOut] for writing\n[$!]");
    
    if( defined($nativeStats{'noSeqs'}[0]) && defined($randomStats{'noSeqs'}[0]) && scalar( @{ $randomStats{'noSeqs'} }  )){
	printf ORFUT "orf.ZnoSeqs:%0.5f\t ", computeZscore(\@{ $randomStats{'noSeqs'}}, $nativeStats{'noSeqs'}[0]);
    }
    
    if( defined($nativeStats{'totRes'}[0]) && defined($randomStats{'totRes'}[0]) && scalar( @{ $randomStats{'totRes'} }  )){
	printf ORFUT "orf.ZtotRes:%0.5f\t ", computeZscore(\@{ $randomStats{'totRes'}}, $nativeStats{'totRes'}[0]);
    }
    
    if( defined($nativeStats{'aveLen'}[0]) && defined($randomStats{'aveLen'}[0]) && scalar( @{ $randomStats{'aveLen'} }  )){
	printf ORFUT "orf.ZaveLen:%0.5f\t ",  computeZscore(\@{ $randomStats{'aveLen'}}, $nativeStats{'aveLen'}[0]);
    }
    
    # my @aa = qw(A C D E F G H I K L M N P Q R S T V W Y);
    # my $entropy=0;
    # foreach my $aa (@aa){
    # 	next if (not defined($nativeStats{'res' . $aa}[0] ) or not defined($randomStats{'res' . $aa}[0]));
    # 	next if (not isNumeric($nativeStats{'res' . $aa}[0]) or not isNumeric($randomStats{'res' . $aa}[0]) or
    # 		 $nativeStats{'res' . $aa}[0] == 0 or $randomStats{'res' . $aa}[0] == 0);
	
    # 	my $meanResFreq=Math::NumberCruncher::Mean(             \@{ $randomStats{'res' . $aa}   });
    # 	my $stdvResFreq=Math::NumberCruncher::StandardDeviation(\@{ $randomStats{'res' . $aa}   });
	
    # 	$entropy+= $nativeStats{'res' . $aa}[0]*log( $nativeStats{'res' . $aa}[0] / $meanResFreq )/log(2);
    # 	printf ORFUT "z.$aa:%0.5f\t ", ($meanResFreq-$nativeStats{'res' . $aa}[0])/ $stdvResFreq;	
    # }
    # printf ORFUT "AARelEnTropy:%0.5f\t", $entropy; 	
    
    printf ORFUT "\n";
    close(ORFUT); 
    
    return 0; 
}

######################################################################
#parseEslSeqStat: parse the statistics from esl-seqstat output from multi-protein fasta files.
#                 Could be multiple concatenated seqstat outputs -- returns a hash of arrays:
sub parseEslSeqStat {
    
    my $inFile = shift; 
    my %stats;
    
    open(STATIN, "< $inFile") or die("FATAL: failed to open [$inFile] for reading\n[$!]");
    
    while(my $in = <STATIN>){
	
	if($in =~ /Number of sequences:\s+(\d+)/){
	    push( @{ $stats{'noSeqs'} }, $1 ); 
	}
	elsif($in =~ /Total # residues:\s+(\d+)/){
	    push( @{ $stats{'totRes'} }, $1 ); 
	}
	elsif($in =~ /Average length:\s+(\S+)/){
	    push( @{ $stats{'aveLen'} }, $1 ); 
	}
	elsif($in =~ /residue: (\w)\s+\d+\s+(\S+)/){
	    push( @{ $stats{'res' . $1} }, $2 ); 
	}
    }
    close(STATIN); 
    
    return \%stats;
    
}

######################################################################
#printNucFreqStats: print the length, mono-, di- & tri- frequency statistics:
sub printNucFreqStats {
    my ($nativeFreqs, $randomFreqs, $rName, $nut) = @_;
    
    my %nativeFreqs = %{$nativeFreqs};
    my %randomFreqs = %{$randomFreqs};
    
    my $absDiffRat    = abs($nativeFreqs{0}{'len'} - $randomFreqs{0}{'len'})/$nativeFreqs{0}{'len'} if (defined($nativeFreqs{0}{'len'}));
    my $lenRelEntropy = 0.0;
    $lenRelEntropy = $absDiffRat * log( $absDiffRat )/log(2) if($absDiffRat>0);
    printf { $nut } "$rName\t0:len:%d/%d:%0.5f\t", $nativeFreqs{0}{'len'}, $randomFreqs{0}{'len'}, $lenRelEntropy;  
    for (my $kmerLen=1; $kmerLen<=3; $kmerLen++){
	my $entropy=0;
	printf { $nut } "$kmerLen"; 
	foreach my $kmer (sort keys %{$nativeFreqs{$kmerLen}}) {
	    #printf "$kmerLen\t$kmer\t%0.4f\t%0.4f\n", $nativeFreqs{$kmerLen}{$kmer}, $randomFreqs{$kmerLen}{$kmer};
	    next if (not defined($nativeFreqs{$kmerLen}{$kmer}) or not defined($randomFreqs{$kmerLen}{$kmer}));
	    #printf ":$kmer:%0.5f/%0.5f", $nativeFreqs{$kmerLen}{$kmer}, $randomFreqs{$kmerLen}{$kmer};
	    next if (not isNumeric($nativeFreqs{$kmerLen}{$kmer}) or not isNumeric($randomFreqs{$kmerLen}{$kmer}) or
		     $nativeFreqs{$kmerLen}{$kmer} == 0 or $randomFreqs{$kmerLen}{$kmer} == 0);
	    $entropy+= $nativeFreqs{$kmerLen}{$kmer}*log( $nativeFreqs{$kmerLen}{$kmer} / $randomFreqs{$kmerLen}{$kmer} )/log(2);
	}
	printf { $nut } ":%0.5f\t", $entropy; 
    }
    
    printf { $nut } "cg:%0.5f\n", $nativeFreqs{1}{'C'} + $nativeFreqs{1}{'G'};
    #printf { $nut } "\n";

    return 0; 
}



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
    $seq =~ tr/a-z/A-Z/;#Must be uppercase!
    my $len = length($seq);
    
    my %kmerFreq;
    
    $kmerFreq{0}{'len'} = $len;
    
    
    #printf "len:%d\n", $kmerFreq{0}{'len'}; 
    
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
	    #printf "kmer[$kmer] kmerLen[$kmerLen] %0.2f\n", $kmerFreq{$kmerLen}{$kmer};
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

######################################################################
sub help {
    print STDERR <<EOF;

$0: compare statistics for null (randomised) and native sequences for each set of null sequences using a range of different stats.

Usage:   $0 -i <seqfile> 
  Options:
    -N|--natdir  [str]            Directory containing fasta files for the native sequences 
    -r|--randdir [str]            Directory containing fasta files for the null sequences (filename prefixes should match native sequences) 
    -h|--help                     Show this help.
    -v|--verbose                  Print lots of stuff.

  Dependencies:
  TODO:

EOF
}





