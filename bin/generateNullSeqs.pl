#!/usr/bin/perl 

#generateNullSeqs.pl: Generate null (randomised) sequences for each fasta file in a directory 



use strict;
use warnings;
use Getopt::Long;

my ($help, $verbose); 
my ($num,$inDir,$outDir) = (10,"./","./"); 

&GetOptions(
    "n|num=s"              => \$num,
    "d|indir=s"            => \$inDir,
    "o|outdir=s"           => \$outDir,
    "h|help"             => \$help,
    "v|verbose"          => \$verbose
    );

if( $help ) {
    &help();
    exit(1);
}

#read sequence file names in:
my @seqFiles = glob("$inDir/*.fasta  $inDir/*.fa  $inDir/*.fna"); 
print "seqFiles: [@seqFiles]\n" if defined($verbose);

my $timeOut = "$outDir/times.txt"; 
unlink($timeOut);

foreach my $seqFile (@seqFiles){
    my $fileRoot = $$ . "-" . rand(10000); 
    if($seqFile =~ /$inDir\/(.*?)\.(fasta|fa|fna)$/){
	$fileRoot = $1;
    }
    
    print "seqFile: [$seqFile] fileRoot: [$fileRoot]\n" if defined($verbose);

    #EASEL ESL-SHUFFLE OPTIONS: explore shuffling in a local window???
    my $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.mononuc.fasta"; 
    runEslShuffle($seqFile, $num, $outFile, "-m", $timeOut);
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.dinuc.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-d", $timeOut); 
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.kmer6.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-k 6", $timeOut); 
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.markov0.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-0", $timeOut); 
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.markov1.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-1", $timeOut); 

    ######################################################################
    #EMBOSS SHUFFLESEQ
    $outFile = $outDir . "/" . $fileRoot . ".shuffleseq.fasta" ;
    runShuffleSeq($seqFile, $num, $outFile, $timeOut);     

    ######################################################################
    #USHUFFLE
    my $seq = parseFastaSeq($seqFile);     
    $outFile = $outDir . "/" . $fileRoot . ".ushuffle-k1.fasta" ;
    runUshuffle($seq, $num, $outFile, "-k 1", $timeOut);
    $outFile = $outDir . "/" . $fileRoot . ".ushuffle-k2.fasta" ;
    runUshuffle($seq, $num, $outFile, "-k 2", $timeOut);
    $outFile = $outDir . "/" . $fileRoot . ".ushuffle-k10.fasta" ;
    runUshuffle($seq, $num, $outFile, "-k 6", $timeOut);

    ######################################################################
    #SHUFFLER
    #3 main modes: 1. Eulerian (DEFAULT preserve k-mer counts with a random Eulerian walk through k-mer graph),
    #              2. Linear (split sequence into k-mers and shuffle these)
    #              3. Markov (emit k-mers based upon frequency in input sequence) 
    #./bin/shuffler -f -k 1    -i seq.fasta
    #./bin/shuffler -f -k 2    -i seq.fasta
    #./bin/shuffler -f -k 6    -i seq.fasta
    #./bin/shuffler -f -k 1 -l -i seq.fasta
    #./bin/shuffler -f -k 2 -l -i seq.fasta
    #./bin/shuffler -f -k 6 -l -i seq.fasta
    #./bin/shuffler -f -k 1 -m -i seq.fasta
    #./bin/shuffler -f -k 2 -m -i seq.fasta
    #./bin/shuffler -f -k 6 -m -i seq.fasta
    
    ######################################################################
    #simuG: randomnly add SNPs/INDELs to a sequence 
    #simuG.pl -refseq ~/projects/genome-randomisation-survey/data/sequences/AB370334.1.fasta -snp_count 1000 -indel_count 10 -cnv_count 2 -inversion_count 1 -translocation_count 1  -prefix output_prefix
    
    
    ######################################################################
    #ROSE
    $outFile = $outDir . "/" . $fileRoot . ".rose-001.fasta" ;
    runRose($seq, $num, $outFile,  0.1, $timeOut);
    $outFile = $outDir . "/" . $fileRoot . ".rose-010.fasta" ;
    runRose($seq, $num, $outFile,  1.0, $timeOut);
    $outFile = $outDir . "/" . $fileRoot . ".rose-100.fasta" ;
    runRose($seq, $num, $outFile, 10.0, $timeOut);
    
    ######################################################################
    #PhastSim

    #No codon model:
    #phastSim --outpath data/sequences-null/ --reference ~/projects/genome-randomisation-survey/data/sequences/AB370334.1.fasta --treeFile bin/rose.tree.1        --scale 0.3333333333 --createPhylip --createFasta --omegaCategoryProbs 0.1 0.2 0.3 0.4 --omegaCategoryRates 1.0 5.0 2.0 2.3 --hyperMutProbs 0.01 0.04 --hyperMutRates 100 10 --mutationRate JC69 --indels --insertionRate GAMMA 0.1 1.0 --deletionRate CONSTANT 0.1 --insertionLength GEOMETRIC 0.9 --deletionLength NEGBINOMIAL 2 0.95 --rootGenomeFrequencies 0
    #Minimal params:
    #phastSim --outpath data/sequences-null/ --reference ~/projects/genome-randomisation-survey/data/sequences/AB370334.1.fasta --treeFile bin/rose.tree.1  --createFasta --hyperMutProbs 0.01 0.04 --hyperMutRates 100 10 --indels --insertionRate GAMMA 0.1 1.0 --deletionRate CONSTANT 0.1 --insertionLength GEOMETRIC 0.9 --deletionLength NEGBINOMIAL 2 0.95

    #Codon model:
    #phastSim --outpath data/sequences-null/ --reference ~/projects/genome-randomisation-survey/data/sequences/AB370334.1.fasta --treeFile bin/rose.tree.1 --codon --scale 0.3333333333 --createPhylip --createFasta --omegaCategoryProbs 0.1 0.2 0.3 0.4 --omegaCategoryRates 1.0 5.0 2.0 2.3 --hyperMutProbs 0.01 0.04 --hyperMutRates 100 10 --mutationRate JC69 --indels --insertionRate GAMMA 0.1 1.0 --deletionRate CONSTANT 0.1 --insertionLength GEOMETRIC 0.9 --deletionLength NEGBINOMIAL 2 0.95 --rootGenomeFrequencies 0
    
    #Fasta written to "sars-cov-2_simulation_output.fasta"?!!!

    ######################################################################
    
    
    
    last; 
}



exit(0); 

######################################################################
#runEslShuffle: run esl-shuffle on an input sequence file
# Usage: esl-shuffle    [options] <seqfile>  (shuffles individual sequences)
#   -N <n> : generate <n> samples (per input seq/msa)  [1]  (n>0)
#   -o <f> : direct output data to file <f>
#   -m     : shuffle preserving monoresidue composition  [default]
#   -d     : shuffle preserving mono- and di-residue composition
#   -k <n> : shuffle nonoverlapping <n>-mers  (n>0)
#   -0     : generate with 0th order Markov properties per input
#   -1     : generate with 1st order Markov properties per input
#   -w <n> : regionally shuffle inputs in window size <n>  (n>0)

sub runEslShuffle {
    my ($seqfile, $n, $outfile, $options, $timeOut) = @_;    
    my $exe = "esl-shuffle -N $n -o $outfile $options $seqfile";
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42  $exe 2>> $timeOut") and die "FATAL: [$exe] failed!\n[$!]";    
    return 0;
}

######################################################################
#runShuffleSeq: run shuffleseq on an input sequence file
#  Standard (Mandatory) qualifiers:
# [-sequence]          seqall     Sequence(s) filename and optional format, or
#                                 reference (input USA)
# [-outseq]            seqoutall  [<sequence>.<format>] Sequence set(s)
#                                 filename and optional format (output USA)
#  -shuffle            integer    [1] Number of shuffles (Any integer value)
sub runShuffleSeq {
    my ($seqfile, $n, $outfile, $timeOut) = @_;    
    my $exe = "shuffleseq  -shuffle $n -sequence $seqfile -outseq $outfile ";
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and die "FATAL: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]";    
    return 0;
}

######################################################################
#runUshuffle: run ushuffle on an input sequence file
# uShuffle: a useful tool for shuffling biological sequences while preserving the k-let counts
# Options:
#   -s <string>     specifies the sequence
#   -n <number>     specifies the number of random sequences to generate
#   -k <number>     specifies the let size
#   -seed <number>  specifies the seed for random number generator
#   -b              benchmark
sub runUshuffle {
    my ($seq, $n, $outfile, $options, $timeOut) = @_;
    my $exe = "ushuffle -n $n $options -s $seq > $outfile\.temp ";
    print "Running: [ushuffle -n $n $options -s <seq> > $outfile\.temp]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and die "FATAL: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]";

    #convert ushuffle format to fasta!!!
    open( USH,"< $outfile\.temp") or die("FATAL: failed to open [$outfile\.temp]\n[$!]");
    open(FUSH,"> $outfile")       or die("FATAL: failed to open [$outfile]\n[$!]");
    my $cnt=0;
    while(my $in = <USH>){
	if($in =~ /^(\S+)$/){
	    $cnt++;
	    print FUSH ">ushuffle.$cnt\n$1\n";	    
	}		
    }
    close( USH);
    close(FUSH);
    unlink("$outfile\.temp");
    print "WARNING: sequence-count[$cnt] not equal to expectation [$n] for ushuffle [$outfile]\n" if ($cnt != $n);
    return 0;
}

######################################################################
#runRose: run rose on an input sequence file
#Usage: rose-1.3.1/src/rose -v -h -I <include path> <inputfile> | -
sub runRose {
    my ($seq, $n, $outfile, $len, $timeOut) = @_;
    
    #write rose parameter file: 
    open(ROS,"> $outfile\.rosepar") or die("FATAL: failed to open [$outfile\.rosepar]\n[$!]");
    printf ROS "SequenceNum = %d
ChooseFromLeaves = True
TheSequence = \42$seq\42
TheTree = ((", $n+1;
    for (my $i = 0; $i<$n-1; $i++){
	print ROS "rose$i:$len,";
    }
	print ROS "rose$n:$len):0.001);
InputType = 4 // DNA
TheAlphabet = \42ACGT\42
TheFreq = [0.25,0.25,0.25,0.25]
TheInsertThreshold = 0.001
TheDeleteThreshold = 0.001
TheInsFunc = [.2,.2,.2,1,1,1,1]
TheDelFunc = [.2,.2,.2,1,1,1,1]
TheDNAmodel = \42F84\42
";
    close(ROS);

    #run rose: 
    my $exe = "rose $outfile\.rosepar > $outfile\.temp ";
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and die "FATAL: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]";
    
    #convert rose output to fasta:
    open(ROSIN,"< $outfile\.temp") or die("FATAL: failed to open [$outfile\.temp]\n[$!]");
    open(ROSUT,"> $outfile"      ) or die("FATAL: failed to open [$outfile]\n[$!]");
    my $printme=0;
    while(my $in = <ROSIN>){
	if($in =~ /^>(rose\d+)/){
	    $printme++;
	}
	elsif($in =~ /^Alignment/){
	    $printme=0;
	    last;
	}
	print ROSUT $in if ( $printme && length($in)); 
    }
    close(ROSUT);
    close(ROSIN);    
    unlink("$outfile\.rosepar");
    unlink("$outfile\.temp");    
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
sub help {
    print STDERR <<EOF;

$0: generate sequences null (randomised) sequences for each of the input sequences using a range of different models.

Usage:   $0 -i <seqfile> 
  Options:
    -n|--num [num]              
    -d|--dir [str]             
    -h|--help                     Show this help.
    -v|--verbose                  Print lots of stuff.

  Dependencies:
    esl-shuffle (easel package https://github.com/EddyRivasLab/easel)
    shuffleseq  (EMBOSS:6.6.0.0 http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/)
TODO:

EOF
}
