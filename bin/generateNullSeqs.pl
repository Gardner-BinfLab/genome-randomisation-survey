#!/usr/bin/perl 

#generateNullSeqs.pl: Generate null (randomised) sequences for each fasta file in a directory 

use strict;
use warnings;
use Getopt::Long;
use Time::HiRes 'gettimeofday', 'tv_interval';
my ($help, $verbose); 
my ($num,$inDir,$outDir) = (10,"./","./"); 

&GetOptions(
    "n|num=s"            => \$num,
    "d|indir=s"          => \$inDir,
    "o|outdir=s"         => \$outDir,
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

my $timeOut = "$outDir/temp.times.txt"; 
unlink($timeOut);

foreach my $seqFile (@seqFiles){
    my $fileRoot = $$ . "-" . rand(10000); 
    if($seqFile =~ /$inDir\/(.*?)\.(fasta|fa|fna)$/){
	$fileRoot = $1;
    }
    
    print "seqFile: [$seqFile] fileRoot: [$fileRoot]\n" if defined($verbose);
    
    #Some methods require sequence info: 
    my $seq = parseFastaSeq($seqFile);     
    #my $seqLength = length($seq);
    
    my ($aFreq, $cFreq, $gFreq, $tFreq, $seqLength) = seqStat($seq);     
    printf "len:%d A:%0.1f C:%0.1f G:%0.1f T:%0.1f\n", $seqLength, $aFreq/$seqLength, $cFreq/$seqLength, $gFreq/$seqLength, $tFreq/$seqLength if defined($verbose);
    
    ######################################################################
    #EASEL ESL-SHUFFLE OPTIONS: explore shuffling in a local window???
    my $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.mononuc.fasta"; 
    runEslShuffle($seqFile, $num, $outFile, "-m", $timeOut)   if (not -s $outFile);
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.dinuc.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-d", $timeOut)   if (not -s $outFile); 
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.kmer6.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-k 6", $timeOut) if (not -s $outFile); 
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.markov0.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-0", $timeOut)   if (not -s $outFile); 
    $outFile = $outDir . "/" . $fileRoot . ".eslshuffle.markov1.fasta" ;
    runEslShuffle($seqFile, $num, $outFile, "-1", $timeOut)   if (not -s $outFile); 

    ######################################################################
    #EMBOSS SHUFFLESEQ
    $outFile = $outDir . "/" . $fileRoot . ".shuffleseq.fasta" ;
    runShuffleSeq($seqFile, $num, $outFile, $timeOut) if (not -s $outFile);     

    ######################################################################
    #USHUFFLE
#    $outFile = $outDir . "/" . $fileRoot . ".ushuffle-k1.fasta" ;
#    runUshuffle($seq, $num, $outFile, "-k 1", $timeOut);
#    $outFile = $outDir . "/" . $fileRoot . ".ushuffle-k2.fasta" ;
#    runUshuffle($seq, $num, $outFile, "-k 2", $timeOut);
#    $outFile = $outDir . "/" . $fileRoot . ".ushuffle-k6.fasta" ;
#    runUshuffle($seq, $num, $outFile, "-k 6", $timeOut);
#ALSO: ######################################################################
    #https://meme-suite.org/meme/doc/fasta-shuffle-letters.html
    #"The underlying implementation uses uShuffle."
    
    ######################################################################
    #SHUFFLER
    #3 main modes: 1. Eulerian (DEFAULT preserve k-mer counts with a random Eulerian walk through k-mer graph),
    #              2. Linear (split sequence into k-mers and shuffle these)
    #              3. Markov (emit k-mers based upon frequency in input sequence) 

    #Eulerian:
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.eulerian1.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 1",    $timeOut) if (not -s $outFile);    
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.eulerian2.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 2",    $timeOut) if (not -s $outFile);    
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.eulerian6.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 6",    $timeOut) if (not -s $outFile);    
    #Linear:
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.linear1.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 1 -l", $timeOut) if (not -s $outFile);    
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.linear2.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 2 -l", $timeOut) if (not -s $outFile);    
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.linear6.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 6 -l", $timeOut) if (not -s $outFile);    
    #Markov:
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.markov1.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 1 -m", $timeOut) if (not -s $outFile);    
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.markov2.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 2 -m", $timeOut) if (not -s $outFile);    
    $outFile = $outDir . "/" . $fileRoot . ".shuffler.markov6.fasta";
    runShuffler($seqFile, $num, $outFile, "-k 6 -m", $timeOut) if (not -s $outFile);    

    ######################################################################
    #BACTOME RANDOMSEQ
    #python randomseq.py FLS --length=100 --n=10 --allow_start=False --allow_stop=False --start_codons="TTG,CTG,ATG" --stop_codons="TAA,TAG,TGA" --cap_start=True --cap_stop=True --selection="A,250;T,250;G,250;C,250" --fasta=True --prefix="Test"
    
    $outFile = "$outDir/" . $fileRoot . ".randomseq.unconstrained.fasta";
    runRandomSeq($outFile, "--allow_start=True --allow_stop=True --cap_start=False --cap_stop=False --length=$seqLength --n=$num --prefix=\42$fileRoot\42 --selection=\42A,$aFreq;T,$tFreq;G,$gFreq;C,$cFreq\42 ", $timeOut) if (not -s $outFile);
    
    #TOO SLOW -- POSSIBLY IMPOSSIBLE TO SATISFY CONSTRAINTS IN LONG SEQS...
#    $outFile = "$outDir/" . $fileRoot . ".randomseq.constrained";
#    runRandomSeq($outFile, "--allow_start=False --allow_stop=False --cap_start=False --cap_stop=False --length=$seqLength --n=$num --prefix=\42$fileRoot\42 --selection=\42A,$aFreq;T,$tFreq;G,$gFreq;C,$cFreq\42 ", $timeOut);
    
    ######################################################################
    #simuG: randomly add SNPs/INDELs to a sequence
    #Make snp_count proportional to length, not using indel option as the code hangs.
    
    #requests: a "-n" option to generate $num sequences...  
    
    #This code may be too buggy! The method fails on some runs.   :-/
    #One percent:
    $outFile = $outDir . "/" . $fileRoot . ".simuG.01percent.fasta";
    my $options = sprintf("-snp_count %d ", int(0.01*$seqLength) ); 
    runSimuG($seqFile, $num, $outFile, $options, $timeOut) if (not -s $outFile);
    #Ten percent:
    #$outFile = $outDir . "/" . $fileRoot . ".simuG.10percent.fasta";
    #$options = sprintf("-snp_count %d ", int(0.1*$seqLength) ); 
    #runSimuG($seqFile, $num, $outFile, $options, $timeOut);
    #Fifty percent:
    #$outFile = $outDir . "/" . $fileRoot . ".simuG.50percent.fasta";
    #$options = sprintf("-snp_count %d ", int(0.5*$seqLength) ); 
    #runSimuG($seqFile, $num, $outFile, $options, $timeOut);

    #Not able to run "-coding_partition_for_snp_simulation coding -gene_gff seq.gff3"
    #Appears to require a highly specific GFF format...
    #simuG.pl -coding_partition_for_snp_simulation coding -gene_gff U00096-short.3.gff -refseq U00096-short.3.fasta -prefix sequence.simuG -snp_count 50
    
    
    ######################################################################
    #ROSE
    $outFile = $outDir . "/" . $fileRoot . ".rose-001.fasta" ;
    runRose($seq, $num, $outFile,  0.1, $timeOut) if (not -s $outFile);
    $outFile = $outDir . "/" . $fileRoot . ".rose-010.fasta" ;
    runRose($seq, $num, $outFile,  1.0, $timeOut) if (not -s $outFile);
    $outFile = $outDir . "/" . $fileRoot . ".rose-100.fasta" ;
    runRose($seq, $num, $outFile, 10.0, $timeOut) if (not -s $outFile);
    
    ######################################################################
    #PhastSim
    $outFile = $fileRoot . ".phastsim-001" ;
    runPhastSim($seqFile, $num, "$outDir/", $outFile,  0.1, "",        $timeOut) if (not -s "$outDir/$outFile\.fasta");
    $outFile = $fileRoot . ".phastsim-010" ;
    runPhastSim($seqFile, $num, "$outDir/", $outFile,  1.0, "",        $timeOut) if (not -s "$outDir/$outFile\.fasta");

    $outFile = $fileRoot . ".phastsim-codon-001" ;
    runPhastSim($seqFile, $num, "$outDir/", $outFile,  0.1, "--codon", $timeOut) if (not -s "$outDir/$outFile\.fasta");
    $outFile = $fileRoot . ".phastsim-codon-010" ;
    runPhastSim($seqFile, $num, "$outDir/", $outFile,  1.0, "--codon", $timeOut) if (not -s "$outDir/$outFile\.fasta");

    
#    last; 
}

my $timeOutNew = "$outDir/times.txt"; 
timeOut2TSV($timeOut,$timeOutNew,$outDir);

#Fisher-Yates check! 
print "CHECK FOR OVER REPRESENTED SEQS/IMBALANCE!!!!!\n";

#Randomnly sample an 8-mer from $seq & 0:$seqLength - 8; 



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
    my $exe = "esl-shuffle -N $n $options -o $outfile $seqfile";
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    my $start = [ gettimeofday() ];
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42  $exe 2>> $timeOut") and die "FATAL: [$exe] failed!\n[$!]\n";
    #Print time spent in system loop:
    my $total_secs = tv_interval($start);
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);

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
    my $exe = "shuffleseq  -shuffle $n -sequence $seqfile -outseq $outfile.temp ";
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    my $start = [ gettimeofday() ];
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and die "FATAL: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]\n";
    #Print time spent in system loop:
    my $total_secs = tv_interval($start);
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);
    
    #convert shuffleseq output format to uniquely named fasta sequences!!!
    open( SHQ,"< $outfile\.temp") or die("FATAL: failed to open [$outfile\.temp]\n[$!]\n");
    open(FSHQ,"> $outfile")       or die("FATAL: failed to open [$outfile]\n[$!]\n");
    my $cnt=0;
    while(my $in = <SHQ>){
	if($in =~ /^>\s*(\S+)$/){
	    $cnt++;
	    print FSHQ ">shuffleseq.$cnt $1\n";	    
	}
	else{
	    print FSHQ "$in";	    
	}
    }
    close( SHQ);
    close(FSHQ);
    unlink("$outfile\.temp");
    print "WARNING: sequence-count[$cnt] not equal to expectation [$n] for shuffleseq [$outfile]\n" if ($cnt != $n);
    
    return 0;
}


######################################################################
#runRandomSeq: run randomseq on an input sequence composition
#Usage: python randomseq.py FLS --length=100 --n=10 --allow_start=False --allow_stop=False --start_codons="TTG,CTG,ATG" --stop_codons="TAA,TAG,TGA" --cap_start=True --cap_stop=True --selection="A,250;T,250;G,250;C,250" --fasta=True --prefix="Test"
sub runRandomSeq {
    my ($outfile, $options, $timeOut)=@_;
    my $exe = "randomseq.py FLS --fasta=True $options > $outfile\.temp ";
    
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    my $start = [ gettimeofday() ];
    my $timedExe = "/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut";
    system($timedExe) and die "FATAL: [$timedExe] failed!\n[$!]\n";
    #Print time spent in system loop:
    my $total_secs = tv_interval($start);
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);
    
    #convert randomseq output format to uniquely named fasta sequences!!!
    open( RAS,"< $outfile\.temp") or die("FATAL: failed to open [$outfile\.temp]\n[$!]\n");
    open(FRAS,"> $outfile")       or die("FATAL: failed to open [$outfile]\n[$!]\n");
    my $cnt=0;
    while(my $in = <RAS>){
	if($in =~ /^>\s*(\S+)$/){
	    $cnt++;
	    print FRAS ">randomseq.$cnt $1\n";	    
	}
	else{
	    print FRAS "$in";	    
	}
    }
    close( RAS);
    close(FRAS);
    unlink("$outfile\.temp");
    print "WARNING: sequence-count[$cnt] not equal to expectation [$n] for randomseq [$outfile]\n" if ($cnt != $n);
    
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
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile >> $timeOut");
    my $start = [ gettimeofday() ];
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and die "FATAL: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]\n";
    #Print time spent in system loop:
    my $total_secs = tv_interval($start);
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);
    
    #convert ushuffle format to fasta!!!
    open( USH,"< $outfile\.temp") or die("FATAL: failed to open [$outfile\.temp]\n[$!]\n");
    open(FUSH,"> $outfile")       or die("FATAL: failed to open [$outfile]\n[$!]\n");
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
#runShuffler: run shuffler on an input sequence file
# shuffler v1.3  Copyright (C) 2019  Benjamin Jean-Marie Tremblay                 
# Usage:  shuffler [options] -i [filename] -o [filename]                                                                                                          
#  -i <str>   Input filename. All white space will be removed.                    
#  -o <str>   Output filename. 
#  -k <int>   K-let size. Defaults to 1.                                          
#  -m         Use the markov shuffling method. Defaults to euler.                 
#  -l         Use the linear shuffling method. Defaults to euler.                 
#  -f         Indicate the input is fasta formatted. 
sub runShuffler {
    my ($seqfile, $n, $outfile, $options, $timeOut) = @_;    
    my $tempOut = $outfile . ".temp";
    my $exe = "shuffler -n $n -f $options -i $seqfile -o $tempOut";
    print "Running: [$exe]\n" if defined($verbose);
    
    system("echo $outfile >> $timeOut");
    
    #Execute & record run times: 
    my $start = [ gettimeofday() ];
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and die "FATAL: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]\n";
    my $total_secs = tv_interval($start);

    #Rename sequences:
    open(SIN, "< $tempOut") or die "FATAL: failed to open [$tempOut]!\n[$!]\n";
    open(SUT, "> $outfile") or die("FATAL: failed to open [$outfile]!\n[$!]\n");
    my $i=1;
    while( my $sin = <SIN>){
	if($sin =~ /^>(.*)/){
	    print SUT ">shuffler$i $1\n";
	    $i++;
	}
	else{
	    print SUT $sin;
	}
    }
    close(SIN);
    close(SUT);
    unlink($tempOut);
    
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);
    return 0;
}


######################################################################
#runSimuG: run shuffler on an input sequence file
#simuG.pl -refseq sequence.fasta -snp_count 1000 -indel_count 5 -cnv_count 1 -inversion_count 10 -translocation_count 1  -prefix sequence.simuG
sub runSimuG {
    my ($seqfile, $n, $outfile, $options, $timeOut) = @_;    
    my $exe = "simuG.pl -refseq $seqfile -prefix sequence.simuG $options > /dev/null ";
    print "Running: [$exe]\n" if defined($verbose);
    
    open(UT, "> $outfile") or die("FATAL: failed to open [$outfile]\n[$!]\n");
    system("echo $outfile >> $timeOut");
    
    #PRINT TIME...
    my $start = [ gettimeofday() ];
    for (my $i =0; $i<$n; $i++){
	system($exe) and die "FATAL: failed to execute [$exe]!\n$!";
	open(SIN, "< sequence.simuG.simseq.genome.fa") or die("FATAL: failed to open [sequence.simuG.simseq.genome.fa]\n[$!]\n");
	while( my $sin = <SIN>){
	    if($sin =~ /^>(.*)/){
		print UT ">simuG$i $1\n";
	    }
	    else{
		print UT $sin;
	    }
	}
	close(SIN);
    }
    close(UT);

    #Print time spent in system loop:
    my $total_secs = tv_interval($start);
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);
    
    unlink("sequence.simuG.simseq.genome.fa", "sequence.simuG.refseq2simseq.SNP.vcf", "sequence.simuG.refseq2simseq.map.txt");
    
    return 0;
}


######################################################################
#runRose: run rose on an input sequence file
#Usage: rose-1.3.1/src/rose -v -h -I <include path> <inputfile> | -
sub runRose {
    my ($seq, $n, $outfile, $len, $timeOut) = @_;
    
    #write rose parameter file: 
    open(ROS,"> $outfile\.rosepar") or die("FATAL: failed to open [$outfile\.rosepar]\n[$!]\n");
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
    my $start = [ gettimeofday() ];
    #system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and die "FATAL: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]\n";
    my $runFail=0; 
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and ($runFail=1);
    if($runFail){
	print "WARNING: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]\n";
	return 1;
    }


    #Print time spent in system loop:
    my $total_secs = tv_interval($start);
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);
    
    #convert rose output to fasta:
    open(ROSIN,"< $outfile\.temp") or die("FATAL: failed to open [$outfile\.temp]\n[$!]\n");
    open(ROSUT,"> $outfile"      ) or die("FATAL: failed to open [$outfile]\n[$!]\n");
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
#runPhastSim: run phastSim on an input sequence file
#Usage: phastSim [-h] [--outpath OUTPATH] [--reference REFERENCE]
#                 [--rootGenomeLength ROOTGENOMELENGTH]
#                 [--rootGenomeFrequencies ROOTGENOMEFREQUENCIES [ROOTGENOMEFREQUENCIES ...]]
#                 [--treeFile TREEFILE] [--scale SCALE] [--seed SEED]
#                 [--alpha ALPHA] [--invariable INVARIABLE]
#                 [--mutationRates MUTATIONRATES [MUTATIONRATES ...]] [--codon]
#                 [--omegaAlpha OMEGAALPHA] [--omegaBeta OMEGABETA]
#                 [--omegaCategoryProbs OMEGACATEGORYPROBS [OMEGACATEGORYPROBS ...]]
#                 [--omegaCategoryRates OMEGACATEGORYRATES [OMEGACATEGORYRATES ...]]
#                 [--categoryProbs CATEGORYPROBS [CATEGORYPROBS ...]]
#                 [--categoryRates CATEGORYRATES [CATEGORYRATES ...]]
#                 [--hyperMutProbs HYPERMUTPROBS [HYPERMUTPROBS ...]]
#                 [--hyperMutRates HYPERMUTRATES [HYPERMUTRATES ...]]
#                 [--noHierarchy] [--noNormalization] [--verbose]
#                 [--outputFile OUTPUTFILE] [--alternativeOutputFormat]
#                 [--createNewick] [--createFasta] [--createPhylip]
#                 [--createInfo] [--createMAT] [--indels]
#                 [--insertionRate INSERTIONRATE [INSERTIONRATE ...]]
#                 [--deletionRate DELETIONRATE [DELETIONRATE ...]]
#                 [--insertionLength INSERTIONLENGTH [INSERTIONLENGTH ...]]
#                 [--deletionLength DELETIONLENGTH [DELETIONLENGTH ...]]
#                 [--mutationsTSVinput MUTATIONSTSVINPUT]

# Efficiently simulate sequence evolution along phylogenies with short branches.

sub runPhastSim {
    my ($seq, $n, $outpath, $outfile, $len, $options, $timeOut) = @_;
    
    #write phastSim tree file:
    my $treeFile = $outpath . $outfile . '.phastsimtree';
    open(PHAS,"> $treeFile") or die("FATAL: failed to open [$treeFile]\n[$!]\n");
    printf PHAS "((";
    for (my $i = 0; $i<$n-1; $i++){
	print PHAS "phastsim$i:$len,";
    }
    print PHAS "phastsim$n:$len):0.001);";
    close(PHAS);
    
    #run phastsim:
    my $exe = "phastSim --outpath $outpath --outputFile $outfile\.temp --reference $seq --treeFile $treeFile  --createFasta $options --hyperMutProbs 0.01 0.04 --hyperMutRates 100 10 --indels --insertionRate GAMMA 0.1 1.0 --deletionRate CONSTANT 0.1 --insertionLength GEOMETRIC 0.9 --deletionLength NEGBINOMIAL 2 0.95  > /dev/null";
    print "Running: [$exe]\n" if defined($verbose);
    system("echo $outfile\.fasta >> $timeOut");
    my $start = [ gettimeofday() ];
    my $runFail=0; 
    system("/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut") and ($runFail=1);
    if($runFail){
	print "WARNING: [/usr/bin/time -f \42real \%e\\nuser \%U\\nsys \%S\\ncpu \%P\\n\42 $exe 2>> $timeOut] failed!\n[$!]\n";
	return 1;
    }

    
    #Print time spent in system loop:
    my $total_secs = tv_interval($start);
    open(TO, ">> $timeOut") or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    print TO "Perl time: [$total_secs] (secs)\n\n";
    close(TO);

    #Remove gap symbols:
    open(PHIN, "< $outpath$outfile\.temp.fasta") or die("FATAL: failed to open [$outpath$outfile\.temp.fasta]\n[$!]\n");
    open(PHUT, "> $outpath$outfile\.fasta")      or die("FATAL: failed to open [$outpath$outfile\.fasta]     \n[$!]\n");
    while(my $in = <PHIN>){
	if($in =~ /^>/){
	    print PHUT $in;
	}
	else{
	    $in =~ s/-//g;
	    print PHUT $in;
	}
    }
    close(PHUT);
    close(PHIN);
    
    unlink(("$outpath$outfile\.temp.fasta", "$outpath$outfile\.temp.txt"));    
    unlink("$treeFile");
    return 0;
}
    
######################################################################
#parseFastaSeq: return first biosequence in a given fasta file: 
sub parseFastaSeq {
    my $inFile = shift;
    open(IN, "< $inFile") or die("FATAL: failed to open [$inFile]\n[$!]\n");
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
#timeOut2TSV: convert timeOut to a tsv format:
sub timeOut2TSV {
    my ($timeOut,$timeOutNew,$outDir) = @_;
    
    open(TIN, "< $timeOut"   ) or die("FATAL: failed to open [$timeOut]\n[$!]\n");
    open(TUT, "> $timeOutNew") or die("FATAL: failed to open [$timeOutNew]\n[$!]\n");
    print TUT "ID\treal\tuser\tsys\tcpu\tperlTime\n";
    my ($id,$real,$user,$sys,$cpu,$perlTime);
    while(my $tin = <TIN>){
	chomp($tin);
	if($tin =~ /.*fasta$/){
	    $tin =~ s/$outDir\///;
	    $tin =~ s/\.fasta$//;
	    $id = $tin;
	}
	elsif($tin =~ /^real\s+(\S+)/){
	    $real=$1;
	}
	elsif($tin =~ /^user\s+(\S+)/){
	    $user=$1;
	}
	elsif($tin =~ /^sys\s+(\S+)/){
	    $sys=$1;
	}
	elsif($tin =~ /^cpu\s+(\S+)\%/){
	    $cpu=$1;
	}
	elsif($tin =~ /^Perl time: \[(\S+)\] \(secs\)/){
	    $perlTime=$1;
	    print TUT "$id\t$real\t$user\t$sys\t$cpu\t$perlTime\n";
	    ($id,$real,$user,$sys,$cpu,$perlTime)=('blah','NA','NA','NA','NA','NA','NA');
	}
    }
    close(TUT);
    close(TIN);
    
    #unlink($timeOut);
    
    return 0;
}

######################################################################
#seqStat: return # of A, C, G & T in a 
sub seqStat {
    my $seq = shift;
    $seq =~ tr/a-z/A-Z/;
    my @seq = split(//,$seq);
    my $seqLength=0;
    my %counter;
    $counter{'A'}=$counter{'C'}=$counter{'G'}=$counter{'T'}=0;
    foreach my $nuc (@seq){
	if(isNucleotide($nuc)){
	    $counter{$nuc}++;
	    $seqLength++; 
	}
    }
    
    return ($counter{'A'}, $counter{'C'}, $counter{'G'}, $counter{'T'}, $seqLength); 
}

######################################################################
#isNucleotide: returns true if input character is a nucleotide (IUPAC codes):
sub isNucleotide {
    my $a = shift;
    
    if (defined($a)){
        $a =~ tr/a-z/A-Z/;
    }
    
    if (defined($a) && length($a) && ($a =~ /[ACGT]/) ){   #($a =~ /[ACGUTRYWSMKBDHVN]/) ){ STICK TO UNAMBIGUOUS NUCS FOR NOW...
        return 1;
    }
    else {
        return 0;
    }
    
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
    ushuffle 
    shuffler
    simuG.pl
    rose
    phastSim
    
    System tools: 
    /usr/bin/time
    echo
  TODO:

EOF
}
