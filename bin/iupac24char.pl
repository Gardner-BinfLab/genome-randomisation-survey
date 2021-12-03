#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
my ($help, $verbose); 
my ($inFile) = (""); 

&GetOptions(
    "i|inFile=s"         => \$inFile,
    "h|help"             => \$help,
    "v|verbose"          => \$verbose
    );

if( $help ) {
    &help();
    exit(1);
}
elsif(not -s $inFile){
    print "FATAL: input file required [-i fileName]!\n";
    &help();
    exit(1);
}

stripFastaSeq($inFile);

exit(); 

######################################################################
#stripFastaSeq: print biosequences in a given fasta file, replacing non-ACGT characters: 
sub stripFastaSeq {
    my $inFile = shift;
    open(IN, "< $inFile") or die("FATAL: failed to open [$inFile]\n[$!]\n");
    my $name;
    my $seq='';
    while(my $in = <IN>){
	if($in =~ /^>\s*(\S+)/ && not defined($name)){
	    $name=$1;
	    print $in;
	}
	else {
	    chomp($in);
	    #$seq .= $in;
	    $seq = replaceNonACGT($in);
	    print "$seq\n";
	}
    }
    close(IN);
    
    
    
    return 0; 
}



######################################################################
#replaceNonACGT: work through a sequence, replace each non-ACGT character with a randomly selected nucleotide from IUPAC code 
sub replaceNonACGT {
    my $seq = shift;
    my @seq = split(//, $seq);

    for (my $i=0; $i<scalar(@seq); $i++){
	next if isNucleotide($seq[$i]);
	
	$seq[$i] = randomIUPACreplacement($seq[$i]); 
	
    }
    $seq = join('', @seq);
    
    return $seq; 
}
######################################################################
#randomIUPACreplacement: replace IUPAC redundancy nucs with a randomly selected alternative:
sub randomIUPACreplacement {
    
    my $nuc = shift; 
    $nuc =~ tr/a-z/A-Z/;
    
    my %IUPAC2counts = (
    A => {
	A => 1,
	C => 0,
	G => 0,
	T => 0,
    },
    C => {
	A => 0,
	C => 1,
	G => 0,
	T => 0,
    },
    G => {
	A => 0,
	C => 0,
	G => 1,
	T => 0,
    },
    T => {
	A => 0,
	C => 0,
	G => 0,
	T => 1,
    },
    U => {
	A => 0,
	C => 0,
	G => 0,
	T => 1,
    },
    R => {
	A => 0.5,
	C => 0,
	G => 0.5,
	T => 0,
    },
    Y => {
	A => 0,
	C => 0.5,
	G => 0,
	T => 0.5,
    },
    S => {
	A => 0,
	C => 0.5,
	G => 0.5,
	T => 0,
    },
    W => {
	A => 0.5,
	C => 0,
	G => 0,
	T => 0.5,
    },
    M => {
	A => 0.5,
	C => 0.5,
	G => 0,
	T => 0,
    },
    K => {
	A => 0,
	C => 0,
	G => 0.5,
	T => 0.5,
    },
    B => {
	A => 0,
	C => 1/3,
	G => 1/3,
	T => 1/3,
    },
    D => {
	A => 1/3,
	C => 0,
	G => 1/3,
	T => 1/3,
    },
    H => {
	A => 1/3,
	C => 1/3,
	G => 0,
	T => 1/3,
    },
    V => {
	A => 1/3,
	C => 1/3,
	G => 1/3,
	T => 0,
    },
    N => {
	A => 0.25,
	C => 0.25,
	G => 0.25,
	T => 0.25,
    },
    );
    
    my $randNum = rand(1);
    my $cumSum  = 0.0; 
    if(not defined($IUPAC2counts{$nuc})){
	print "WARNING: [$nuc] is a non-IUPAC character!\n";
	return $nuc; 	    
    }
    else {
	foreach my $acgt (keys %{$IUPAC2counts{$nuc}}) {
	    return $acgt if ($cumSum <= $randNum && $randNum <= $cumSum+$IUPAC2counts{$nuc}{$acgt});
	    $cumSum+=$IUPAC2counts{$nuc}{$acgt};
	}
    }
    
    print "WARNING: something went wrong, returning [$nuc]! cumSum=[$cumSum], randNum=[$randNum]\n";
    return $nuc;
    
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

$0: take a fasta file as input, replace non-ACGT IUPAC characters with randomly selected alternatives.

Usage:   $0 -i <seqfile> 
  Options:
    -i|--inFile [fastaFile]       Input filename 
    -h|--help                     Show this help.
    -v|--verbose                  Print lots of stuff.

  Dependencies:
    
  TODO:

EOF
}


