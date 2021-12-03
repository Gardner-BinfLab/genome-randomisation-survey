# genome-randomisation-survey

* A survey of different genome randomisation methods.


## Files:

* README.md: this readme file

* LICENSE: software licence 

* bin: directory of software scripts

* data: data directory

* docs: documents, including draft manuscript

## Work flow:

1. Fetch sequences for testing different tools. Selected sequences
cover a range of lengths and C+G contents. These include single gene
sequences (both ncRNA and protein-coding), viroid, viral, bacterial,
archaeal and eukaryotic chromosomes and genomes.

-- All sequences are sourced from the EMBL nucleotide archive (ENA),
   selected to uniformally cover a range of lengths and C+G contents
   (docs/figures/seqStats.pdf)

-- See data/sequenceList.tsv for a list of sequence accessions and
   descriptions

-- To fetch the sequences from ENA, run:

```
./bin/fetchSeqs.pl
```

#Attempted to use for simuG to randomise sequence, accounting for coding regions -- fails to read GFF:  

#embl2gff.pl    ./data/sequences/AB000109.1.embl >  ./data/sequences/AB000109.1.gff

#Check and strip off non-ACGT characters -- causes phastSim to fail...
ls ./data/sequences/*fasta | awk '{print "echo "$1" && esl-seqstat -c -a "$1}' | sh > ./data/sequences/esl-seqstat.out
egrep  'fasta|^residue|^Total'  ./data/sequences/esl-seqstat.out | egrep -v 'residue: A|residue: C|residue: G|residue: T' | grep -B 2 ^resi
egrep  'fasta|^residue|^Total'  ./data/sequences/esl-seqstat.out | egrep -v 'residue: A|residue: C|residue: G|residue: T' | grep -B 2 ^resi | grep fasta$ | awk '{print "./bin/iupac24char.pl -i "$1"> blah && mv blah "$1}' | sh


2. Run & time different randomisation methods. Selected methods include:

* shuffle methods

 * EASEL  esl-shuffle 

 * EMBOSS shuffleseq 

 * Shuffler

 * USHUFFLE
 

* Markov methods & simulated evolutionary methods

 * EASEL esl-shuffle -0|1 (-w)

 * Shuffler (-m) 

 * ROSE

 * PhastSim

* Add random mutations or sequencing errors (SNPs, INDELs, structural variation etc.)

 * simuG

* Monte-Carlo methods

 * ?????????

* Quick run:


```
./bin/generateNullSeqs.pl -n 10 -d ./data/sequences  -o ./data/sequences-null -v  && cat ./data/sequences-null/times.txt
```


3. Evaluate shuffled sequences

 * similarity with input sequence

 * C+G content and shared k-mers

 * # ORFs, ncRNAs (Rfam), domains (Pfam), ...



 Combine statistics:
https://www.nature.com/articles/s41598-021-86465-y
 	 Fisher's method:
	 -2*\sum log(p)  ~ \chi^2


4. Generate figure(s)

5. Write manuscript






 






