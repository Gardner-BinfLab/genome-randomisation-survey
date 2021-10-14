# genome-randomisation-survey

* A survey of different genome randomisation methods.


## Files:

* README.md: this readme file

* LICENSE: software licence 

* bin: directory of software scripts

* data: data directory

* docs: documents, including draft manuscript

## Work flow:

1. Fetch sequences for testing different tools (these cover a range of lengths, C+G contents, "complexities" (genes, repeats, etc).

-- All sequences are sourced from the EMBL nucleotide archive, selected to 

Rfam (14.5) & RNAcentral ()


2. Run & time different randomisation methods

* shuffle methods

 * EASEL  esl-shuffle -m|d (-w)

 * EMBOSS shuffleseq 

 * USHUFFLE

* Markov methods & simulated evolution

 * EASEL esl-shuffle -0|1 (-w) 

 * SQUID shuffle -0|1 (-w) 

 * ROSE

 * EXCLUDED: INDELible (2009 MAC/Windows binaries...)

 * EXCLUDED: SISSIz  -- alignment input

* Monte-Carlo methods

 * D-Tailor, BiDAS

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






 






