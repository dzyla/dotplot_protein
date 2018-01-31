# dotplot_protein
Script for creating a dotplot

Requires matplotlib, numpy, biopython, pandas and tqdm.
It scores the sequences with blosum62 matrix and shows the dotplot above some threshold. Unfortunately, it is not the fastest script and for long sequences (~2000 aa) it takes even 30 min with window size of 20 aa. Small sequences work quite well (1 - 10 s) dependingon the word size.

How to use the script: 

sequence_dotplot.py -f fasta.file -s1 1 -s2 0 -w 10 -t 50

-f fasta file format with sequences
-s1 / -s2 sequence number in the file (if the same sequence use the same number)
-w windo size of calulated aligement
-t treshold of measured values in % of 100 corresponding aminoacids in sequence

![alt text](https://github.com/dzyla/dotplot_protein/blob/master/dotplot.png
)
