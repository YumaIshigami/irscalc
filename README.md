# irscalc

This repository contains analysis scripts for:<br>
Ishigami et al. (2021) A single m<sup>6</sup>A modification in U6 snRNA diversifies exon sequence at the 5’ splice site. (Nature Communications)<br>

Environment:<br>
- python 3.8.5 with packages pysam, Biopython and numpy installed<br>

Script gtfbamcounts.py counts the reads mapped on splice junctions and classify them depending on their splicing pattern.<br>
- Usage: python3 gtfbamcount.py [BAM file name] [GTF file name] [output bamcount file name]<br>

Script irsseq.py calculates the IRS/PCS scores of each intron based on the count file generated by gtfbamcounts.py and exports them with their sequence based on genomic DNA sequence FASTA file.<br>
- Usage: python3 gtfbamcount.py [bamcount file list: control] [bamcount file list: treated] [genome FASTA file name] [output file name]<br>
- Bamcount filenames are delimited by commas without spaces.<br>
- Example: python3 gtfbamcount.py wt1.bam.txt,wt2.bam.txt,wt3.bam.txt ko1.bam.txt,ko2.bam.txt,ko3.bam.txt genome.fa outfile.txt<br>

```python3 ```
