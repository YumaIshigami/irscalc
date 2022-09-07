This repository contains analysis scripts for:<br>
Ishigami et al. (2021) A single m<sup>6</sup>A modification in U6 snRNA diversifies exon sequence at the 5’ splice site. (Nature Communications)<br>

Environment:<br>
- python 3.8.5 with packages pysam, Biopython and numpy<br>

Script gtfbamcounts.py counts the reads mapped on splice junctions and classify them depending on their splicing pattern.<br>
The BAM files need to be sorted by position.<br>
- Usage: python3 gtfbamcount.py [BAM file name] [GTF file name] [output bamcount file name] [OPTIONAL:-m if exons on the minus strand are sorted in ascending order]<br>

Script irsseq.py calculates the IRS/PCS scores of each intron based on the count file generated by gtfbamcounts.py and exports them with their sequence based on genomic DNA sequence FASTA file.<br>
- Usage: python3 irsseq.py [bamcount file list: control] [bamcount file list: treated] [genome FASTA file name] [output file name]<br>
- Bamcount filenames are delimited by commas without spaces.<br>

Example code:

```
python3 gtfbamcount.py WT1.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf WT1_count.txt
python3 gtfbamcount.py WT2.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf WT2_count.txt
python3 gtfbamcount.py WT3.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf WT3_count.txt
python3 gtfbamcount.py WT4.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf WT4_count.txt
python3 gtfbamcount.py KO1.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf KO1_count.txt
python3 gtfbamcount.py KO2.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf KO2_count.txt
python3 gtfbamcount.py KO3.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf KO3_count.txt
python3 gtfbamcount.py KO4.bam Schizosaccharomyces_pombe.ASM294v2.46.gtf KO4_count.txt

python3 irsseq.py WT1_count.txt,WT2_count.txt,WT3_count.txt,WT4_count.txt \
                  KO1_count.txt,KO2_count.txt,KO3_count.txt,KO4_count.txt \
                  Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa pseseqcount.txt
```
