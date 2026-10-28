Project 1: Quick BED Sort

Objective
The goal of this project is to filter and sort BED records according to a chromosome selection file.

Input Files
selection.tsv
shuf.a.bed
shuf.b.bed

Requirements
Python 3

Create Combined BED File
cat shuf.a.bed/shuf.a.bed shuf.b.bed/shuf.b.bed > combined.bed

Run
python3 project1_script6.py selection.tsv < combined.bed > output.bed

Algorithm
Read chromosome names from selection.tsv.
Store chromosome order in a list.
Create a separate container for each chromosome.
Read BED records from standard input.
Filter chromosomes not present in the selection file.
Store records chromosome-wise.
Sort records within each chromosome using genomic start coordinates.
Print chromosomes in the same order as specified in the selection file.

Output
output.bed

Performance
Runtime on the provided dataset:

real 1m20.757s

Validation
Chromosome order can be checked using:
cut -f1 output.bed | uniq

and compared with:
cat selection.tsv

Author
Gayathry Vijayan

M.Tech Structural and Computational Biology

IIT Roorkee
