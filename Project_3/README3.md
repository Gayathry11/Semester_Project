Project 3: V-Plot Generation

Overview
This project generates a V-plot heatmap from mapped cfDNA fragment data.

The pipeline reads fragment and motif data from standard input, calculates fragment-center offsets relative to motif centers, generates a long-form matrix, and visualizes the matrix as a heatmap.

The long-form matrix stores counts for each combination of offset and fragment length.

Input Data
Input file:
mapped.bed.gz
The input contains motif coordinates together with mapped cfDNA fragments.

Files
Matrix Generation Script
project3_final.py
Reads data from standard input and creates a long-form matrix.

Heatmap Generation Script
plot_heatmap_v2.py
Reads the matrix file and generates the V-plot heatmap.

Matrix Output
matrix_long_form.tsv
Format:
offset    fragment_length    count

Example:

-500    50    24
-500    51    39
-500    52    27

Heatmap Output
vplot_heatmap_v2.png

Methodology

Step 1
Read intersected fragment data from standard input.

Step 2
Calculate motif midpoint:
motif_mid = (motif_start + motif_end) // 2

Step 3
Calculate fragment midpoint:
frag_mid = (frag_start + frag_end) // 2

Step 4
Calculate offset:
offset = frag_mid - motif_mid

Step 5
Count occurrences for each:
(offset, fragment_length)
pair.

Step 6
Write counts to a long-form matrix file.

Step 7
Generate a heatmap using Pandas, Matplotlib, and Seaborn.

The heatmap uses logarithmic count scaling to improve visualization and highlight enrichment patterns across fragment lengths and motif offsets.

Generate Matrix File
zcat mapped.bed.gz | python3 project3_final.py > matrix_long_form.tsv

Generate Heatmap
python3 plot_heatmap_v2.py

Output:
vplot_heatmap_v2.png

View Heatmap
explorer.exe vplot_heatmap_v2.png

Output Files
project3_final.py
plot_heatmap_v2.py
matrix_long_form.tsv
vplot_heatmap_v2.png
readme3.md

Dependencies
Install required libraries:
python3 -m pip install pandas matplotlib seaborn

Results
A long-form matrix containing 126,124 offset-length observations was generated and visualized as a heatmap.

The plot displays the distribution of cfDNA fragment lengths relative to motif-centered offsets and provides a graphical representation of fragment occupancy patterns around motif locations.

Author
Gayathry Vijayan

Semester Project

Project 3 – V-Plot Generation
