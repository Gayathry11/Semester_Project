Project 3 - V Plot Generation

Description
This project generates a V-plot heatmap from mapped cfDNA fragment data.

Input
mapped.bed.gz

Files
project3_final.py
plot_heatmap.py
vplot_heatmap.png

Generate long-form matrix
zcat mapped.bed.gz | python3 project3_final.py > matrix_long_form.tsv

Generate heatmap
python plot_heatmap.py

Output
vplot_heatmap.png
