Project 2: Rescale Query Fragment Length Distribution to a Reference Distribution

Overview
This project rescales the fragment length distribution of a query BED file so that it more closely matches a given reference distribution. The method uses probabilistic downsampling, where fragments are retained according to probabilities derived from the ratio of reference and query frequencies.

The output consists of:

A rescaled BED file
A comparison plot showing the original query, reference, and rescaled distributions

Input Files
1. reference.hist
A tab-separated file containing the reference fragment length distribution.

Example:

35 0.00012
36 0.00015
37 0.00017
...
Column 1: Fragment length
Column 2: Normalized frequency

2. query.bed.gz
Compressed BED file containing genomic fragments.

Example:

chr1    1000    1170
chr1    2500    2665

Fragment length is calculated as:

Length = End - Start

Methodology

Step 1: Read Reference Distribution
The reference fragment length frequencies are loaded from reference.hist.

Step 2: Compute Query Fragment Length Distribution
Fragment lengths are extracted from query.bed.gz and converted into a frequency distribution.

Step 3: Calculate Retention Probability
For each fragment length:

keep_probability = min(1.0, reference_frequency / query_frequency)

This ensures that overrepresented lengths are downsampled while underrepresented lengths are retained.

Step 4: Resample Query Fragments
Each fragment is independently retained according to its calculated probability.

Step 5: Generate Outputs
The retained fragments are written to a new compressed BED file and the resulting distribution is plotted against the original query and reference distributions.

Requirements

Install dependencies:
pip install numpy matplotlib

How to Run
Place the following files in the same directory:

reference.hist
query.bed.gz
rescale_script2.py

Run:
python3 rescale_script2.py

Output Files

Rescaled BED File
rescaled_query.bed.gz
Contains all fragments retained after resampling.

Distribution Plot
rescaled_distribution.png

Visual comparison of:
Original query distribution
Reference distribution
Rescaled distribution

Example Run
Original fragments : 13,835,737
Kept fragments     : 7,198,965
BED output         : rescaled_query.bed.gz
Plot output        : rescaled_distribution.png

Repository Structure
Project_2/
│
├── rescale_script2.py
├── README.md
├── reference.hist
├── query.bed.gz
├── rescaled_query.bed.gz
└── rescaled_distribution.png

Notes
The rescaling method uses downsampling only.
No artificial fragments are created.
Since only existing fragments can be retained, the rescaled distribution may not perfectly overlap the reference distribution at every fragment length.
The generated plot is used to visually evaluate the quality of the rescaling.

Author
Gayathry Vijayan

M.Tech Structural and Computational Biology
Indian Institute of Technology Roorkee
