Project 2: Fragment Length Distribution Rescaling
Overview
This project performs fragment length distribution rescaling of a query BED dataset to match a given reference fragment length distribution.

The objective is to normalize fragment-size biases between datasets by selecting fragments from the query dataset such that the resulting fragment length distribution closely follows the reference distribution.

The implementation uses exact length-wise sampling rather than probabilistic downsampling, resulting in an almost perfect agreement with the reference profile.

Input Files
Query Dataset
query.bed.gz
Compressed BED file containing genomic fragments.

Example:

chr1    1000    1150
chr1    2500    2660
chr2    8000    8150
Fragment length is calculated as:

Length = End - Start
Reference Distribution
reference.hist
Histogram containing fragment lengths and their frequencies.

Example:

50   0.0012
51   0.0015
52   0.0018
Format:

Fragment_Length Frequency
Methodology
Step 1: Load Reference Distribution
The reference histogram is read and converted into a NumPy vector.

The distribution is normalized using:

ref_vec = ref_vec / ref_vec.sum()
This ensures all frequencies sum to 1.

Step 2: Read Query Fragments
The query BED file is scanned and fragment lengths are calculated:

frag_len = end - start
Fragments are grouped according to their lengths.

Step 3: Compute Query Distribution
The frequency of each fragment length is calculated using:

np.bincount()
The counts are then normalized to obtain the query distribution.

Step 4: Determine Maximum Scaling Factor
Since new fragments cannot be created, only existing fragments can be retained.

The script computes the largest possible scaling factor:

scale_factor = np.min(
    available[valid] / ref_vec[valid]
)
This guarantees that:

Target Count ≤ Available Count
for every fragment length.

Step 5: Exact Length-Wise Sampling
For each fragment length, the exact required number of fragments is selected using:

random.sample()
This eliminates the variability caused by probabilistic sampling and produces a much closer match to the reference distribution.

Step 6: Generate Rescaled Dataset
Selected fragments are written to:

rescaled_query.bed.gz
Step 7: Visualization
Three distributions are displayed in a comparison plot:

Original Query Distribution
Reference Distribution
Rescaled Distribution
Output:

rescaled_distribution.png
Running the Script
Execute:

python3 rescale_script6.py
Output Files
Rescaled BED File
rescaled_query.bed.gz
Contains the final rescaled fragment dataset.

Distribution Plot
rescaled_distribution.png
Provides a visual comparison of:

Query distribution
Reference distribution
Rescaled distribution
Example Output
Original fragments : 13,835,737
Kept fragments     : 1,220,621
Correlation        : 1.000000
RMSE               : 0.00000097

BED output         : rescaled_query.bed.gz
Plot output        : rescaled_distribution.png
Performance Metrics
Correlation
Measures similarity between the rescaled and reference distributions.

1.0 = Perfect Agreement
Observed:

Correlation = 1.000000
RMSE (Root Mean Square Error)
Measures the difference between the reference and rescaled distributions.

0 = Perfect Match
Observed:

RMSE = 0.00000097
Results
The original query dataset exhibited a fragment length distribution that differed substantially from the reference distribution.

After exact length-wise rescaling:

The rescaled distribution closely matched the reference distribution.
Correlation reached 1.000000.
RMSE approached zero.
Visual inspection confirmed near-perfect overlap between the reference and rescaled distributions.
These results demonstrate successful fragment length normalization.

Requirements
Install dependencies:

pip install numpy matplotlib
Files Generated
query.bed.gz
reference.hist
rescale_script6.py
rescaled_query.bed.gz
rescaled_distribution.png
README.md
Author
Gayathry Vijayan

Semester Project

Fragment Length Distribution Rescaling Using Exact Length-Wise Sampling
