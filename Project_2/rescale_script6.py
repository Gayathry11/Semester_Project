#!/usr/bin/env python3

import gzip
import random
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# ==================================================
# Reproducibility
# ==================================================

random.seed(42)

# ==================================================
# Load reference distribution
# ==================================================

ref_map = {}

with open("reference.hist") as fh:
    for line in fh:
        length, freq = line.strip().split()
        ref_map[int(length)] = float(freq)

max_len = max(ref_map)

ref_vec = np.zeros(max_len + 1)

for length, freq in ref_map.items():
    ref_vec[length] = freq

# Normalize reference distribution
ref_vec = ref_vec / ref_vec.sum()

# ==================================================
# Read query fragments
# ==================================================

qlens = []
fragments_by_length = defaultdict(list)

with gzip.open("query.bed.gz", "rt") as fh:

    for line in fh:

        cols = line.split()

        if len(cols) < 3:
            continue

        start = int(cols[1])
        end = int(cols[2])

        frag_len = end - start

        if 0 < frag_len <= max_len:
            qlens.append(frag_len)
            fragments_by_length[frag_len].append(line)

qlens = np.array(qlens)

# ==================================================
# Query distribution
# ==================================================

qcount = np.bincount(
    qlens,
    minlength=max_len + 1
)

qfreq = qcount / qcount.sum()

# ==================================================
# Determine maximum scaling factor
# ==================================================

available = qcount.astype(float)

valid = (
    (ref_vec > 0) &
    (available > 0)
)

scale_factor = np.min(
    available[valid] / ref_vec[valid]
)

target_counts = np.floor(
    ref_vec * scale_factor
).astype(int)

# ==================================================
# Exact per-length sampling
# ==================================================

kept_lines = []
rescaled_lengths = []

for frag_len, lines in fragments_by_length.items():

    desired = target_counts[frag_len]

    desired = min(desired, len(lines))

    if desired > 0:

        sampled = random.sample(
            lines,
            desired
        )

        kept_lines.extend(sampled)

        rescaled_lengths.extend(
            [frag_len] * desired
        )

# ==================================================
# Write BED output
# ==================================================

with gzip.open(
    "rescaled_query.bed.gz",
    "wt"
) as out:
    out.writelines(kept_lines)

# ==================================================
# Distributions
# ==================================================

rescaled_lengths = np.array(
    rescaled_lengths
)

h_query = np.bincount(
    qlens,
    minlength=max_len + 1
)

h_rescaled = np.bincount(
    rescaled_lengths,
    minlength=max_len + 1
)

h_query = h_query / h_query.sum()
h_rescaled = h_rescaled / h_rescaled.sum()

# ==================================================
# Statistics
# ==================================================

correlation = np.corrcoef(
    ref_vec,
    h_rescaled
)[0, 1]

rmse = np.sqrt(
    np.mean(
        (ref_vec - h_rescaled) ** 2
    )
)

# ==================================================
# Plot
# ==================================================

x = np.arange(max_len + 1)

plt.figure(
    figsize=(12, 6),
    dpi=150
)

# Original Query
plt.plot(
    x,
    h_query,
    color="purple",
    linewidth=2,
    label=f"Query (n={len(qlens):,})"
)

# Reference as red circles
plt.plot(
    x,
    ref_vec,
    'o',
    color='red',
    markersize=3,
    alpha=0.8,
    label='Reference'
)

# Rescaled as blue dashed line
plt.plot(
    x,
    h_rescaled,
    color='blue',
    linewidth=2,
    linestyle='--',
    label=f"Rescaled (n={len(rescaled_lengths):,})"
)

plt.xlabel(
    "Fragment Length"
)
plt.ylabel(
    "Normalized Frequency"
)

plt.title(
    "Fragment Length Distribution Rescaling"
)

plt.legend()

plt.grid(
    alpha=0.3
)

plt.tight_layout()

plt.savefig(
    "rescaled_distribution.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ==================================================
# Report
# ==================================================

print()
print(f"Original fragments : {len(qlens):,}")
print(f"Kept fragments     : {len(rescaled_lengths):,}")
print(f"Correlation        : {correlation:.6f}")
print(f"RMSE               : {rmse:.8f}")
print()
print("BED output         : rescaled_query.bed.gz")
print("Plot output        : rescaled_distribution.png")
