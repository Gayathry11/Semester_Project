#!/usr/bin/env python3

import gzip
import random
import numpy as np
import matplotlib.pyplot as plt

# Reproducibility
random.seed(42)

# --------------------------------------------------
# Load reference distribution
# --------------------------------------------------
ref_map = {}

with open("reference.hist") as fh:
    for line in fh:
        length, freq = line.strip().split()
        ref_map[int(length)] = float(freq)

max_len = max(ref_map)

ref_vec = np.zeros(max_len + 1)

for length, freq in ref_map.items():
    ref_vec[length] = freq

# --------------------------------------------------
# First pass: compute query fragment lengths
# --------------------------------------------------
qlens = []

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

qlens = np.array(qlens)

# --------------------------------------------------
# Query distribution
# --------------------------------------------------
uniq, counts = np.unique(qlens, return_counts=True)

qcount = np.zeros(max_len + 1)
qcount[uniq] = counts

qfreq = qcount / qcount.sum()

# --------------------------------------------------
# Sampling probabilities
# --------------------------------------------------
keep_prob = np.zeros(max_len + 1)

mask = qfreq > 0

keep_prob[mask] = np.minimum(
    1.0,
    ref_vec[mask] / qfreq[mask]
)

# --------------------------------------------------
# Second pass: sample BED entries
# --------------------------------------------------
rescaled_lengths = []
kept_lines = []

with gzip.open("query.bed.gz", "rt") as fh:
    for line in fh:

        cols = line.split()

        if len(cols) < 3:
            continue

        start = int(cols[1])
        end = int(cols[2])

        frag_len = end - start

        if frag_len <= max_len:

            if random.random() <= keep_prob[frag_len]:
                kept_lines.append(line)
                rescaled_lengths.append(frag_len)

# --------------------------------------------------
# Write rescaled BED file
# --------------------------------------------------
with gzip.open("rescaled_query.bed.gz", "wt") as out:
    out.writelines(kept_lines)

# --------------------------------------------------
# Calculate final distribution
# --------------------------------------------------
rescaled_lengths = np.array(rescaled_lengths)

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

# --------------------------------------------------
# Plot
# --------------------------------------------------
x = np.arange(max_len + 1)

plt.figure(figsize=(12,6), dpi=150)

plt.plot(
    x,
    h_query,
    color="purple",
    lw=2,
    label=f"Query (n={len(qlens):,})"
)

plt.plot(
    x,
    ref_vec,
    color="green",
    lw=2,
    label="Reference"
)

plt.plot(
    x,
    h_rescaled,
    color="cyan",
    lw=2,
    label=f"Rescaled (n={len(rescaled_lengths):,})"
)

plt.xlabel("Fragment Length")
plt.ylabel("Normalized Frequency")
plt.title("Fragment Length Distribution Rescaling")

plt.legend()
plt.tight_layout()

plt.savefig("rescaled_distribution.png", dpi=300)
plt.show()

print(f"Original fragments : {len(qlens):,}")
print(f"Kept fragments     : {len(rescaled_lengths):,}")
print("BED output         : rescaled_query.bed.gz")
print("Plot output        : rescaled_distribution.png")
