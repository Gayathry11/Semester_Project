import sys
from collections import defaultdict

counts = defaultdict(int)

for line in sys.stdin:
    f = line.strip().split()

    if len(f) < 12:
        continue

    motif_start = int(f[2])
    motif_end = int(f[3])

    frag_start = int(f[8])
    frag_end = int(f[9])

    frag_len = int(f[11])

    motif_mid = (motif_start + motif_end) // 2
    frag_mid = (frag_start + frag_end) // 2

    offset = frag_mid - motif_mid

    if -500 <= offset <= 500 and 50 <= frag_len <= 175:
        counts[(offset, frag_len)] += 1

for (offset, frag_len), count in sorted(counts.items()):
    print(offset, frag_len, count, sep="\t")
