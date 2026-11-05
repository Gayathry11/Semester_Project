import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv(
    "matrix_long_form.tsv",
    sep="\t",
    header=None,
    names=["offset", "length", "count"]
)

pivot = df.pivot(
    index="length",
    columns="offset",
    values="count"
).fillna(0)

pivot = pivot.sort_index(
    ascending=False
)

pivot = np.log1p(pivot)

plt.figure(figsize=(12,8))

sns.heatmap(
    pivot,
    cmap="magma",
    cbar_kws={"label":"log(Fragment count)"}
)

plt.xlabel("Distance from motif center (bp)")
plt.ylabel("Fragment Length (bp)")
plt.title("cfDNA V-Plot Heatmap")

plt.tight_layout()

plt.savefig(
    "vplot_heatmap_v2.png",
    dpi=300
)

plt.show()
