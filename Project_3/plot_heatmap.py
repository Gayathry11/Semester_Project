import pandas as pd
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

plt.figure(figsize=(12,8))

sns.heatmap(
    pivot,
    cmap="Purples",
    cbar_kws={"label":"Fragment center count"}
)

plt.xlabel("Distance from motif center (bp)")
plt.ylabel("cfDNA fragment length (bp)")
plt.title("V-Plot Heatmap")

plt.tight_layout()
plt.savefig("vplot_heatmap.png", dpi=300)
