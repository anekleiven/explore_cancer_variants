# ============================================================
# Germline Proximity Analysis
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load variant data 
print("Loading somatic variants..")
somatic = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t", low_memory=False
)
print(f"Loaded {len(somatic):,} variants.\n")

# ── 2. Rebuild germline_per_gene ──────────────────────────────
print("Loading germline variants...")
germline_cols = ["Assembly", "Origin", "Name", "GeneID", "ClinicalSignificance"]
somatic_gene_ids = set(somatic["Entrez_Gene_Id"].dropna().unique())

chunks = []
for chunk in pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/data/variant_summary.txt.gz",
    sep="\t", usecols=germline_cols, chunksize=100000, low_memory=False
):
    chunks.append(chunk[chunk["GeneID"].isin(somatic_gene_ids)])

germline = pd.concat(chunks, ignore_index=True)

# filter to pathogenic germline variants from GRCh37
germline = germline[
    (germline["Assembly"] == "GRCh37") &
    (germline["Origin"] == "germline") &
    (germline["ClinicalSignificance"].isin(["Pathogenic", "LikelyPathogenic"]))
]

# extract amino acid position
germline["AA_Position"] = (
    germline["Name"]
    .str.split(" ", n=1).str[1]
    .str.strip("()")
    .str.extract(r"p\.[A-Za-z]{3}(\d+)[A-Za-z]{3}")
    .astype(float)
)

germline = germline.dropna(subset=["AA_Position"])
germline = germline.rename(columns={"GeneID": "Entrez_Gene_Id"})

# build dictionary: gene_id -> list of germline AA positions
germline_per_gene = (
    germline.groupby("Entrez_Gene_Id")["AA_Position"]
    .apply(list)
    .to_dict()
)
print(f"Built germline dictionary for {len(germline_per_gene):,} genes.\n")

# ── 3. Protein lengths (approximated from max observed position) ──
protein_lengths = (
    somatic.groupby("Entrez_Gene_Id")["Protein_position"]
    .max()
    .to_dict()
)

# ── 4. Compute N(d) and K(d) ─────────────────────────────────
max_dist = 300
N = np.zeros(max_dist + 1)  # background: all residues
K = np.zeros(max_dist + 1)  # observed: oncogenic variants

oncogenic_variants = somatic[somatic["ONCOGENIC"] == "Oncogenic"]

print("Computing enrichment...")
for gene_id, germline_positions in germline_per_gene.items():

    L = protein_lengths.get(gene_id)
    if L is None or pd.isna(L):
        continue

    L = int(L)
    germline_arr = np.array(germline_positions)

    # N(d): background distribution across all residues
    all_positions = np.arange(1, L + 1)
    dists = np.abs(all_positions[:, None] - germline_arr[None, :]).min(axis=1)
    dists = np.clip(dists, 0, max_dist).astype(int)
    np.add.at(N, dists, 1)

    # K(d): oncogenic variants in this gene
    gene_onco = oncogenic_variants[
        (oncogenic_variants["Entrez_Gene_Id"] == gene_id) &
        (oncogenic_variants["Protein_position"].notna())
    ]["Protein_position"].values

    if len(gene_onco) == 0:
        continue

    onco_dists = np.abs(gene_onco[:, None] - germline_arr[None, :]).min(axis=1)
    onco_dists = np.clip(onco_dists.astype(int), 0, max_dist)
    np.add.at(K, onco_dists, 1)



# ── 5. Compute enrichment E(d) = K(d) / N(d) ─────────────────
with np.errstate(divide='ignore', invalid='ignore'):
    E = np.where(N > 0, K / N, np.nan)

# ── 6. Plot ───────────────────────────────────────────────────
distances = np.arange(max_dist + 1)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].plot(distances, K, color="salmon")
axes[0].set_title("Raw distribution\n(K(d): oncogenic variants per distance)")
axes[0].set_xlabel("Distance to nearest pathogenic germline variant (aa)")
axes[0].set_ylabel("Count")

axes[1].plot(distances, E, color="steelblue")
axes[1].axhline(y=np.nanmean(E), color="gray", linestyle="--", label="Mean")
axes[1].set_title("Normalised enrichment\n(E(d) = K(d) / N(d))")
axes[1].set_xlabel("Distance to nearest pathogenic germline variant (aa)")
axes[1].set_ylabel("Enrichment (observed/expected)")
axes[1].legend()

plt.tight_layout()
plt.savefig("plots/germline_proximity_enrichment.png", dpi=150)
plt.show()

print("Done! 🎉")