"""
====================================================================
Variants MAVEs Analysis Script 
====================================================================

Script: variants_with_maves.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants with different oncogenicity distribute across 
Multiplexed Assays of Variant Effect (MAVEs)

Major outputs:
--------------
1. Number of variants with MAVEdb scores
2. Number of MaveDB_scores for each oncogenicity class 
3. Top genes with MaveDB_score
4. Filter out wanted MAVE experiments (by number of variants with MAVE score)
5. Box plot of MaveDB score distributions within top MAVE genes 
6. Density plot og MaveDB score distributions within top MAVE genes

All plots are saved in:
   plots/maves

"""

print("-"*50)
print("MAVE score analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------------------------------------
# Load variant data 
# ------------------------------------------------------------

print("Loading variant data..")

variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves_expanded.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} somatic variants.")
print("-"*30)

# ------------------------------------------------------------
# Explore number of variants with MAVE scores
# ------------------------------------------------------------

variants_with_mave = variants[variants["MaveDB_score"].notna()]

print(f"Number of variants with MAVE scores: {len(variants_with_mave):,} variants.") 
print(f"Percentage of variants with MAVE scores: {(len(variants_with_mave)/len(variants)*100):.2f}%.")
print("-"*30)

# ------------------------------------------------------------
# Explore number of variants with MAVE scores for
# different oncogenicity classes 
# ------------------------------------------------------------

# create summary table 
oncogenicity_summary = variants_with_mave.groupby("ONCOGENIC").size().reset_index(name="Variant_Count")

# rename column
oncogenicity_summary.columns = ["Oncogenicity_Class", "Variant_Count"]

# sort by count 
oncogenicity_summary = oncogenicity_summary.sort_values("Variant_Count", ascending = False) 

print("Number of MAVE scores per oncogenicity class:")
print(oncogenicity_summary,"\n") 

# plot summary 
print("Plotting MAVE score per oncogenicity class...")

oncogenicity_palette = {
    "Oncogenic": "#c4314a",
    "Likely Oncogenic": "#D98C6A",
    "Likely Neutral": "#88aed1",
    "Inconclusive": "#f9c74f",
    "Unknown": "#848a8e",
    "Resistance": "#ba7ad4"
}

plt.figure(figsize=(8,5))
sns.barplot(data=oncogenicity_summary, 
            x="Oncogenicity_Class",
            y="Variant_Count",
            hue="Oncogenicity_Class", 
            palette=oncogenicity_palette, 
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Variants with MAVE Scores per Oncogenicity Class", fontsize=14, pad=10)
plt.xlabel("Oncogenicity Class", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12) 
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.tight_layout() 
plt.savefig("plots/maves/oncogenicity_classes_maves.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'plots/maves/oncogenicity_classes_maves.png'\n")
print("-"*30)

# ------------------------------------------------------------
# Explore top genes with MAVEdb_score 
# ------------------------------------------------------------

# summary per gene 
gene_summary = variants_with_mave.groupby("Hugo_Symbol").size().reset_index(name="Variant_Count")

# rename columne
gene_summary.columns = ["Gene", "Variant_Count"]

# sort by count 
gene_summary = gene_summary.sort_values("Variant_Count", ascending = False) 

print("Number of MaveDB_scores per gene:")
print(gene_summary, "\n") 

# plot summary 
print("Plotting top genes with MaveDB score (by variant count)..")

plt.figure(figsize=(8,5))
sns.barplot(data=gene_summary.head(10), 
            x="Gene",
            y="Variant_Count",
            color="#c4314a",
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Top Genes with MAVE Scores", fontsize=14, pad=10)
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12) 
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.tight_layout() 
plt.savefig("plots/maves/top_genes_with_maves.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'plots/maves/top_genes_with_maves.png'")
print("-"*30)

# ------------------------------------------------------------
# Find oncogenicity distribution within genes with MaveDB_score 
# ------------------------------------------------------------

oncogenicity_classes = ["Likely Neutral", "Oncogenic"] 
variants_with_mave_filtered = variants_with_mave[variants_with_mave["ONCOGENIC"].isin(oncogenicity_classes)]

mave_summary = variants_with_mave_filtered.groupby(["Hugo_Symbol", "MaveDB_urn", "ONCOGENIC"]).size()
print("Variant count for different MAVE experiments in different oncogenicity classes:")
print(mave_summary)

# -------------------------------------------------------------
# Explore MAVE genes 
# -------------------------------------------------------------

print("Filtering out experiments from the top MAVE genes..")
# filter to wanted MAVE experiments 
wanted_urns = ["urn:mavedb:00000068-a-1", "urn:mavedb:00000013-a-1", "urn:mavedb:00000115-a-7", "urn:mavedb:00000081-a-2"] 

mave_filtered = variants_with_mave_filtered[variants_with_mave_filtered["MaveDB_urn"].isin(wanted_urns)]

mave_filtered["MaveDB_score"] = pd.to_numeric(
    mave_filtered["MaveDB_score"], errors="coerce"
)
# Drop rows without MAVE score
mave_filtered = mave_filtered.dropna(subset=["MaveDB_score"])

# -------------------------------------------------------------
# Plot score distributions for wanted MAVE experiments
# -------------------------------------------------------------

print("Plotting score distributions for wanted MAVE experiments..")

palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

# --- BOXPLOT ---

g = sns.FacetGrid(mave_filtered, col="Hugo_Symbol", sharey=False, sharex=False)
g.map_dataframe(sns.boxplot, x="ONCOGENIC", y="MaveDB_score", 
                order=oncogenicity_classes,
                palette=palette,
                hue="ONCOGENIC",
                hue_order=oncogenicity_classes)

# Titles and axis labels
g.set_titles(col_template="{col_name}") 
g.set_axis_labels(x_var="Oncogenicity", y_var="MAVE Functional Score")
g.figure.suptitle("MAVE Score Distribution by Oncogenicity Class", 
                   fontsize=14, y=1.02)

plt.tight_layout()
plt.savefig("plots/maves/mave_score_by_oncogenicity.png", dpi=300, bbox_inches="tight")


# --- DENSITY PLOT ---

g = sns.FacetGrid(mave_filtered, col="Hugo_Symbol", sharey=False, sharex=False, 
                  hue="ONCOGENIC", hue_order=oncogenicity_classes, palette=palette)
g.map(sns.kdeplot, "MaveDB_score", fill=True, alpha=0.5)

# Add sample size annotations per panel
for ax, gene in zip(g.axes.flat, g.col_names):
    gene_data = mave_filtered[mave_filtered["Hugo_Symbol"] == gene]
    n_neutral = (gene_data["ONCOGENIC"] == "Likely Neutral").sum()
    n_oncogenic = (gene_data["ONCOGENIC"] == "Oncogenic").sum()
    
    ax.annotate(
        f"n(Likely Neutral) = {n_neutral}\nn(Oncogenic) = {n_oncogenic}",
        xy=(0.05, 0.95),           # top-left corner of panel
        xycoords="axes fraction",
        fontsize=6,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.5)
    )

# Titles and axis labels
g.set_titles(col_template="{col_name}")
g.set_axis_labels(x_var="MAVE Functional Score", y_var="Density")
g.figure.suptitle("MAVE Score Distribution by Oncogenicity Class",
                   fontsize=14, y=1.02)
g.add_legend(title="Oncogenicity", fontsize=9)

plt.tight_layout()
plt.savefig("plots/maves/mave_score_density_by_oncogenicity.png", dpi=300, bbox_inches="tight")

print("Plotting complete. Plots saved in folder plots/maves/.")


print("MAVE exploratory analyses complete!🥳🥳")