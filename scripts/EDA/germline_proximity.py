
# ====================================================================
# Germline Proximity Analysis 
# ====================================================================

"""
Script: germline_proximity.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants with different oncogenicity distribute in relation to
known pathogenic germline variants 

Major outputs:
--------------
1. Number of variants with a germline distance
2. Descriptive statistics 
3. Distribution of germline distance between classes 
5. Germline distance distributions per gene (genes with most oncogenic variants)

All plots are saved in:
   plots/germline_proximity/

"""

print("-"*50)
print("Variant Germline Proximity Analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import ks_2samp
import os
import argparse
from pathlib import Path

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Explore distance to pathogenic germline variants in somatic variant data."
    ) 

    parser.add_argument(
        "--variants", 
        type=Path, 
        required=False, 
        default="/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_tsg_og.tsv",
        help="Path to the input file with variant data."
    )

    parser.add_argument(
        "--output", 
        type=Path, 
        required=False, 
        default="/home/anekl/git/master/explore_cancer_variants/output/germline_dist_filtered.tsv",
        help="Path to the output file path."
    )

    return parser.parse_args() 


args = getargs() 

# ------------------------------------------------------------
# Create output directory 
# ------------------------------------------------------------

save_dir = "plots/germline_proximity"
os.makedirs(save_dir, exist_ok=True) 

# ------------------------------------------------------------
# Load variant data 
# ------------------------------------------------------------

print(f"Loading variant file:\n{args.variants}\n")

variants = pd.read_csv(
    args.variants,
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} somatic variants.")
print("-"*50) 

# ------------------------------------------------------------
# Find the number of variants with a germline distance 
# ------------------------------------------------------------

variants_with_dist = variants[variants["Germline_Proximity"].notna()]
variants_with_dist["Germline_Proximity"] = pd.to_numeric(variants_with_dist["Germline_Proximity"], errors='coerce')

count_with_distance = len(variants_with_dist)
print(f"Number of variants with germline distances: {count_with_distance:,}")

percent_with_distance = (count_with_distance / len(variants)) * 100
print(f"Percent of variants with germline distances: {percent_with_distance:.2f}%\n")

# ------------------------------------------------------------
# Number of variants with germline distance within each class 
# ------------------------------------------------------------

has_dist_summary = (
  variants_with_dist
  .groupby("ONCOGENIC")["Germline_Proximity"]
  .count()
)

print("Number of variants with a germline distance per class:")
print(has_dist_summary, "\n")

# ------------------------------------------------------------
# Stats by ONCOGENICITY classes 
# ------------------------------------------------------------

stats_summary = variants.groupby("ONCOGENIC")["Germline_Proximity"].agg(
    count="count",
    median="median",
    mean="mean",
    std="std"
).reset_index()

print("-"*50)
print("Descriptive statistics of germline distances:")
print(stats_summary)
print("-"*50)

# ------------------------------------------------------------
# Filtered plot data
# ------------------------------------------------------------

wanted_classes = ["Oncogenic", "Likely Neutral"]

variants_plot = (
    variants_with_dist[variants_with_dist["ONCOGENIC"].isin(wanted_classes)]
    .copy()
)

variants_plot["log_dist"] = np.log10(variants_plot["Germline_Proximity"] + 1)

sns.set_theme(style="whitegrid")
palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

# ------------------------------------------------------------
# Density plot with both classes
# ------------------------------------------------------------

print("\nPlotting distribution of germline proximity (comparison plot)..\n")

plt.figure(figsize=(8, 5))
sns.kdeplot(
    data=variants_plot,
    x="log_dist",
    hue="ONCOGENIC",
    fill=True,
    common_norm=False,
    palette=palette,
    alpha=0.5
)

plt.title("Distribution of Germline Distances (Comparison)", fontsize=14)
plt.xlabel("Distance to nearest pathogenic germline variant (Log10 bp + 1)", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.savefig(f"{save_dir}/combined_dist.png", bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/combined_dist.png'.\n")

# ------------------------------------------------------------
# Histogram with distance proportions per class 
# ------------------------------------------------------------

print("-"*50)
print("Plotting germline distances per class..\n")

sns.set_style("white")

g = sns.FacetGrid(
    variants_plot, col="ONCOGENIC", hue="ONCOGENIC",
    palette=palette, height=4, aspect=1.2, sharey=True
)

g.map(sns.histplot, "log_dist", kde=True, bins=30, alpha=0.4, stat="proportion")

for ax, label in zip(g.axes.flat, wanted_classes):
    subset = variants_plot[variants_plot["ONCOGENIC"] == label]
    median_val = subset["log_dist"].median()
    n_count = len(subset)

    ax.axvline(median_val, color="black", linestyle="--", alpha=0.7)
    actual_median = (10**median_val) - 1
    ax.set_title(f"{label} (n={n_count})\n(Median: {actual_median:.0f} bp)", fontsize=12)

g.set_axis_labels("Log10(Distance + 1)", "Proportion")

sns.despine() 
plt.tight_layout()
plt.savefig(f"{save_dir}/dists_per_class.png", bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/dists_per_class.png'.\n")

# ------------------------------------------------------------
# Boxplot of germline distances 
# ------------------------------------------------------------

print("-"*50)
print("Plotting boxplot of germline distance data..\n")

plt.figure(figsize=(8, 5))
sns.boxplot(
    data=variants_plot,
    x="ONCOGENIC",
    y="log_dist",
    palette=palette
)

plt.title("Boxplot of Distances (Log-scaled)")
plt.ylabel("Log10(Distance + 1)")
plt.tight_layout()
plt.savefig(f"{save_dir}/boxplot_germline_dist.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Boxplot saved as '{save_dir}/boxplot_germline_dist.png'.\n")


# ------------------------------------------------------------
# Find the top genes (genes with most oncogenic variants) with germline distances
# ------------------------------------------------------------

oncogenic_var = variants_plot[variants_plot["ONCOGENIC"] == "Oncogenic"].copy()
top_genes_full = oncogenic_var["Hugo_Symbol"].value_counts()

print("\nThe top 10 genes (with most oncogenic variants) are:")
print(top_genes_full.head(10), "\n")

# ------------------------------------------------------------
# Germline distances in top genes (original data)
# ------------------------------------------------------------

for gene in top_genes_full.head(20).index:
    gene_data = variants_plot[variants_plot["Hugo_Symbol"] == gene]

    if gene_data["ONCOGENIC"].nunique() < 2:
        print(f"Skipping {gene}: Not enough groups to compare.\n")
        continue
    else: 
        print(f"Plotting {gene} germline distances..")

    n_onco = len(gene_data[gene_data["ONCOGENIC"] == "Oncogenic"])
    n_neut = len(gene_data[gene_data["ONCOGENIC"] == "Likely Neutral"])

    plt.figure(figsize=(8, 5))
    sns.set_style("white") 

    ax = sns.kdeplot(
        data=gene_data,
        x="log_dist",
        hue="ONCOGENIC",
        fill=True,
        common_norm=False,
        palette=palette,
        alpha=0.3,
        linewidth=2
    )

    sns.move_legend(ax, "upper right", frameon=False)
    sns.despine() 

    plt.xlim(-1, 5)
    plt.title(f"Germline Proximity: {gene}\n(Oncogenic n={n_onco}, Neutral n={n_neut})", fontsize=14, pad=15)
    plt.xlabel("Distance to nearest pathogenic germline variant (Log10 bp + 1)", fontsize=12, labelpad=10)
    plt.ylabel("Density", fontsize=12, labelpad=10)
    plt.savefig(f"{save_dir}/dist_{gene}.png", bbox_inches="tight")
    plt.show()

    print(f"\nPlotting complete! Plot saved as '{save_dir}/dist_{gene}.png'.\n")

# ------------------------------------------------------------
# Save variants with germline distances in top genes to .tsv 
# ------------------------------------------------------------

print("Saving variants with germline distances in gene with most oncogenic variants to .tsv file..")
top_20_genes = top_genes_full.head(20).index.tolist() 
print(f"The top 20 genes: {top_20_genes}")
top_20_variants = variants_plot[variants_plot["Hugo_Symbol"].isin(top_20_genes)]

# save as .tsv 
top_20_variants.to_csv(args.output, sep="\t", index=False) 
print(f"Filtered variant file saved as: \n {args.output}")
print("-"*30)


print("\nGermline distance exploratory analyses complete!🥳🧬")
