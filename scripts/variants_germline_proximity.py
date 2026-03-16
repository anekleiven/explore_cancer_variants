"""
Variants Germline Proximity Analysis 
------------------------------------------------------

Script: variants_with_germline_proximity.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants with different oncogenicity distribute in relation to
known pathogenic germline variants 

Major outputs:
--------------
1. Number of variants with a germline distance
2. Percent of variants with a germline distance
3. Percent of variants with a germline distance (per class) 
4. Simple descriptive statistics 
5. Distribution of germline distance between classes 
6. Boxplot to spot outliers and look at distribution 
7. Bi-modal distribution study oncogenic variants
8. Germline distances per gene (top oncogenic genes)

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

# ------------------------------------------------------------
# Load variant data 
# ------------------------------------------------------------

print("Loading variant data..")

variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} somatic variants.")

# ------------------------------------------------------------
# Find the number of variants with a germline distance 
# ------------------------------------------------------------

count_with_distance = variants["Germline_Proximity"].notna().sum()
print(f"Number of variants with germline distances: {count_with_distance:,}")

percent_with_distance = (count_with_distance / len(variants)) * 100
print(f"Percent of variants with germline distances: {percent_with_distance:.2f}%\n")

# ------------------------------------------------------------
# Number of variants with germline distance within each class 
# ------------------------------------------------------------

has_dist_summary = variants.groupby("ONCOGENIC")["Germline_Proximity"].apply(
    lambda x: x.notna().mean() * 100
)
print("Percent of variants with a germline distance per class:")
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

print("Descriptive statistics germline distances:")
print(stats_summary)

# ------------------------------------------------------------
# Prepare filtered plot data (used across all figures below)
# ------------------------------------------------------------

wanted_classes = ["Oncogenic", "Likely Neutral"]

variants_plot = (
    variants[variants["ONCOGENIC"].isin(wanted_classes)]
    .dropna(subset=["Germline_Proximity"])
    .copy()
)

variants_plot["log_dist"] = np.log10(variants_plot["Germline_Proximity"] + 1)

sns.set_theme(style="whitegrid")
palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

# ------------------------------------------------------------
# FIGURE 1: Density plot with both classes
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
plt.savefig("plots/germline_proximity/combined_dist.png", bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/germline_proximity/combined_dist.png'.\n")

# ------------------------------------------------------------
# FIGURE 2: Histogram with distance proportions per class 
# ------------------------------------------------------------

print("-"*50)
print("Plotting germline distances per class..\n")

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
    ax.set_title(f"{label} (n={n_count})\n(Median: 10^{median_val:.1f} bp)", fontsize=12)

g.set_axis_labels("Log10(Distance + 1)", "Proportion")

plt.tight_layout()
plt.savefig("plots/germline_proximity/dists_per_class.png", bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/germline_proximity/dists_per_class.png'.\n")

# ------------------------------------------------------------
# FIGURE 3: Boxplot of germline distances 
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
plt.savefig("plots/germline_proximity/boxplot_germline_dist.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Boxplot saved as 'plots/germline_proximity/boxplot_germline_dist.png'.\n")

# ------------------------------------------------------------
# Statistical test: Kolmogorov-Smirnov 
# ------------------------------------------------------------

# Model assumptions:
#       The samples should be independent.       
#       The dependent variable must be measured on an ordinal or continuous scale. 
#       The distributions should be continuous to avoid ties. 

# Hypotheses: 
#       H0: The samples come from the same distribution.
#       H1: The samples come from different distributions. 

print("-"*50)
print("Running Kolmogorov-Smirnov test on germline distance data (log-transformed)..\n")

oncogenic = variants_plot[variants_plot["ONCOGENIC"] == "Oncogenic"]["log_dist"]
neutral = variants_plot[variants_plot["ONCOGENIC"] == "Likely Neutral"]["log_dist"]

ks_statistic, p_value = ks_2samp(oncogenic, neutral, alternative="two-sided", mode="auto")

alpha = 0.05

print("Results KS-test:")
print(f"KS-statistic: {ks_statistic}")
print(f"p-value: {p_value:.4f}")

if p_value < alpha:
    print("Reject the null hypothesis. The oncogenic and neutral variants come from two different distributions.")
else:
    print("Failed to reject the null hypothesis.")

# ------------------------------------------------------------
# Investigate the bi-modal distribution for
# germline distances among oncogenic variants 
# ------------------------------------------------------------

print("\nInvestigating the bi-modal distribution for germline distances among oncogenic variants..")

oncogenic_var = variants_plot[variants_plot["ONCOGENIC"] == "Oncogenic"].copy()
peak_a_variants = oncogenic_var[oncogenic_var["log_dist"] <= 1].copy()
peak_b_variants = oncogenic_var[oncogenic_var["log_dist"] >= 1.5].copy()

top_genes_full = oncogenic_var["Hugo_Symbol"].value_counts()
genes_peak_a   = peak_a_variants["Hugo_Symbol"].value_counts()
genes_peak_b   = peak_b_variants["Hugo_Symbol"].value_counts()

print("\nThe top 10 oncogenic genes total are:")
print(top_genes_full.head(10), "\n")
print("The top 10 oncogenic genes in peak A are:")
print(genes_peak_a.head(10), "\n")
print("The top 10 oncogenic genes in peak B are:")
print(genes_peak_b.head(10), "\n")

# ------------------------------------------------------------
# Proportion plot germline distances
# ------------------------------------------------------------

print("Plotting distribution of variants in peak A and B for top oncogenic genes..")

peak_a_variants["Peak"] = "Peak A (<10bp)"
peak_b_variants["Peak"] = "Peak B (>30bp)"

combined_df = pd.concat([peak_a_variants, peak_b_variants])
top_15_genes = combined_df["Hugo_Symbol"].value_counts().head(15).index
filtered_df = combined_df[combined_df["Hugo_Symbol"].isin(top_15_genes)]

peak_counts = (
    filtered_df
    .groupby(["Hugo_Symbol", "Peak"])
    .size()
    .unstack(fill_value=0)
)

peak_proportions = peak_counts.div(peak_counts.sum(axis=1), axis=0)
peak_proportions = peak_proportions.sort_values("Peak A (<10bp)", ascending=False)

x_labels = [f"{gene} (n={int(peak_counts.loc[gene].sum())})" for gene in peak_proportions.index]

ax = peak_proportions.plot(
    kind="bar",
    stacked=True,
    figsize=(8, 5),
    color=("#d16a58", "#d6b4be"),
    edgecolor="white"
)

plt.title("Proportion of Oncogenic Variants in Distance Peaks (Top 15 Genes)", fontsize=14, pad=20)
plt.xlabel("Gene (Total variants in peaks)", fontsize=12)
plt.ylabel("Proportion of variants", fontsize=12)
plt.xticks(range(len(x_labels)), x_labels, rotation=45, ha="right", fontsize=10)
plt.legend(title="Distance to Germline", loc="upper left", bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig("plots/germline_proximity/germline_proportions.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/germline_proximity/germline_proportions.png'\n")

# ------------------------------------------------------------
# Germline distances in top genes (original data)
# ------------------------------------------------------------

for gene in top_genes_full.head(15).index:
    gene_data = variants_plot[variants_plot["Hugo_Symbol"] == gene]

    if gene_data["ONCOGENIC"].nunique() < 2:
        print(f"Skipping {gene}: Not enough groups to compare.\n")
        continue
    else: 
        print(f"Plotting {gene} germline distances..")

    n_onco = len(gene_data[gene_data["ONCOGENIC"] == "Oncogenic"])
    n_neut = len(gene_data[gene_data["ONCOGENIC"] == "Likely Neutral"])

    plt.figure(figsize=(8, 5))

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

    sns.move_legend(ax, "upper right")
    plt.xlim(-1, 5)
    plt.title(f"{gene}: Germline Proximity\n(Oncogenic n={n_onco}, Neutral n={n_neut})", fontsize=14)
    plt.xlabel("Distance to nearest pathogenic germline variant (Log10 bp + 1)", fontsize=12)
    plt.ylabel("Density", fontsize=12)
    plt.savefig(f"plots/germline_proximity/dist_{gene}.png", bbox_inches="tight")
    plt.show()

    print(f"Plotting complete! Plot saved as 'plots/germline_proximity/dist_{gene}.png'.\n")

# ------------------------------------------------------------
# Germline distances in top genes
# Original data combined with synonymous variant data
# ------------------------------------------------------------

# Load neutral data 
neutral_ref = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/synonymous_dataset.tsv",
    sep="\t",
    low_memory=False
)

palette_3 = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
    "Synonymous": "#a5a8aa"
}

print("Plotting distribution with extended synonymous variant group..")

for gene in top_genes_full.head(15).index:

    onco_data = variants_plot[
        (variants_plot["Hugo_Symbol"] == gene) &
        (variants_plot["ONCOGENIC"] == "Oncogenic")
    ].copy()
    onco_data["Group"] = "Oncogenic"

    likely_neut_data = variants_plot[
        (variants_plot["Hugo_Symbol"] == gene) &
        (variants_plot["ONCOGENIC"] == "Likely Neutral")
    ].copy()
    likely_neut_data["Group"] = "Likely Neutral"

    synonymous_data = neutral_ref[neutral_ref["Hugo_Symbol"] == gene].copy()
    synonymous_data["Group"] = "Synonymous"

    gene_plot_df = (
        pd.concat([onco_data, likely_neut_data, synonymous_data])
        .dropna(subset=["Germline_Proximity"])
        .copy()
    )

    if gene_plot_df.empty or gene_plot_df["Group"].nunique() < 2:
        print(f"Skipping {gene}: Too few groups to compare.\n")
        continue
    else: 
        print(f"Plotting germline distances for {gene}..")

    counts = gene_plot_df["Group"].value_counts()
    n_onco   = counts.get("Oncogenic", 0)
    n_l_neut = counts.get("Likely Neutral", 0)
    n_s_neut = counts.get("Synonymous", 0)

    plt.figure(figsize=(10, 6))

    ax = sns.kdeplot(
        data=gene_plot_df,
        x="log_dist",
        hue="Group",
        fill=True,
        common_norm=False,
        palette=palette_3,
        alpha=0.3,
        linewidth=2
    )

    sns.move_legend(ax, "upper right", title="Variant Category")
    plt.xlim(-1, 5)
    plt.title(
        f"{gene}: Germline Proximity Comparison\n"
        f"Onco={n_onco}, L.Neut={n_l_neut}, Syn={n_s_neut}",
        fontsize=13
    )
    plt.xlabel("Distance to nearest pathogenic germline variant (Log10 bp + 1)", fontsize=12)
    plt.ylabel("Density", fontsize=12)

    plt.tight_layout()
    plt.savefig(f"plots/germline_proximity/dist_extended_{gene}.png", dpi=300, bbox_inches="tight")
    plt.show()

    print(f"Plotting complete! Plot saved as 'plots/germline_proximity/dist_extended_{gene}.png'.\n")
