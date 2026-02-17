"""
====================================================================
Variants Germline Proximity Analysis 
====================================================================

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

All plots are saved in:
   plots/

"""


print("\n========================================================")
print("VARIANT GERMLINE PROXIMITY ANALYSIS")
print("========================================================")

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ------------------------------------------------------------
# Load variant data 
# ------------------------------------------------------------

print("\nLoading variant data...\n")

variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} somatic variants.\n")

# ------------------------------------------------------------
# Find the number of variants with a germline distance 
# ------------------------------------------------------------

count_with_distance = variants["Germline_Proximity"].notna().sum()
print(f"Number of variants with germline distances: {count_with_distance:,}")

percent_with_distance = (count_with_distance/len(variants))*100 
print(f"Percent of variants with germline distances: {percent_with_distance:.2f}%\n")

# ------------------------------------------------------------
# Number of variants with germline distance within each class 
# ------------------------------------------------------------

variants["has_distance"] = variants["Germline_Proximity"].notna() 

has_dist_summary = variants.groupby('ONCOGENIC')['has_distance'].mean() * 100
print("Percent of variants with a germline distance per class:")
print(has_dist_summary, "\n") 

# ------------------------------------------------------------
# Stats by ONCOGENICITY classes 
# ------------------------------------------------------------

# Group by oncogenicity and calculate key numbers 
stats_summary = variants.groupby("ONCOGENIC")["Germline_Proximity"].agg(
    count='count',
    median='median',
    mean='mean',
    std='std'
).reset_index()

print("Descriptive statistics:")
print(stats_summary)

# ------------------------------------------------------------
# Distribution of germline distance per class
# (excluding variants with missing distance)
# ------------------------------------------------------------

# filter to only include the wanted classes (oncogenic, likely oncogenic and likely neutral)
wanted_classes = ["Oncogenic", "Likely Oncogenic", "Likely Neutral"]
variants_plot = variants[variants["ONCOGENIC"].isin(wanted_classes)].copy() 

# remove rows with missing germline distances 
variants_plot = variants_plot.dropna(subset=["Germline_Proximity"])

# create log column (dist + 1 to avoid issues with 0 values) 
variants_plot["log_dist"] = np.log10(variants_plot["Germline_Proximity"] + 1)

# plot distributions 
sns.set_theme(style="whitegrid")
palette={
  "Oncogenic": "#C4473B",
  "Likely Oncogenic": "#D98C6A",
  "Likely Neutral": "#7e8aa2"
  }

# FIGURE 1: Density plot with all three classes 
print("Plotting distribution of germline proximity (comparison plot)...\n")

plt.figure(figsize=(8,5)) 
sns.kdeplot(
  data=variants_plot,
  x="log_dist",
  hue="ONCOGENIC",
  fill=True,
  common_norm=False,
  palette=palette, 
  alpha=0.5
)

plt.title("Distribution of Germline Proximity (Comparison)", fontsize=14)
plt.xlabel("Distance to nearest pathogenic germline variant (Log10 bp + 1)", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.savefig("plots/combined_dist.png", bbox_inches='tight')
plt.show()

print("Plotting complete! Plot saved as 'plots/combined_dist.png'.\n")


# FIGURE 2: Density plot per class 
print("Plotting germline distances per class...\n")

g = sns.FacetGrid(
    variants_plot, col="ONCOGENIC", hue="ONCOGENIC", 
    palette=palette, height=4, aspect=1.2, sharey=False
)
g.map(sns.histplot, "log_dist", kde=True, bins=30, alpha=0.4)

# Add vertical line for median of each plot 
for ax, label in zip(g.axes.flat, wanted_classes):
    median_val = variants_plot[variants_plot["ONCOGENIC"] == label]["log_dist"].median()
    ax.axvline(median_val, color="black", linestyle="--", alpha=0.7)
    ax.set_title(f"{label}\n(Median: 10^{median_val:.1f} bp)", fontsize=12)

g.set_axis_labels("Log10(Distance + 1)", "Count")
plt.tight_layout()
plt.savefig("plots/dists_per_class.png", bbox_inches='tight')
plt.show()

print("Plotting complete! Plot saved as 'plots/dists_per_class.png'.\n")


# FIGURE 3: Boxplot 

print("Plotting boxplot of germline distance data...\n")

plt.figure(figsize=(8,5))

sns.boxplot(
  data=variants_plot, 
  x="ONCOGENIC", 
  y="log_dist", 
  palette=palette)

plt.title("Boxplot of Distances (Log-scaled)")
plt.ylabel("Log10(Distance + 1)")

plt.tight_layout() 
plt.savefig("plots/boxplot_germline_dist.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Boxplot saved as 'plots/boxplot_germline_dist.png'.\n")

print("Germline proximity analysis complete!🎉🥳\n")