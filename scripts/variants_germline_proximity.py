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
7. Statistical test: Mann-Whitney U and rank-biserion correlation 

All plots are saved in:
   plots/

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

print("Descriptive statistics germline distances:")
print(stats_summary)

# ------------------------------------------------------------
# Distribution of germline distance per class
# (excluding variants with missing distance)
# ------------------------------------------------------------

# filter to only include the wanted classes (oncogenic, likely oncogenic and likely neutral)
wanted_classes = ["Oncogenic", "Likely Neutral"]
variants_plot = variants[variants["ONCOGENIC"].isin(wanted_classes)].copy() 

# remove rows with missing germline distances 
variants_plot = variants_plot.dropna(subset=["Germline_Proximity"])

# create log column (dist + 1 to avoid issues with 0 values) 
variants_plot["log_dist"] = np.log10(variants_plot["Germline_Proximity"] + 1)

# plot distributions 
sns.set_theme(style="whitegrid")
palette={
  "Oncogenic": "#C4473B",
  "Likely Neutral": "#7e8aa2"
  }

# FIGURE 1: Density plot with both classes
print("Plotting distribution of germline proximity (comparison plot)..\n")

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


# FIGURE 2: Histogram with distance proportions per class 
print("-"*30)
print("Plotting germline distances per class..\n")

g = sns.FacetGrid(
    variants_plot, col="ONCOGENIC", hue="ONCOGENIC", 
    palette=palette, height=4, aspect=1.2, sharey=True
)

g.map(sns.histplot, "log_dist", kde=True, bins=30, alpha=0.4, stat="proportion")

# Add vertical line for median of each plot 
for ax, label in zip(g.axes.flat, wanted_classes):
    subset = variants_plot[variants_plot["ONCOGENIC"] == label]
    median_val = subset["log_dist"].median()
    n_count = len(subset)

    ax.axvline(median_val, color="black", linestyle="--", alpha=0.7)
    ax.set_title(f"{label} (n={n_count})\n(Median: 10^{median_val:.1f} bp)", fontsize=12)

    g.set_axis_labels("Log10(Distance + 1)", "Proportion")

plt.tight_layout()
plt.savefig("plots/dists_per_class.png", bbox_inches='tight')
plt.show()

print("Plotting complete! Plot saved as 'plots/dists_per_class.png'.\n")


# FIGURE 3: Boxplot 
print("-"*30)
print("Plotting boxplot of germline distance data..\n")

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

# ============================================================
# Statistical test: Mann-Whitney U test 
# ============================================================

# Statistical test to check if there is a significant difference in the distribution 
# of germline distances between oncogenic and likely neutral variants. 

# Hypotheses
#       H0 The distribution of distances to the nearest germline variant is 
#       the same for oncogenic and likely neutral variants.
#       H1 The distribution of distances to the nearest germline variant is 
#       not the same for oncogenic and likely neutral variants.

# Model assumptions: 
# 1.  The variable (germline proximity) is continuous 
# 2.  The data is assumed to be non-normal
# 3.  The data in both groups have similar distributions 
# 4.  The samples should be independent 

# import library 
from scipy.stats import mannwhitneyu

# define the data, drop NA values 
oncogenic = variants[variants["ONCOGENIC"] == "Oncogenic"]["Germline_Proximity"].dropna()
neutral = variants[variants["ONCOGENIC"] == "Likely Neutral"]["Germline_Proximity"].dropna()

# perform Mann-Whitney U test 
print("-"*30)
print("Running Mann-Whitney U test on the germline proximity data..\n")

alpha = 0.05 

stat, p = mannwhitneyu(oncogenic, neutral, alternative="two-sided") 
print("Results:")
print(f"Mann-Whitney U: {stat:.3f}, p-value: {p:.4f}")

if p < alpha:
    print("Reject the null hypothesis. The germline distances for oncogenic and likely neutral variants come from different distributions.\n")
else: 
    print("Failed to reject the null hypothesis.\n") 

# calculate rank-biserial correlation 
# (effect size for mann-whitney u) 
n1 = len(oncogenic)
n2 = len(neutral)
r = (2 * stat) / (n1 * n2) - 1
print(f"Rank-biserial correlation: {r:.3f}")

# calculate probability 
probability = (1+r)/2 
print(f"The probability of a random oncogenic variant having a higher germline distance than a neutral variant is: {probability*100:.2f}%.")

print("-"*30)


# ============================================================
# Statistical test: Kolmogorov-Smirnov 
# ============================================================

# Model assumptions:
#       The samples should be independent.       
#       The dependent variable must be measured on an ordinal or continuous scale. 
#       The distributions should be continuous to avoid ties. 

# Hypotheses: 
#       H0: The samples come from the same distribution.
#       H1: The samples come from different distributions. 


# import library
from scipy.stats import ks_2samp

# run statistic
print("Running Kolmogorov-Smirnov test on germline distance data..\n")

ks_statistic, p_value = ks_2samp(oncogenic, neutral, alternative='two-sided', mode='auto')

alpha = 0.05

print("Results KS-test:")
print(f"KS-statistic: {ks_statistic}")
print(f"p-value: {p_value:.4f}")

if p_value < alpha:
    print("Reject the null hypothesis. The oncogenic and neutral variants come from two different distributions.")
else:
    print("Failed to reject the null hypothesis.") 


# ============================================================
# Investigate the bi-modal distribution for
# germline distances among oncogenic variants 
# ============================================================

# Define the two peaks: 
# peak A: log distance < 1.0 (< 10bp) 
# peak B: log distance > 1.5 (> 30bp) 

# Define the data 
oncogenic_var  = variants_plot[variants_plot["ONCOGENIC"] == "Oncogenic"].copy() 
peak_a_variants = oncogenic_var[oncogenic_var["log_dist"] <= 1]
peak_b_variants = oncogenic_var[oncogenic_var["log_dist"] >= 1.5]

# Extract top genes (based on the number of variants)
top_genes_full = oncogenic_var["Hugo_Symbol"].value_counts().head(10)
genes_peak_a   = peak_a_variants["Hugo_Symbol"].value_counts().head(10)
genes_peak_b   = peak_b_variants["Hugo_Symbol"].value_counts().head(10)

# Print results
print("The top 10 oncogenic genes total are:")
print(top_genes_full)
print("The top 10 oncogenic genes in peak A are:")
print(genes_peak_a) 
print("The top 10 oncogenic genes i peak B are:")
print(genes_peak_b)


# create proportion plot 

# label the peaks
peak_a_variants["Peak"] = "Peak A (<10bp)"
peak_b_variants["Peak"] = "Peak B (>30bp)"

# combine df
combined_df = pd.concat([peak_a_variants, peak_b_variants])

# only include genes from peak a and b 
genes_a_b = list(set(genes_peak_a.index) | set(genes_peak_b.index))
combined_df = combined_df[combined_df["Hugo_Symbol"].isin(genes_a_b)]

# Count variants per gene per peak
peak_counts = (
    combined_df
    .groupby(["Hugo_Symbol", "Peak"])
    .size()
    .unstack(fill_value=0)  # long format to pivot 
)

# Calculate proportion of a gene's variants in each peak 
peak_proportions = peak_counts.div(peak_counts.sum(axis=1),axis=0)

# Proportion plot 
peak_proportions.plot(
    kind="bar", 
    stacked=True, 
    figsize=(10,5),
    color=("#d6b4be", "#d16a58")
    )

plt.title("Distribution of oncogenic variants across distance peaks (top genes)", fontsize=14)
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Proportion of variants", fontsize=12)
plt.xticks(rotation=45, fontsize=9)
plt.legend(
    title="Peak",
    loc="upper right",
)
plt.tight_layout()
plt.savefig("plots/germline_proportions.png", dpi=300, bbox_inches="tight")
plt.show()


print("\nGermline proximity analysis complete!🥳🥳\n")