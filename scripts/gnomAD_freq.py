"""
# ============================================================
# ANALYSIS: gnomAD allele frequencies by Oncogenicity category
# ============================================================

Script: gnomAD_freq.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across gnomAD allele frequencies. 

Script content:
--------------
1. Function to analyze and plot gnomAD allele frequencies for a given oncogenicity class. 
2. gnomAD AF analysis for all oncogenicity classes 
3. Log-scaled KDE-comparison between oncogenic and likely neutral variants. 
4. Statistical analysis: Mann-Whitney U test
5. Statistical analysis: Chi-Squared test 

All plots are saved in:
    plots/

"""

print("\n========================================================")
print("gnomAD ALLELE FREQUENCY ANALYSIS")
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

print("\n------------------------------------------------------")
print("LOAD VARIANT DATA")
print("------------------------------------------------------\n")

print("Loading variant data...\n")

variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.\n")


# ============================================================
# Function to analyze gnomAD allele frequencies
# ============================================================


def analyze_gnomad_af(variants: pd.DataFrame, status: str, plotname: str, color: str = "teal"):
    """Analyze and plot gnomAD_AF distribution for a given ONCOGENIC category."""

    print(f"Extracting gnomAD allele frequencies for variants with '{status}' oncogenicity...\n")

    subset = variants[variants["ONCOGENIC"] == status].copy()
    total = len(subset)
    print(f"Found {total:,} variants with '{status}' oncogenicity.\n")

    # Convert to numeric and separate valid frequencies
    subset["gnomAD_AF"] = pd.to_numeric(subset["gnomAD_AF"], errors="coerce")
    subset_af = subset[subset["gnomAD_AF"].notna() & subset["gnomAD_AF"] > 0]
    missing_af = subset["gnomAD_AF"].isna().sum()

    print(f"{missing_af:,} of {total:,} '{status}' variants "
          f"({100 * missing_af / total:.1f}%) lack gnomAD allele frequency data.\n")

    # Summary statistics
    print(f"Summary statistics for gnomAD_AF among '{status}' variants with available data:\n")
    print(subset_af["gnomAD_AF"].describe())

    # Plot distribution
    print("\nPlotting gnomAD allele frequency distribution (log scale)...\n")
    plt.figure(figsize=(8,5))

    sns.histplot(
        data=subset_af,
        x="gnomAD_AF",
        log_scale=True,
        bins=50,
        color=color, 
        edgecolor="0.1",
        linewidth=0.3
    )

    plt.axvline(0.001, color="red", linestyle="--", label="Rare/common cutoff (0.001)")
    plt.axvline(0.01, color="orange", linestyle="--", label="Polymorphism threshold (0.01)")
    plt.title(f"Distribution of gnomAD AF for '{status}' Variants", fontsize=14, pad=10)
    plt.xlabel("gnomAD_AF (log10 scale)", fontsize=12)
    plt.ylabel("Number of variants", fontsize=12)
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(f"plots/{plotname}", dpi=300, bbox_inches="tight")
    plt.show()
  

    # Count rare vs common
    common = (subset_af["gnomAD_AF"] > 0.01).sum()
    rare = (subset_af["gnomAD_AF"] <= 0.01).sum()

    print(f"Common '{status}' variants (gnomAD_AF > 0.01): {common:,}")
    print(f"Rare   '{status}' variants (gnomAD_AF ≤ 0.01): {rare:,}")
    print(f"Total with AF available: {len(subset_af):,}\n")
    print("-" * 60 + "\n")


# ============================================================
# Run analysis for all oncogenic categories 
# ============================================================

print("\n------------------------------------------------------")
print("CALL gnomAD AF ANALYSIS FUNCTION")
print("------------------------------------------------------\n")

print("Running gnomAD AF analysis for all given oncogenicity classes...\n")
print("-" * 60 + "\n")

# 1. Unknown
analyze_gnomad_af(variants, "Unknown", color="#848a8e", plotname="gnomAD_unknown.png")

# 2. Likely Oncogenic
analyze_gnomad_af(variants, "Likely Oncogenic", color="#D98C6A", plotname="gnomAD_likely_onco.png")

# 3. Oncogenic
analyze_gnomad_af(variants, "Oncogenic", color="#C4473B", plotname="gnomAD_onco.png")

# 4. Inconclusive
analyze_gnomad_af(variants, "Inconclusive", color="#f9c74f", plotname="gnomAD_inconclusive.png")

# 5. Likely Neutral
analyze_gnomad_af(variants, "Likely Neutral", color="#7e8aa2", plotname="gnomAD_likely_neutral.png")

print("gnomAD frequency analysis completed successfully for all oncogenicity classes.")
print("Plots saved in folder 'plots/'")


# ============================================================
# Log-scaled KDE comparison 
# ============================================================

print("\n------------------------------------------------------")
print("LOG-SCALED KDE COMPARISON (ONCOGENIC VS LIKELY NEUTRAL)")
print("------------------------------------------------------\n")

# select classes 
wanted = ["Oncogenic", "Likely Neutral"]
filtered = variants[variants["ONCOGENIC"].isin(wanted)].copy()

# convert gnomAD_AF to numeric
filtered["gnomAD_AF"] = pd.to_numeric(filtered["gnomAD_AF"], errors="coerce")

# drop NA and zero values
filtered = filtered[filtered["gnomAD_AF"].notna() & (filtered["gnomAD_AF"] > 0)]

# define colors for the given classes
palette = {
    "Oncogenic": "#C4473B",
    "Likely Neutral": "#7e8aa2",
    }

# create log-scaled density plot 
print("Plotting density of gnomAD allele frequencies for oncogenic vs likely neutral...\n")

plt.figure(figsize=(8,5)) 

sns.kdeplot(
    data=filtered,
    x="gnomAD_AF",                                      
    hue="ONCOGENIC", 
    palette=palette,
    log_scale=True,
    linewidth=1.5,
    common_norm=False
)

plt.title("gnomAD AF Distribution Across Oncogenicity Classes", fontsize=14, pad=10) 
plt.xlabel("gnomAD_AF (log10 scale)", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.tight_layout()

plt.savefig("plots/gnomAD_combined_KDE.png", dpi=300, bbox_inches="tight") 

plt.show() 

print("Plotting complete! Plot saved as 'plots/gnomAD_combined_KDE.png'\n")


# ============================================================
# Statistical analyses: Mann-Whitney U test 
# ============================================================

# Hypotheses
# H0 The distribution of gnomAD allele frequencies is the same for oncogenic and likely neutral variants.
# H1 The distribution of gnomAD allele frequencies is not the same for oncogenic and likely neutral variants.

# Model assumptions: 
# The variable (gnomAF) is continuous 
# The data is assumed to be non-normal
# The data in both groups have similar distributions 
# The samples should be independent 

from scipy.stats import mannwhitneyu

# define the data 
oncogenic = filtered[filtered["ONCOGENIC"] == "Oncogenic"]["gnomAD_AF"]
neutral = filtered[filtered["ONCOGENIC"]=="Likely Neutral"]["gnomAD_AF"]

# perform Mann-Whitney U test 
stat, p = mannwhitneyu(oncogenic, neutral, alternative="two-sided") 
print(f"Mann-Whitney U: {stat:.3f}, p-value: {p:.4f}")

# calculate rank-biserial correlation 
# (effect size for mann-whitney u) 
n1 = len(oncogenic)
n2 = len(neutral)
r = 1 - (2 * stat) / (n1 * n2)
print(f"Rank-biserial correlation: {r:.3f}")

probability = (1+r)/2 
print(f"The probability of a random oncogene variant having a higher gnomAD_AF than a neutral is: {probability*100:.2f}%.\n")

# ============================================================
# Statistical analysis: Chi-Square
# ============================================================

# Check whether there is an significant association between variant oncogenicity and  
# the precence of gnomAD allele frequency data  

# Model assumptions: 
# Both variables are categorical.
# All observations are independent.
# Each observation must only contribute to one cell.
# The expected frequency in each cell should be at least five.

# Hypotheses: 
# H0:   There is no association between variant pathogenicity (oncogenic vs. likely neutral) 
#       and the presence of gnomAD allele frequency data
# H1:   There is a statistically significant association between variant pathogenicity and 
#       the presence of gnomAD allele frequency data.

from scipy.stats import chi2_contingency 
from scipy.stats.contingency import odds_ratio

# make contingency table for observed values 
oncogenic = variants[variants["ONCOGENIC"] == "Oncogenic"]["gnomAD_AF"]
neutral = variants[variants["ONCOGENIC"]=="Likely Neutral"]["gnomAD_AF"]

oncogenic_missing = (pd.to_numeric(oncogenic, errors="coerce").isna()).sum() 
neutral_missing = (pd.to_numeric(neutral, errors="coerce").isna()).sum() 

table = [
    [len(oncogenic) - oncogenic_missing, oncogenic_missing],
    [len(neutral) - neutral_missing, neutral_missing]
]

observed_table = pd.DataFrame(
    table, 
    columns=["Reported (AF>0)", "Missing (NaN)"],
    index=["Oncogenic", "Likely Neutral"]
)

print("Contingency table (observed values):")
print(observed_table) 
print("-"*30)

print("Running Chi-square test on gnomAD_AF variant data...\n")

# run chi-square test 
chi2, p, dof, expected = chi2_contingency(observed_table)

# contingency table expected values (under H0) 
expected_table = pd.DataFrame(
    expected,
    columns=["Reported (AF>0)", "Missing (NaN)"],
    index=["Oncogenic", "Likely Neutral"]
)

print("Contingency table (expected values):")
print(expected_table.round(2))
print("-"*30)

# print results from Chi-square 
print("Results:")
print(f"Chi-square: {chi2:.3f}, p-value: {p:.4f}")

# cramers v (effect size for Chi-square) 
n = observed_table.values.sum()
cramers_v = np.sqrt(chi2 / (n * (min(2, 2) - 1)))
print(f"Cramér's V (effect size Chi-square): {cramers_v:.3f}\n")
print("-"*30)

# calculate the odds-ratio 
result = odds_ratio(observed_table)
ci = result.confidence_interval(confidence_level=0.95)

print(f"Odds-ratio: {result.statistic:.2f}")
print(f"95% CI: [{ci.low:.2f}, {ci.high:.2f}]")

print("\nVariant gnomAD allele frequency analysis complete!🥳🥳\n")

