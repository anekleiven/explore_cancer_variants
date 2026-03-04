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
4. Statistical analysis: Mann-Whitney U test with rank-biserion correlation 
5. Statistical analysis: Chi-Squared test with Cramer's V and Odds-Ratio 

All plots are saved in:
    plots/

"""
print("-"*50)
print("gnomAD Allele Frequency Analysis🤓")
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

print(f"Loaded {len(variants):,} variants.")


# ============================================================
# Function to analyze gnomAD allele frequencies
# ============================================================

def analyze_gnomad_af(variants: pd.DataFrame, status: str, plotname: str, color: str = "teal"):
    """Analyze and plot gnomAD_AF distribution for a given ONCOGENIC category."""

    print(f"Extracting gnomAD allele frequencies for variants with '{status}' oncogenicity..")

    subset = variants[variants["ONCOGENIC"] == status].copy()
    total = len(subset)
    print(f"Found {total:,} variants with '{status}' oncogenicity.\n")

    # Convert to numeric and separate valid frequencies
    subset["gnomAD_AF"] = pd.to_numeric(subset["gnomAD_AF"], errors="coerce")
    subset_af = subset[subset["gnomAD_AF"].notna() & subset["gnomAD_AF"] > 0]
    missing_af = subset["gnomAD_AF"].isna().sum()

    print(f"{missing_af:,} of {total:,} '{status}' variants "
          f"({100 * missing_af / total:.1f}%) lack gnomAD allele frequency data.")

    # Plot distribution
    print("\nPlotting gnomAD allele frequency distribution (log-scaled)..\n")
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

    print("\nDistribution of common vs. rare variants:")
    print(f"Common '{status}' variants (gnomAD_AF > 0.01): {common:,}")
    print(f"Rare   '{status}' variants (gnomAD_AF ≤ 0.01): {rare:,}")
    print(f"Total with AF available: {len(subset_af):,}")
    print("-"*30)

# ============================================================
# Run analysis for all oncogenic categories 
# ============================================================

print("-"*30)
print("Running gnomAD AF analysis for all given oncogenicity classes..\n")

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

print("-"*30)
print("Log-scaled KDE comparison plot (oncogenic vs likely neutral)")


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
print("Plotting density of gnomAD allele frequencies for oncogenic vs likely neutral..\n")

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
# Statistical test: Mann-Whitney U test 
# ============================================================

# Hypotheses
# H0 The distribution of gnomAD allele frequencies is the same for oncogenic and likely neutral variants.
# H1 The distribution of gnomAD allele frequencies is not the same for oncogenic and likely neutral variants.

# Model assumptions: 
# 1.  The variable (gnomAF) is continuous 
# 2.  The data is assumed to be non-normal
# 3.  The data in both groups have similar distributions 
# 4.  The samples should be independent 

# import library 
from scipy.stats import mannwhitneyu

print("-"*30)
print("Running Mann-Whitney U test on gnomAD_AF variant data..\n")

# Group by oncogenicity and calculate descriptive statistics 
stats_summary = filtered.groupby("ONCOGENIC")["gnomAD_AF"].agg(
    count='count',
    median='median',
    mean='mean',
    std='std'
).reset_index()

print("Decriptive statistics:")
print(stats_summary,"\n")

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
r = (2 * stat) / (n1 * n2) - 1 
print(f"Rank-biserial correlation: {r:.3f}")

# calculate probability 
probability = (1+r)/2 
print(f"The probability of a random oncogenic variant having a higher gnomAD_AF than a neutral is: {probability*100:.2f}%.")
print("-"*30)

# ============================================================
# Statistical test: Chi-Square
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

print("Running Chi-square test on gnomAD_AF variant data..\n")

# import libraries 
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


# run chi-square test 
chi2, p, dof, expected = chi2_contingency(observed_table)

# contingency table expected values (under H0) 
expected_table = pd.DataFrame(
    expected,
    columns=["Reported (AF>0)", "Missing (NaN)"],
    index=["Oncogenic", "Likely Neutral"]
)

print("\nContingency table (expected values):")
print(expected_table.round(2))

# print results from Chi-square 
print("\nResults:")
print(f"Chi-square: {chi2:.3f}, p-value: {p:.4f}")

# cramers v (effect size for Chi-square) 
n = observed_table.values.sum()
k = min(observed_table.shape) - 1 
cramers_v = np.sqrt(chi2 / (n*k))
print(f"Cramér's V (effect size Chi-square): {cramers_v:.3f}")

# calculate the odds-ratio 
result = odds_ratio(observed_table)
ci = result.confidence_interval(confidence_level=0.95)

print(f"Odds-ratio: {result.statistic:.2f}")
print(f"95% CI: [{ci.low:.2f}, {ci.high:.2f}]")
print("-"*30)

print("\nVariant gnomAD allele frequency analysis complete!🥳🥳\n")

