# ===============================================================
# gnomAD Allele Frequency Analysis
# ===============================================================

"""
Script: gnomAD_freq.py
Author: Ane Kleiven

This script performs a multi-step explorative data analysis on 
gnomAD allele frequencies. 

Script content:
--------------
1. Load variant data 
2. Filter and clean variant data 
3. Perform descriptive statistics per oncogenicity class 
4. gnomAD AF analysis per oncogenicity class
5. Log-scaled KDE-comparison between oncogenic and likely neutral variants. 


All plots are saved in:
    plots/gnomad

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
import argparse
from pathlib import Path
import os

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Explore gnomAD AF in variant data."
    ) 

    parser.add_argument(
        "--variants", 
        type=Path, 
        required=False, 
        default="/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_tsg_og.tsv",
        help="Path to the input file with variant data."
    )

    return parser.parse_args() 


args = getargs() 

# ------------------------------------------------------------
# Create output directory 
# ------------------------------------------------------------

save_dir = "plots/gnomad"
os.makedirs(save_dir, exist_ok=True) 

# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------

print(f"Loading variant file:\n{args.variants}")

variants = pd.read_csv(
    args.variants,
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.\n")

# ------------------------------------------------------------
# Filter and clean data 
# ------------------------------------------------------------

print("Filter and clean gnomAD data..")

# Select classes 
wanted = ["Oncogenic", "Likely Neutral"]
filtered = variants[variants["ONCOGENIC"].isin(wanted)].copy()

# Convert gnomAD_AF to numeric
filtered["gnomAD_AF"] = pd.to_numeric(filtered["gnomAD_AF"], errors="coerce")

# Drop NA and zero values
filtered = filtered[filtered["gnomAD_AF"].notna() & (filtered["gnomAD_AF"] > 0)]

print(f"{len(filtered)} oncogenic and neutral variants have gnomAD allele frequencies.")
print("-"*50)

# ------------------------------------------------------------
# Descriptive statistics per oncogenicity class
# ------------------------------------------------------------

# Group by oncogenicity and calculate descriptive statistics 
stats_summary = filtered.groupby("ONCOGENIC")["gnomAD_AF"].agg(
    count='count',
    median='median',
    mean='mean',
    std='std'
).reset_index()

print("Decriptive statistics by oncogenicity:")
print(stats_summary)
print("-"*50)

# ------------------------------------------------------------
# Function to analyze gnomAD allele frequencies
# ------------------------------------------------------------

def analyze_gnomad_af(df: pd.DataFrame, status: str, plotname: str, color: str = "teal"):
    """
    Analyze and plot gnomAD_AF for a given oncogencity class.
    Log-scaled distribution. 

    """
    print(f"Extracting gnomAD allele frequencies for variants with '{status}' oncogenicity..")

    subset = df[df["ONCOGENIC"] == status].copy()
    print(f"Found {len(subset):,} variants with '{status}' oncogenicity.")

    # Plot distribution
    print("Plotting gnomAD allele frequency distribution (log-scaled)..\n")
    plt.figure(figsize=(8,5))
    sns.set_style("whitegrid")

    sns.histplot(
        data=subset,
        x="gnomAD_AF",
        log_scale=True,
        bins=25,
        color=color, 
        edgecolor="0.1",
        linewidth=0.3
    )

    plt.axvline(0.01, color="red", linestyle="--", label="Rare/common cutoff (0.01)")
    plt.axvline(0.1, color="orange", linestyle="--", label="Polymorphism threshold (0.1)")
    plt.title(f"gnomAD AF for '{status}' Variants", fontsize=14, pad=10)
    plt.xlabel("gnomAD_AF (log10 scale)", fontsize=12)
    plt.ylabel("Number of variants", fontsize=12)
    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(f"{save_dir}/{plotname}", dpi=300, bbox_inches="tight")
    plt.show()
  
    # Count rare vs common
    common = (subset["gnomAD_AF"] > 0.01).sum()
    rare = (subset["gnomAD_AF"] <= 0.01).sum()

    print("\nDistribution of common vs. rare variants:")
    print(f"Common '{status}' variants (gnomAD_AF > 0.01): {common:,}")
    print(f"Rare   '{status}' variants (gnomAD_AF ≤ 0.01): {rare:,}")
    print(f"Total with AF available: {len(subset):,}")
    print("-"*50)

# ------------------------------------------------------------
# Run analysis for all oncogenic classes
# ------------------------------------------------------------

print("Running gnomAD AF analysis for all given oncogenicity classes..\n")

# Oncogenic
analyze_gnomad_af(filtered, "Oncogenic", color="#c4314a", plotname="gnomAD_onco.png")

# Likely Neutral
analyze_gnomad_af(filtered, "Likely Neutral", color="#88aed1", plotname="gnomAD_likely_neutral.png")

print("gnomAD frequency analysis successfully completed for all oncogenicity classes.")
print(f"Plots saved in folder '{save_dir}'")


# ------------------------------------------------------------
# Log-scaled KDE comparison 
# ------------------------------------------------------------

print("-"*50)
print("Log-scaled KDE comparison plot (Oncogenic vs. Likely Neutral)")

# Define colors for the given classes
palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
    }

# Create log-scaled density plot 
print("Plotting density of gnomAD allele frequencies for Oncogenic vs. Likely Neutral..\n")

plt.figure(figsize=(8,5)) 
sns.set_style("white")

sns.kdeplot(
    data=filtered,
    x="gnomAD_AF",                                      
    hue="ONCOGENIC", 
    palette=palette,
    log_scale=True,
    linewidth=1.5,
    common_norm=False
)

sns.move_legend(plt.gca(), "upper left", bbox_to_anchor=(1.05, 1), title="Oncogenicity Class", frameon=False)
sns.despine()

plt.title("gnomAD Population Allele Frequencies", fontsize=14, pad=15) 
plt.xlabel("Allele frequency (log10 scale)", fontsize=12, labelpad=10)
plt.ylabel("Density", fontsize=12, labelpad=10)

plt.tight_layout()
plt.savefig(f"{save_dir}/gnomAD_combined_KDE.png", dpi=300, bbox_inches="tight") 
plt.show() 

print(f"Plotting complete! Plot saved as '{save_dir}/gnomAD_combined_KDE.png'")
print("-"*50)


print("\nVariant gnomAD allele frequency analysis complete!🥳🥳\n")
