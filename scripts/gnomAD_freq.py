# ============================================================
# ANALYSIS: gnomAD allele frequencies by Oncogenicity category
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------
print("\nLoading variant data...\n")

variants = pd.read_csv(
    "annotation_pipeline/output/variants_with_func_sites.tsv",
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
    subset_af = subset[subset["gnomAD_AF"].notna()]
    missing_af = subset["gnomAD_AF"].isna().sum()

    print(f"{missing_af:,} of {total:,} '{status}' variants "
          f"({100 * missing_af / total:.1f}%) lack gnomAD allele frequency data.\n")

    # Summary statistics
    print(f"Summary statistics for gnomAD_AF among '{status}' variants with available data:\n")
    print(subset_af["gnomAD_AF"].describe())

    # Plot distribution
    print("\nPlotting gnomAD allele frequency distribution (log scale)...\n")
    plt.figure(figsize=(7, 4))
    sns.histplot(
        data=subset_af,
        x="gnomAD_AF",
        log_scale=True,
        bins=50,
        color=color
    )

    plt.axvline(0.001, color="red", linestyle="--", label="Rare/common cutoff (0.001)")
    plt.axvline(0.01, color="orange", linestyle="--", label="Polymorphism threshold (0.01)")
    plt.title(f"Distribution of gnomAD AF for '{status}' Variants", fontweight="bold")
    plt.xlabel("gnomAD_AF (log10 scale)")
    plt.ylabel("Number of variants")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"visualize_variants/plots/{plotname}", dpi=300, bbox_inches="tight")
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

# 1. Unknown
analyze_gnomad_af(variants, "Unknown", color="#848a8e", plotname="gnomAD_unknown.png")

# 2. Likely Oncogenic
analyze_gnomad_af(variants, "Likely Oncogenic", color="#d49392", plotname="gnomAD_likely_onco.png")

# 3. Oncogenic
analyze_gnomad_af(variants, "Oncogenic", color="#ef6f6c", plotname="gnomAD_onco.png")

# 4. Inconclusive
analyze_gnomad_af(variants, "Inconclusive", color="#f9c74f", plotname="gnomAD_inconclusive.png")

# 5. Likely Neutral
analyze_gnomad_af(variants, "Likely Neutral", color="#5b8fdc", plotname="gnomAD_likely_neutral.png")

print("gnomAD frequency analysis completed successfully for all oncogenicity classes.\n")


# ============================================================
# Log-scaled KDE comparison 
# ============================================================

# select classes 
wanted = ["Oncogenic", "Likely Neutral"]
filtered = variants[variants["ONCOGENIC"].isin(wanted)].copy()

# convert AF to numeric
filtered["gnomAD_AF"] = pd.to_numeric(filtered["gnomAD_AF"], errors="coerce")

# drop NA and zeros 
filtered = filtered.dropna(subset=["gnomAD_AF"])
filtered = filtered[filtered["gnomAD_AF"] > 0]

# define colors 
palette = {
    "Oncogenic": "#ef6f6c",
    "Likely Neutral": "#5b8fdc",
    }

# plot density 

plt.figure(figsize=(10,6)) 

ax = sns.kdeplot(
    data=filtered,
    x="gnomAD_AF",                                      
    hue="ONCOGENIC", 
    palette=palette,
    log_scale=True,
    linewidth=2,
    common_norm=False
)

plt.title("gnomAD AF Distribution Across Oncogenicity Classes", 
          fontsize=16, fontweight="bold") 
plt.xlabel("gnomAD_AF (log10 scale)")
plt.ylabel("Density")
plt.tight_layout()

plt.savefig("visualize_variants/plots/gnomAD_combined_KDE.png", dpi=300, bbox_inches="tight") 

plt.show() 
