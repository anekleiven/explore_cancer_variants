"""
Functional Sites Analysis 
-----------------------------------------------

Script: functional_sites.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across different functional protein features.

Major outputs:
--------------
1. Counts and fractions of Oncogenic vs Likely Neutral variants inside/outside functional sites
2. Plot: Variant counts per functional site type per class
3. Plot: Fraction of variants per functional site type per class 
4. Plot: Top genes by oncogenic variant burden per functional site
5. Plot: Heatmap of top oncogenic genes (variant burden) in filtered functional sites  
6. Plot: Oncogenic-to-neutral ratios in top oncogenic genes per functional site type

All plots are saved in:
    plots/functionalsites

"""

print("-"*50)
print("Variant functional site analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
from pathlib import Path

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Explore functional sites in variant data."
    ) 

    parser.add_argument(
        "--variants", 
        type=Path, 
        required=False, 
        default="/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
        help="Path to the input file with variant data."
    )

    return parser.parse_args() 


args = getargs() 

# ------------------------------------------------------------
# Import Variant Data
# ------------------------------------------------------------

print("Loading variant data..")

variants = pd.read_csv(
    args.variants, sep="\t", low_memory=False
)

print(f"Loaded {len(variants):,} variants.")

# ------------------------------------------------------------
# Filter to Oncogenic and Likely Neutral variants
# ------------------------------------------------------------

print("-"*50)
print("Filtering data oncogenic & likely oncogenic variants..")

classes = ["Oncogenic", "Likely Neutral"]
variants = variants[variants["ONCOGENIC"].isin(classes)]
print(f"Remaining variants after filtering: {len(variants):,}")

# Ensure boolean data type for 'IN_FUNC_SITE' 
variants["IN_FUNC_SITE"] = variants["IN_FUNC_SITE"].astype(bool)

# ------------------------------------------------------------
# Summary: variants inside vs outside functional sites
# ------------------------------------------------------------

print("-"*50)
print("Creating summary of variants inside/outside functional sites..")
summary = []
for c in classes:
    subset = variants[variants["ONCOGENIC"] == c]
    total = len(subset)
    inside = subset["IN_FUNC_SITE"].sum()
    fraction_in = inside / total if total > 0 else 0

    summary.append({
        "Class": c,
        "Total_variants": total,
        "Inside_func_site": inside,
        "Outside_func_site": total - inside,
        "Fraction_inside": fraction_in
    })

summary_df = pd.DataFrame(summary)
print("Distribution of variants inside and outside functional sites:")
print(summary_df)

# ------------------------------------------------------------
# Expand FEATURE_TYPE 
# ------------------------------------------------------------

print("-"*50)
print("Expanding FEATURE_TYPE (one row per variant per functional site)..")

expanded = (
    variants
    .dropna(subset=["FEATURE_TYPE"])
    .assign(FEATURE_TYPE=lambda df: df["FEATURE_TYPE"].str.split(";"))
    .explode("FEATURE_TYPE")
)

expanded["FEATURE_TYPE"] = expanded["FEATURE_TYPE"].str.strip()
print(f"\nExpanded to {len(expanded):,} feature-variant rows.")

# ------------------------------------------------------------
# Count variants per functional site type per class
# ------------------------------------------------------------

print("-"*50)
print("Counting variants per functional site per class..")

counts = (
    expanded
    .groupby(["FEATURE_TYPE", "ONCOGENIC"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Number of variants in each functional site for all classes:")
print(counts)

# ------------------------------------------------------------
# Plot raw counts
# ------------------------------------------------------------

print("-"*50)
print("Plotting raw counts of variants in functional sites..")

palette = {"Oncogenic": "#c4314a", "Likely Neutral": "#88aed1"}

plt.figure(figsize=(8,5))
sns.barplot(
    data=counts,
    x="FEATURE_TYPE",
    y="Variant_Count",
    hue="ONCOGENIC",
    palette=palette, 
    edgecolor="0.1",
    linewidth=0.3
)
plt.title("Variant Counts per Functional Site Type", fontsize=14, pad=10)
plt.xlabel("Functional Site", fontsize=12)
plt.ylabel("Variant Count", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("plots/functionalsites/counts_per_site.png", dpi=300)
plt.show()

print(f"Plotting complete. Saved as 'plots/functionalsites/counts_per_site.png'")

# ------------------------------------------------------------
# Compute fractions
# ------------------------------------------------------------

print("-"*50)
print("Computing fractions of variants in the different functional sites..")

totals = variants["ONCOGENIC"].value_counts().rename("Total")

counts = counts.merge(totals, left_on="ONCOGENIC", right_index=True)
counts["Fraction"] = counts["Variant_Count"] / counts["Total"]

print("\nFraction of variants in each functional site type per class:")
print(counts)

# ------------------------------------------------------------
# Plot fractions
# ------------------------------------------------------------

print("-"*50)
print("Plotting fractions of variants for each functional site type..\n")

plt.figure(figsize=(8,5))
sns.barplot(
    data=counts,
    x="FEATURE_TYPE",
    y="Fraction",
    hue="ONCOGENIC",
    palette=palette, 
    edgecolor="0.1",
    linewidth=0.3
)
plt.title("Variant Fractions per Functional Site Type", fontsize=14, pad=10)
plt.xlabel("Functional Site", fontsize=12)
plt.ylabel("Fraction", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("plots/functionalsites/fraction_per_site.png", dpi=300)
plt.show()

print("Plotting complete. Saved as 'plots/functionalsites/fraction_per_site.png'")

# ------------------------------------------------------------
# Top genes by oncogenic variant contribution 
# per functional site 
# ------------------------------------------------------------

print("-"*50)
print("Computing oncogenic gene contributions per functional site..")

# Oncogenic counts
oncogenic = expanded[expanded["ONCOGENIC"] == "Oncogenic"].copy()

# Number of oncogenic variants per gene and feature type
onco_counts = (
    oncogenic.groupby(["Hugo_Symbol", "FEATURE_TYPE"])
    .size()
    .reset_index(name="Count"))

feature_totals = onco_counts.groupby("FEATURE_TYPE")["Count"].sum()

# Calculate contribution fraction: 
# "How much of this sites total oncogenic burden comes from this gene?"
onco_counts["Contribution_Fraction"] = (
    onco_counts["Count"] / onco_counts["FEATURE_TYPE"].map(feature_totals)
)

print("Example output oncogenic contributions table:")
print(onco_counts.head(5))

# ------------------------------------------------------------
# Plot oncogenic variant contribution by gene
# per functional site 
# ------------------------------------------------------------

print("-"*50)
print("Plotting oncogenic contribution fractions by gene per functional site..")

# Minimum number of variants for a site to be plotted 
minimum_var = 20 

top_onco_genes = set()  

# Loop over all feature types
for ft in onco_counts["FEATURE_TYPE"].unique():

    total_variants = onco_counts[onco_counts["FEATURE_TYPE"] == ft]["Count"].sum()

    # Skip features with too few variants
    if total_variants < minimum_var:
        print(f"Skipping {ft}. Only {total_variants} variants (< {minimum_var}).")
        continue 

    subset = (
        onco_counts[onco_counts["FEATURE_TYPE"] == ft]
        .sort_values("Contribution_Fraction", ascending=False)
        .head(20)
    )

    top_onco_genes.update(subset["Hugo_Symbol"].tolist())
    
    fig, ax = plt.subplots(figsize=(8,5)) 

    ax.bar(subset["Hugo_Symbol"],
           subset["Contribution_Fraction"], 
           color="#c4314a", 
           edgecolor="0.1",
           linewidth=0.3)
    
    ax.set_title(f"Top genes by oncogenic variant burden — {ft}", fontsize=14, pad=10)
    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Fraction of oncogenic variants in site", fontsize=12)
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    # Add total variants in plot 
    ax.annotate(f"n = {total_variants} oncogenic variants",
                xy=(0.99,0.97), xycoords="axes fraction",
                ha="right", va="top", fontsize=9,
                color="0.4")

    plt.tight_layout()
    plt.savefig(f"plots/functionalsites/contributing_genes_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved in 'plots/functionalsites/'")

# ------------------------------------------------------------
# Heatmap of Oncogenic Driver Genes Across Functional Sites
# ------------------------------------------------------------

print("-"*50)
print("Creating heatmap of top contributing genes per functional site type..")

# filter to functional sites with the most variants 
filtered_sites = ["Binding site", "Modified residue", "Region", "Topological domain"]

# find the top 5 genes with the highest amount of oncogenic variants per site (fraction) 
top_genes = (
    onco_counts[onco_counts["FEATURE_TYPE"].isin(filtered_sites)]
    .sort_values(['FEATURE_TYPE', 'Contribution_Fraction'], ascending=[True,False])
    .groupby('FEATURE_TYPE')
    .head(5)
)

# create matrix for heatmap 
pivot_heatmap = top_genes.pivot(
    index='Hugo_Symbol', 
    columns='FEATURE_TYPE', 
    values='Contribution_Fraction'
).fillna(0)

plt.figure(figsize=(8,6))
sns.heatmap(pivot_heatmap, annot=True, cmap="Reds", fmt=".2f")
plt.title("Gene Contribution per Functional Site Type\n(Fraction of oncogenic variants within each functional site type)", fontsize=14, pad=10)
plt.xlabel("Functional Site Type", fontsize=12)
plt.ylabel("Gene", fontsize=12)

plt.tight_layout()
plt.savefig("plots/functionalsites/top_genes_per_functional_site.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/functionalsites/top_genes_functional_site.png'")

# ------------------------------------------------------------
# Oncogenic-Neutral Variant Ratio in the top contributing genes
# ------------------------------------------------------------

print("-"*50)
print("Comparing oncogenic and likely neutral variants in top contributing genes..")

print("Extracting likely neutral variants..")

comparison = (
    expanded.groupby(["Hugo_Symbol", "FEATURE_TYPE", "ONCOGENIC"])
    .size()
    .unstack(fill_value=0)
)
comparison.columns = ["Likely Neutral", "Oncogenic"]
comparison = comparison.reset_index() 

# calculate ration with pseudocount (to handle zero values)
comparison["Ratio"] = (comparison["Oncogenic"] + 1) / (comparison["Likely Neutral"] + 1)

# preview data
print("Oncogenic and likely neutral variant counts per gene and functional site:")
print(comparison.head(),"\n") 

# ------------------------------------------------------------
# Plot oncogenic-neutral ratios
# ------------------------------------------------------------

print("-"*50)
print("Plotting oncogenic-neutral ratio for all functional site types..")

minimum_var = 20
for ft in filtered_sites: 
    subset = comparison[comparison["FEATURE_TYPE"] == ft].copy() 

    # only look at genes with a high number of oncogenic variants 
    subset = subset[subset["Oncogenic"] > 2]
    subset = subset.sort_values("Ratio", ascending=False).head(15) 

    if subset.empty: 
        continue 

    plt.figure(figsize=(8,5)) 
    sns.barplot(data=subset,
                x="Hugo_Symbol",
                y="Ratio",
                color="#c4314a",
                edgecolor="0.1",
                linewidth=0.3)
    
    plt.title(f"Oncogenic-to-Neutral Variant Ratio — {ft}")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Oncogenic / Neutral ratio (pseudocount +1)")
    plt.tight_layout()
    plt.savefig(f"plots/functionalsites/onco_neutral_ratio_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved in 'plots/functionalsites/'")
print("-"*50)

print("\nFunctional site analysis complete!🥳\n")