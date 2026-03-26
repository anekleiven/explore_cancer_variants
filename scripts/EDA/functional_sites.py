"""
Functional Sites Analysis 
-----------------------------------------------

Script: functional_sites.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across different functional protein features.

Major outputs:
--------------
1. Counts of Oncogenic vs Likely Neutral variants inside/outside functional sites
2. Variant counts per functional site type per class
3. Fractions of variants per functional site type per class 
4. Top genes by oncogenic variant contribution per functional site
5. Heatmap of top genes in filtered functional sites  
6. Comparison of Oncogenic vs Neutral variants in the same genes
7. Oncogenic-to-neutral ratio plots for each functional site type

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
print("Computing gene variant contributions per functional site..")

# Oncogenic counts
oncogenic = expanded[expanded["ONCOGENIC"] == "Oncogenic"].copy()

# Number of oncogenic variants per gene and feature type
gene_feature_counts = (
    oncogenic.groupby(["Hugo_Symbol", "FEATURE_TYPE"])
    .size()
    .reset_index(name="Variant_Count")
    .sort_values("Variant_Count", ascending=False)
)

# Total number of oncogenic variants per feature type 
feature_totals = (
    gene_feature_counts.groupby("FEATURE_TYPE")["Variant_Count"]
    .sum()
    .rename("Feature_Total")
)

# Fraction of oncogenic variants contributed by each gene 
gene_feature_fraction = gene_feature_counts.merge(feature_totals, on="FEATURE_TYPE")
gene_feature_fraction["Fraction_of_Feature"] = (
    gene_feature_fraction["Variant_Count"] / gene_feature_fraction["Feature_Total"]
)

print("Example output gene-feature fraction table:")
print(gene_feature_fraction.head(5))

# ------------------------------------------------------------
# Plot oncogenic variant contribution by gene
# per functional site 
# ------------------------------------------------------------

print("-"*50)
print("Plotting gene variant contributions per functional site..")

# Minimum number of variants for a site to be plotted 
minimum_var = 20 

top_genes_per_feature = set()  

# Loop over all feature types
for ft in gene_feature_fraction["FEATURE_TYPE"].unique():

    total_variants = feature_totals.loc[ft]

    # Skip features with too few variants
    if total_variants < minimum_var:
        print(f"Skipping {ft}. Only {total_variants} variants (< {minimum_var}).")
        continue 

    subset = (
        gene_feature_fraction[gene_feature_fraction["FEATURE_TYPE"] == ft]
        .sort_values("Fraction_of_Feature", ascending=False)
        .head(20)
    )

    top_genes_per_feature.update(subset["Hugo_Symbol"].tolist())
    
    fig, ax = plt.subplots(figsize=(8,5)) 

    ax.bar(subset["Hugo_Symbol"],
           subset["Fraction_of_Feature"], 
           color="#c4314a", 
           edgecolor="0.1",
           linewidth=0.3)
    
    ax.set_title(f"Oncogenic variant contribution by gene: {ft}", fontsize=14, pad=10)
    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Fraction of oncogenic feature variants", fontsize=12)
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
print("Extracting top contributing genes per functional site type..")

# filter to functional sites with the most variants 
filtered_sites = ["Binding site", "Modified residue", "Region", "Topological domain"]

# find the top 5 genes with the highest amount of oncogenic variants per site (fraction) 
top_genes_per_feature = (
    gene_feature_fraction[gene_feature_fraction["FEATURE_TYPE"].isin(filtered_sites)]
    .sort_values(['FEATURE_TYPE', 'Fraction_of_Feature'], ascending=[True,False])
    .groupby('FEATURE_TYPE')
    .head(5)
    .reset_index(drop=True)
)

print("Example output:")
print(top_genes_per_feature.head(10),"\n")

print("Plotting top genes per functional site..")

# create matrix for heatmap 
pivot = top_genes_per_feature.pivot(
    index='Hugo_Symbol', 
    columns='FEATURE_TYPE', 
    values='Fraction_of_Feature'
).fillna(0)

plt.figure(figsize=(7,6))
ax = sns.heatmap(pivot, 
            cmap='Reds', 
            annot=True, 
            fmt='.2f',
            cbar_kws={"label":"Fraction of oncogenic variants in functional site"},
            )

ax.set_title("Gene Distribution Across Functional Site Types", fontsize=14, pad=10)
ax.set_xlabel("Functional Site", fontsize=12)
ax.set_ylabel("Gene (Hugo Symbol)", fontsize=12)

plt.tight_layout()
plt.savefig("plots/functionalsites/top_genes_per_functional_site.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/functionalsites/top_genes_functional_site.png'")

# ------------------------------------------------------------
# Compare likely neutral variants to oncogenic variants 
# in top contributing genes (by oncogenic variant count) 
# ------------------------------------------------------------

print("-"*50)
print("Comparing oncogenic and likely neutral variants in top contributing genes..")

print("Extracting likely neutral variants..")

likely_neutral = expanded[expanded["ONCOGENIC"] == "Likely Neutral"].copy()

print("Counting likely neutral variants per functional site type..")
likely_neutral_counts = (
    likely_neutral.groupby(["Hugo_Symbol", "FEATURE_TYPE"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Merging oncogenic and likely neutral count tables..\n")
comparison = gene_feature_counts.merge(
    likely_neutral_counts,
    on=["Hugo_Symbol", "FEATURE_TYPE"],
    how="outer",
    suffixes=("_oncogenic", "_likely_neutral")
).fillna(0)

print("Oncogenic and likely neutral variant counts per gene and functional site:")
print(comparison.head(),"\n") 

print("Filtering to top oncogenic variant contributors..")
comparison_top = comparison[comparison["Hugo_Symbol"].isin(top_genes_per_feature)]

print("Calculating the oncogenic-neutral ratio for the top contributing genes..")
comparison_top["ratio_onco_neutral"] = (
    (comparison_top["Variant_Count_oncogenic"] + 1) /
    (comparison_top["Variant_Count_likely_neutral"] + 1)
)

print("Example output oncogenic-neutral ratio:")
print(comparison_top.head(),"\n")

# ------------------------------------------------------------
# Plot oncogenic-neutral ratios
# ------------------------------------------------------------

print("-"*50)
print("Plotting oncogenic-neutral ratio for all functional site types..")

minimum_var = 20

for ft in comparison_top["FEATURE_TYPE"].unique(): 

    # Functional site totals 
    n_onco_total = comparison_top[comparison_top["FEATURE_TYPE"] == ft]["Variant_Count_oncogenic"].sum()
    n_neutral_total = comparison_top[comparison_top["FEATURE_TYPE"] == ft]["Variant_Count_likely_neutral"].sum()
    total_variants = n_onco_total + n_neutral_total

    # Skip features with too few variants
    if total_variants < minimum_var:
        print(f"Skipping {ft}. Only {total_variants} variants (< {minimum_var}).")
        continue

    subset = (
        comparison_top[comparison_top["FEATURE_TYPE"] == ft]
        .sort_values("ratio_onco_neutral", ascending=False)
    )

    # Only keep genes above the pseudocount baseline
    subset = subset[subset["ratio_onco_neutral"] > 2]

    # Filtered totals for the plot annotation
    n_onco = subset["Variant_Count_oncogenic"].sum()
    n_neutral = subset["Variant_Count_likely_neutral"].sum()

    fig, ax = plt.subplots(figsize=(8,5)) 
    ax.bar(subset["Hugo_Symbol"], 
            subset["ratio_onco_neutral"], 
            color="#c4314a",
            edgecolor="0.1",
            linewidth=0.3)

    ax.set_title(f"Oncogenic vs. Neutral Variant Ratio in '{ft}'", fontsize=14, pad=10)
    ax.set_xlabel("Gene", fontsize=12)
    ax.set_ylabel("Oncogenic-Neutral Variant Ratio", fontsize=12)
    ax.tick_params(axis="x", rotation=45, labelsize=9)
    ax.tick_params(axis="y", labelsize=9)

    # Add number of variants in plot 
    ax.annotate(f"n = {int(n_onco)} oncogenic / {int(n_neutral)} neutral",
            xy=(0.99, 0.97), xycoords="axes fraction",
            ha="right", va="top", fontsize=9,
            color="0.4")

    plt.tight_layout()
    plt.savefig(f"plots/functionalsites/onco_neutral_ratio_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved in 'plots/functionalsites/'")
print("-"*50)

print("\nFunctional site analysis complete!🥳\n")