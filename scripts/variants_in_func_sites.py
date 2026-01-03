"""
====================================================================
Variant Functional-Site Analysis Script
====================================================================

Script: variants_in_func_sites.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across different functional protein features.

Major outputs:
--------------
1. Counts of Oncogenic vs Likely Neutral variants inside/outside functional sites
2. Expansion of multi-feature annotations (e.g., "Binding site; Region")
3. Variant counts per functional site type per class
4. Fractions of variants in each feature type
5. Identification of genes enriched in functional sites
6. Comparison of Oncogenic vs Neutral variants in the same genes
7. Oncogenic-to-neutral ratio plots for each feature type
8. Heatmap to see which genes dominate in each feature type

Key features handled:
---------------------
- Binding site
- Modified residue
- Region
- Topological domain

All plots are saved in:
    visualize_variants/plots/

====================================================================
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ------------------------------------------------------------
# Setup
# ------------------------------------------------------------

# Create output folder if it doesn't exist
os.makedirs("visualize_variants/plots", exist_ok=True)

print("\nLoading variant data...\n")

variants = pd.read_csv(
    "annotation_pipeline/output/variants_with_func_sites.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.\n")

# ------------------------------------------------------------
# Keep only Oncogenic + Likely Neutral variants
# ------------------------------------------------------------

classes = ["Oncogenic", "Likely Neutral"]
variants = variants[variants["ONCOGENIC"].isin(classes)]
print(f"Remaining variants after filtering: {len(variants):,}\n")

# Ensure boolean data type
variants["IN_FUNC_SITE"] = variants["IN_FUNC_SITE"].astype(bool)

# ------------------------------------------------------------
# Summary: variants inside vs outside functional sites
# ------------------------------------------------------------

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
print("Distribution of variants inside and outside functional sites:\n")
print(summary_df)

# ------------------------------------------------------------
# Expand FEATURE_TYPE 
# ------------------------------------------------------------

expanded = (
    variants
    .dropna(subset=["FEATURE_TYPE"])
    .assign(FEATURE_TYPE=lambda df: df["FEATURE_TYPE"].str.split(";"))
    .explode("FEATURE_TYPE")
)

print(" \n\nExpanding FEATURE_TYPE so each type is one row...")
expanded["FEATURE_TYPE"] = expanded["FEATURE_TYPE"].str.strip()
print(f"\nExpanded to {len(expanded):,} feature-variant rows.\n")

# ------------------------------------------------------------
# Count variants per feature type per class
# ------------------------------------------------------------

counts = (
    expanded
    .groupby(["FEATURE_TYPE", "ONCOGENIC"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Number of variants in each functional site for both classes:\n")
print(counts, "\n")

# ------------------------------------------------------------
# Plot raw counts
# ------------------------------------------------------------

print("Plotting raw counts of variants in functional sites...")

palette = {"Oncogenic": "#ef6f6c", "Likely Neutral": "#5b8fdc"}

plt.figure(figsize=(6,4))
sns.barplot(
    data=counts,
    x="FEATURE_TYPE",
    y="Variant_Count",
    hue="ONCOGENIC",
    palette=palette
)
plt.title("Number of Variants per Functional Site Type", fontweight="bold", fontsize=13)
plt.xlabel("Functional Site")
plt.ylabel("Variant Count")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("visualize_variants/plots/counts_per_feature_type.png", dpi=300)
plt.show()

print(f"Plotting complete. Saved as 'visualize_variants/plots/counts_per_feature_type.png'\n")

# ------------------------------------------------------------
# Compute fractions
# ------------------------------------------------------------

print("Computing fractions of variants in the different functional sites...")

totals = variants["ONCOGENIC"].value_counts().rename("Total")

counts = counts.merge(totals, left_on="ONCOGENIC", right_index=True)
counts["Fraction"] = counts["Variant_Count"] / counts["Total"]

print("\nFraction of variants in each feature type per class:\n")
print(counts, "\n")

# ------------------------------------------------------------
# Plot fractions
# ------------------------------------------------------------

print("Plotting fractions of variants for each feature type...")

plt.figure(figsize=(6,4))
sns.barplot(
    data=counts,
    x="FEATURE_TYPE",
    y="Fraction",
    hue="ONCOGENIC",
    palette=palette
)
plt.title("Fraction of Variants per Feature Type", fontweight="bold", fontsize=13)
plt.xlabel("Functional Site Type")
plt.ylabel("Fraction")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("visualize_variants/plots/fraction_per_feature_type.png", dpi=300)
plt.show()

print("Plotting complete. Saved as 'visualize_variants/plots/fraction_per_feature_type.png'\n")

# ------------------------------------------------------------
# Identify genes enriched in functional sites
# ------------------------------------------------------------

print("Identifying oncogenic driver genes enriched in functional sites..\n")

filtered_sites = ["Binding site", "Modified residue", "Region", "Topological domain"]

expanded_filtered = expanded[expanded["FEATURE_TYPE"].isin(filtered_sites)].copy()

# Oncogenic counts
onco = expanded_filtered[expanded_filtered["ONCOGENIC"] == "Oncogenic"].copy()

gene_feature_counts = (
    onco.groupby(["Hugo_Symbol", "FEATURE_TYPE"])
    .size()
    .reset_index(name="Variant_Count")
    .sort_values("Variant_Count", ascending=False)
)

feature_totals = (
    gene_feature_counts.groupby("FEATURE_TYPE")["Variant_Count"]
    .sum()
    .reset_index(name="Feature_Total")
)

gene_feature_fraction = gene_feature_counts.merge(feature_totals, on="FEATURE_TYPE")
gene_feature_fraction["Fraction_of_Feature"] = (
    gene_feature_fraction["Variant_Count"] / gene_feature_fraction["Feature_Total"]
)

print("Example output:\n")
print(gene_feature_fraction.head(5), "\n")

# ------------------------------------------------------------
# Plot top genes per functional site
# ------------------------------------------------------------

print("Plotting top genes per functional site...\n")

for ft in filtered_sites:
    subset = (
        gene_feature_fraction[gene_feature_fraction["FEATURE_TYPE"] == ft]
        .sort_values("Fraction_of_Feature", ascending=False)
        .head(10)
    )
    
    plt.figure(figsize=(6,4))
    plt.bar(subset["Hugo_Symbol"], subset["Fraction_of_Feature"], color="#ef6f6c")
    plt.title(f"Top Driver Genes in '{ft}'", fontweight="bold", fontsize=13)
    plt.xlabel("Gene")
    plt.ylabel("Fraction of variants in feature")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"visualize_variants/plots/topgenes_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved to output folder.\n") 

# ------------------------------------------------------------
# Compare likely neutral variants in top oncogenic driver genes
# ------------------------------------------------------------

print("Extracting likely neutral variants...\n")

likely_neutral = expanded_filtered[expanded_filtered["ONCOGENIC"] == "Likely Neutral"].copy()

print("Counting likely neutral variants per feature type...\n")
likely_neutral_counts = (
    likely_neutral.groupby(["Hugo_Symbol", "FEATURE_TYPE"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Merging oncogenic variants to likely neutral variants...\n")
comparison = gene_feature_counts.merge(
    likely_neutral_counts,
    on=["Hugo_Symbol", "FEATURE_TYPE"],
    how="outer",
    suffixes=("_oncogenic", "_likely_neutral")
).fillna(0)

top_onco_genes = gene_feature_counts["Hugo_Symbol"].unique()
comparison_top = comparison[comparison["Hugo_Symbol"].isin(top_onco_genes)].copy()

print("Calculating the oncogenic-neutral ratio for the top driver genes...\n")
comparison_top["ratio_onco_neutral"] = (
    (comparison_top["Variant_Count_oncogenic"] + 1) /
    (comparison_top["Variant_Count_likely_neutral"] + 1)
)

print("Example output oncogenic-neutral ratio:\n")
print(comparison_top.head(), "\n")

# ------------------------------------------------------------
# Plot oncogenic-neutral ratios
# ------------------------------------------------------------

print("Plotting oncogenic-neutral ratio for all feature types...\n")

for ft in filtered_sites:
    subset = (
        comparison_top[comparison_top["FEATURE_TYPE"] == ft]
        .sort_values("ratio_onco_neutral", ascending=False)
        .head(20)
    )

    plt.figure(figsize=(7,4))
    plt.bar(subset["Hugo_Symbol"], subset["ratio_onco_neutral"], color="#ef6f6c")
    plt.title(f"Oncogenic vs Neutral Variant Ratio in '{ft}'", fontweight="bold", fontsize=13)
    plt.xlabel("Gene")
    plt.ylabel("Ratio")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"visualize_variants/plots/onco-neutral-ratio_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved to output folder.\n")

# ------------------------------------------------------------
# Which genes dominate in each feature type? 
# ------------------------------------------------------------

# Pick top 10 genes from each feature type
top_genes_per_feature = (
    gene_feature_fraction
    .groupby('FEATURE_TYPE')
    .apply(lambda x: x.nlargest(5, 'Fraction_of_Feature'))
    .reset_index(drop=True)
)

pivot = top_genes_per_feature.pivot(
    index='Hugo_Symbol', 
    columns='FEATURE_TYPE', 
    values='Fraction_of_Feature'
).fillna(0)

plt.figure(figsize=(10,8))
sns.heatmap(pivot, cmap='Reds', annot=True, fmt='.2f')
plt.title('Gene Distribution Across Functional Site Types')
plt.tight_layout()
plt.show()

print("\nVariant Functional Site Analysis complete!🥳\n")