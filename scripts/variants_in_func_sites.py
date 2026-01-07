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
    explore_cancer_variants/plots/

====================================================================
"""

print("\n========================================================")
print("VARIANT FUNCTIONAL SITES ANALYSIS")
print("========================================================")

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ------------------------------------------------------------
# Setup for analysis 
# ------------------------------------------------------------

# Create output folder if it doesn't exist
os.makedirs("explore_cancer_variants/plots", exist_ok=True)

print("\n------------------------------------------------------")
print("LOAD VARIANT DATA")
print("------------------------------------------------------\n")

print("Loading variant data...\n")

variants = pd.read_csv(
    "annotation_pipeline/output/variants_with_func_sites.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.")

# ------------------------------------------------------------
# Keep only Oncogenic + Likely Neutral variants
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("FILTER VARIANTS")
print("------------------------------------------------------\n")

print("Filtering data to only contain oncogenic & likely oncogenic variants...\n")

classes = ["Oncogenic", "Likely Neutral"]
variants = variants[variants["ONCOGENIC"].isin(classes)]
print(f"Remaining variants after filtering: {len(variants):,}")

# Ensure boolean data type
variants["IN_FUNC_SITE"] = variants["IN_FUNC_SITE"].astype(bool)

# ------------------------------------------------------------
# Summary: variants inside vs outside functional sites
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("SUMMARY VARIANTS INSIDE VS OUTSIDE FUNCTIONAL SITES")
print("------------------------------------------------------\n")

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

print("\n------------------------------------------------------")
print("Expand FEATURE_TYPE")
print("------------------------------------------------------\n")

expanded = (
    variants
    .dropna(subset=["FEATURE_TYPE"])
    .assign(FEATURE_TYPE=lambda df: df["FEATURE_TYPE"].str.split(";"))
    .explode("FEATURE_TYPE")
)

print("Expanding FEATURE_TYPE so each type is one row...")
expanded["FEATURE_TYPE"] = expanded["FEATURE_TYPE"].str.strip()
print(f"\nExpanded to {len(expanded):,} feature-variant rows.\n")

# ------------------------------------------------------------
# Count variants per feature type per class
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("COUNT VARIANTS PER FEATURE TYPE")
print("------------------------------------------------------\n")

counts = (
    expanded
    .groupby(["FEATURE_TYPE", "ONCOGENIC"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Number of variants in each functional site for both classes:")
print(counts)

# ------------------------------------------------------------
# Plot raw counts
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("PLOT COUNTS OF VARIANTS IN FUNCTIONAL SITES")
print("------------------------------------------------------\n")

print("Plotting raw counts of variants in functional sites...\n")

palette = {"Oncogenic": "#C4473B", "Likely Neutral": "#7e8aa2"}

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
plt.title("Number of Variants per Functional Site Type", fontsize=14, pad=10)
plt.xlabel("Functional Site", fontsize=12)
plt.ylabel("Variant Count", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("explore_cancer_variants/plots/counts_per_feature_type.png", dpi=300)
plt.show()

print(f"Plotting complete. Saved as 'explore_cancer_variants/plots/counts_per_feature_type.png'")

# ------------------------------------------------------------
# Compute fractions
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("COMPUTE FRACTIONS")
print("------------------------------------------------------\n")

print("Computing fractions of variants in the different functional sites...")

totals = variants["ONCOGENIC"].value_counts().rename("Total")

counts = counts.merge(totals, left_on="ONCOGENIC", right_index=True)
counts["Fraction"] = counts["Variant_Count"] / counts["Total"]

print("\nFraction of variants in each feature type per class:\n")
print(counts)

# ------------------------------------------------------------
# Plot fractions
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("PLOT FRACTIONS OF VARIANTS IN FUNCTIONAL SITES")
print("------------------------------------------------------\n")

print("Plotting fractions of variants for each feature type...\n")

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
plt.title("Fraction of Variants per Feature Type", fontsize=14, pad=10)
plt.xlabel("Functional Site Type", fontsize=12)
plt.ylabel("Fraction", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("explore_cancer_variants/plots/fraction_per_feature_type.png", dpi=300)
plt.show()

print("Plotting complete. Saved as 'explore_cancer_variants/plots/fraction_per_feature_type.png'")

# ------------------------------------------------------------
# Identify genes enriched in functional sites
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("FIND GENES ENRICHED IN FUNCTIONAL SITES")
print("------------------------------------------------------\n")

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
print(gene_feature_fraction.head(5))

# ------------------------------------------------------------
# Plot top genes per functional site
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("PLOT TOP GENES PER FUNCTIONAL SITE")
print("------------------------------------------------------\n")

print("Plotting top genes per functional site...\n")

for ft in filtered_sites:
    subset = (
        gene_feature_fraction[gene_feature_fraction["FEATURE_TYPE"] == ft]
        .sort_values("Fraction_of_Feature", ascending=False)
        .head(20)
    )
    
    plt.figure(figsize=(8,5))

    plt.bar(subset["Hugo_Symbol"], 
            subset["Fraction_of_Feature"], 
            color="#C4473B", 
            edgecolor="0.1",
            linewidth=0.3)

    plt.title(f"Top Driver Genes in '{ft}'", fontsize=14, pad=10)
    plt.xlabel("Gene", fontsize=12)
    plt.ylabel("Fraction of variants in feature", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    plt.tight_layout()
    plt.savefig(f"explore_cancer_variants/plots/topgenes_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved in 'explore_cancer_variants/plots'\n")

# ------------------------------------------------------------
# Compare likely neutral variants to oncogenic variants in top oncogenic driver genes
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("COMPARE VARIANTS IN TOP ONCOGENIC DRIVER GENES")
print("------------------------------------------------------\n")


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

print("Example output oncogenic-neutral ratio:")
print(comparison_top.head())

# ------------------------------------------------------------
# Plot oncogenic-neutral ratios
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("PLOT ONCOGENIC-NEUTRAL RATIOS")
print("------------------------------------------------------\n")

print("Plotting oncogenic-neutral ratio for all feature types...\n")

for ft in filtered_sites:
    subset = (
        comparison_top[comparison_top["FEATURE_TYPE"] == ft]
        .sort_values("ratio_onco_neutral", ascending=False)
        .head(20)
    )

    plt.figure(figsize=(8,5))
    plt.bar(subset["Hugo_Symbol"], 
            subset["ratio_onco_neutral"], 
            color="#C4473B",
            edgecolor="0.1",
            linewidth=0.3)

    plt.title(f"Oncogenic vs Neutral Variant Ratio in '{ft}'", fontsize=14, pad=10)
    plt.xlabel("Gene", fontsize=12)
    plt.ylabel("Ratio", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    plt.tight_layout()
    plt.savefig(f"explore_cancer_variants/plots/onco-neutral-ratio_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved in 'explore_cancer_variants/plots'\n")

# ------------------------------------------------------------
# Which genes dominate in each feature type? 
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("FIND DOMINATING GENES FOR EACH FEATURE TYPE")
print("------------------------------------------------------\n")

# Pick top genes from each feature type

print("Extracting top genes per feature type...\n")
top_genes_per_feature = (
    gene_feature_fraction
    .groupby('FEATURE_TYPE')
    .apply(lambda x: x.nlargest(5, 'Fraction_of_Feature'))
    .reset_index(drop=True)
)

print("Example output:")
print(top_genes_per_feature.head(10),"\n")

print("Plotting top genes per feature...\n")

pivot = top_genes_per_feature.pivot(
    index='Hugo_Symbol', 
    columns='FEATURE_TYPE', 
    values='Fraction_of_Feature'
).fillna(0)

plt.figure(figsize=(10,6))
ax = sns.heatmap(pivot, 
            cmap='Reds', 
            annot=True, 
            fmt='.2f',
            cbar_kws={"label":"Fraction of variants in feature"},
            )

ax.set_xlabel("Feature Type", fontsize=12)
ax.set_ylabel("Gene (Hugo Symbol)", fontsize=12)

plt.title('Gene Distribution Across Functional Site Types', fontsize=14, pad=10)

plt.tight_layout()
plt.savefig("explore_cancer_variants/plots/top_genes_per_functional_site.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'explore_cancer_variants/plots/top_genes_per_functional_site.png'\n")



print("========================================================")
print("\nVariant Functional Site Analysis complete!🥳\n")