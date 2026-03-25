
# ====================================================================
# Protein Domain Analysis 
# ====================================================================

"""
Script: protein_domains.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across different protein domains and genes.  

Major outputs:
  1. Variant distribution inside/outside protein domains 
  2. Top domains (by variant count)
  3. Top domains enriched for oncogenic variants 
  4. Heatmap of variants across top domains and top genes 
  5. Driver genes enriched in protein domains 
  6. Heatmap of oncogenic variants across top domains and top genes

All plots are saved in:
    plots/proteindomains

"""

print("-"*50)
print("Protein Domain Analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Explore protein domains in variant data."
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
# Load variant data
# ------------------------------------------------------------

print("Loading variant data..")

variants = pd.read_csv(
    args.variants,
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.")

# ------------------------------------------------------------
# Count variants inside and outside protein domains
# ------------------------------------------------------------

print("-"*50)
print("Counting variants inside/outside protein domains..\n")

oncogenic = variants[variants['ONCOGENIC'] == 'Oncogenic']
neutral = variants[variants['ONCOGENIC'] == 'Likely Neutral']

def count_in_out(df, class_label): 
  inside = df["DOMAIN_NAME"].notna().sum()
  outside = df["DOMAIN_NAME"].isna().sum() 
  total = len(df) 

  print(f"{class_label} variants inside protein domain: {inside:,}")
  print(f"{class_label} variants outside protein domain: {outside:,}")
  print(f"Total {class_label.lower()} variants: {total:,}")
  print(f"Fraction inside domain: {inside/total:.2%}\n")
  return inside, outside, total 

print("ONCOGENIC")
print(count_in_out(oncogenic, "Oncogenic"),"\n")

print("LIKELY NEUTRAL")
print(count_in_out(neutral, "Likely Neutral")) 

# ------------------------------------------------------------
# Expand DOMAIN_NAME for variants with multiple domains 
# ------------------------------------------------------------

print("-"*50)
print("Exploding variants with multiple domains..")

variants_domains = (
  variants
  .dropna(subset=["DOMAIN_NAME"])
  .assign(DOMAIN_NAME = lambda df: df["DOMAIN_NAME"].str.split(";"))
  .explode("DOMAIN_NAME")
  .reset_index(drop=True) 
  )

variants_domains["DOMAIN_NAME"] = variants_domains["DOMAIN_NAME"].str.strip() 

print(f"After exploding: {len(variants_domains):,} domain-variant rows.")

# ------------------------------------------------------------
# Find top protein domains 
# ------------------------------------------------------------

# Extract oncogenic and likely neutral variants
variants_filtered = variants_domains[variants_domains["ONCOGENIC"].isin(["Oncogenic", "Likely Neutral"])]

# Number of variants per domain and class 
domain_counts = pd.crosstab(variants_filtered["DOMAIN_NAME"], variants_filtered["ONCOGENIC"])

# Top 15 domain names by variant count 
top_15_names = domain_counts.sum(axis=1).sort_values(ascending=False).head(15).index

print(f"The top protein domains in this dataset are:")
print(list(top_15_names))

# ------------------------------------------------------------
# Plot Top Domains by Number of Total Variants 
# ------------------------------------------------------------

# Locate the variants within top domains 
domain_plot_data = domain_counts.loc[top_15_names]

# Plot variant counts in top domains 
print("-"*50)
print("Plotting top protein domains by variant counts..\n")

ax = domain_plot_data.plot(
    kind="bar",
    stacked=True,
    figsize=(8,5),
    edgecolor="0.1",
    linewidth=0.3,
    color={"Oncogenic": "#c4314a", "Likely Neutral": "#88aed1"}
)

plt.title("Top Protein Domains by Variant Count", fontsize=14, pad=15)
plt.xlabel("Domain Name", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12)
plt.xticks(rotation=90, ha="right", fontsize=9)
plt.yticks(fontsize=9) 
plt.legend(title="Oncogenicity Class")

plt.tight_layout()
plt.savefig("plots/proteindomains/top_domains.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'plots/proteindomains/top_domains.png'")

# ------------------------------------------------------------
# Oncogenic vs. Neutral Enrichment per Domain 
# ------------------------------------------------------------

print("-"*50)
print("Computing oncogenic vs. neutral enrichment per domain..")

# Count number of variants per class and domain 
domain_class_counts = (
  variants_filtered
  .groupby(["DOMAIN_NAME","ONCOGENIC"])
  .size() 
  .reset_index(name="Count") 
)

# Format into pivot 
domain_pivot = domain_class_counts.pivot(
    index="DOMAIN_NAME",
    columns="ONCOGENIC",
    values="Count"
  ).fillna(0) 

# Calculate ratio 
domain_pivot["Onco_Neutral_Ratio"] = (
  (domain_pivot["Oncogenic"] + 1) / 
  (domain_pivot["Likely Neutral"] + 1) 
)

domain_enrichment = domain_pivot.sort_values("Onco_Neutral_Ratio", ascending=False)

print("\nPreview of the domain enrichment table:\n")
print(domain_enrichment.head(10))


# ------------------------------------------------------------
# Plot Top Domains Enriched for Oncogenic Variants 
# ------------------------------------------------------------

print("-"*50)
print("Plotting top protein domains enriched for oncogenic variants..\n")

plt.figure(figsize=(8,5))
sns.barplot(
  data=domain_enrichment.head(15),
  x=domain_enrichment.head(15).index, 
  y="Onco_Neutral_Ratio",
  palette="Reds_r",
  edgecolor="0.1",
  linewidth=0.3
)

plt.title("Domains Enriched for Oncogenic Variants", fontsize=14, pad=10)
plt.xlabel("Domain Name", fontsize=12)
plt.ylabel("Ratio (Oncogenic-Neutral)", fontsize=12)
plt.xticks(rotation=90, ha="right", fontsize=9)
plt.yticks(fontsize=9) 

plt.tight_layout()
plt.savefig("plots/proteindomains/domain_oncogenic_enrichment.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'plots/proteindomains/domain_oncogenic_enrichment.png'")

# ------------------------------------------------------------
# Combined Oncogenic/Neutral Heatmap 
# ------------------------------------------------------------

# Find the oncogenic fraction of variants in top genes and top domains 

print("-"*50)
print("Plotting combined heatmap (top domains x top genes)..")

# Count oncogenic + neutral variants per domain and gene 
gene_domain_class_counts = (
  variants_filtered
  .groupby(["Hugo_Symbol", "DOMAIN_NAME", "ONCOGENIC"])
  .size() 
  .reset_index(name="Count") 
)
print("Preview of the gene x domain count table:\n")
print(gene_domain_class_counts.head(), "\n")

gene_domain_matrix = gene_domain_class_counts.pivot_table(
    index=["Hugo_Symbol", "DOMAIN_NAME"],
    columns="ONCOGENIC",
    values="Count",
    fill_value=0
).reset_index()

# Compute oncogenic fraction 
gene_domain_matrix["Total"] = (
  gene_domain_matrix["Oncogenic"] + gene_domain_matrix["Likely Neutral"] 
)

gene_domain_matrix["Oncogenic_Fraction"] = (
  gene_domain_matrix["Oncogenic"] / gene_domain_matrix["Total"]
).fillna(0) 

# Filter to top domains 
top_domains = (
    gene_domain_matrix.groupby("DOMAIN_NAME")["Total"]
    .sum() 
    .sort_values(ascending=False) 
    .head(15)
    .index
)

# Filter to top domains
combined_top = gene_domain_matrix[
    gene_domain_matrix["DOMAIN_NAME"].isin(top_domains)
]

# Filter to top genes 
top_genes_combined = (
    combined_top.groupby("Hugo_Symbol")["Total"]
    .sum()
    .sort_values(ascending=False)
    .head(15)
    .index
)

combined_top = combined_top[combined_top["Hugo_Symbol"].isin(top_genes_combined)]

# Pivot for heatmap
heatmap_combined = combined_top.pivot(
    index="Hugo_Symbol",
    columns="DOMAIN_NAME",
    values="Oncogenic_Fraction"
).fillna(0)

# Plot heatmap 
plt.figure(figsize=(7,6))
sns.heatmap(heatmap_combined, 
            cmap="Reds", 
            vmin=0, 
            vmax=1, 
            linewidths=0.2)

plt.title("Oncogenic Fraction per Gene × Domain", fontsize=14, pad=12)
plt.xlabel("Domain Name", fontsize=12)
plt.ylabel("Gene (Hugo Symbol)", fontsize=12)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)

plt.tight_layout()
plt.savefig("plots/proteindomains/heatmap_oncogenic_fraction.png", dpi=300, bbox_inches="tight")
plt.show()

print("Heatmap complete! Saved as 'plots/proteindomains/heatmap_oncogenic_fraction.png'")

# ------------------------------------------------------------
# Identify Driver Genes Enriched in Protein Domains
# ------------------------------------------------------------

# When we look at a specific protein domain, how many of the oncogenic mutations comes from each gene? 

print("-"*50)
print("Identifying oncogenic driver genes enriched in protein domains..\n")

# Extract oncogenic variants from the exploded df 
oncogenic_variants = variants_domains[variants_domains["ONCOGENIC"] == "Oncogenic"].copy()

# Count oncogenic variants per gene x domain 
# "How many oncogenic variants does each gene have in each domain?"
gene_domain_counts = (
    oncogenic_variants.groupby(["Hugo_Symbol", "DOMAIN_NAME"])
    .size()
    .reset_index(name="Variant_Count")
    .sort_values("Variant_Count", ascending=False)
)

# Compute total oncogenic variants per domain 
domain_oncogenic_totals = (
    gene_domain_counts.groupby("DOMAIN_NAME")["Variant_Count"]
    .sum()
    .reset_index(name="Domain_Total")
)

# Compute the fraction contributed per gene 
gene_domain_fraction = gene_domain_counts.merge(
  domain_oncogenic_totals, on="DOMAIN_NAME")

gene_domain_fraction["Fraction_of_Domain"] = (
    gene_domain_fraction["Variant_Count"] / gene_domain_fraction["Domain_Total"]
)

print("Preview of oncogenic driver genes:\n")
print(gene_domain_fraction.head(5), "\n")

# ------------------------------------------------------------
# Top Domains x Top Genes 
# (Oncogenic variants only)
# ------------------------------------------------------------

print("-"*50)
print("Plotting heatmap of top domains x top genes (oncogenic variants)...\n")

n_genes = 15
n_domains = 15

# Top domains by total oncogenic variants 
top_domains = (
  gene_domain_fraction.groupby("DOMAIN_NAME")["Variant_Count"]
  .sum()
  .sort_values(ascending=False) 
  .head(n_domains)
  .index
)

# Filter to top domains
gene_domain_fraction_top = gene_domain_fraction[gene_domain_fraction["DOMAIN_NAME"].isin(top_domains)]

# Pick top genes across top domains 
top_genes = (
  gene_domain_fraction_top.groupby("Hugo_Symbol")["Variant_Count"]
  .sum() 
  .sort_values(ascending=False) 
  .head(n_genes) 
  .index
)

# Filter to top genes 
gene_domain_fraction_top = gene_domain_fraction_top[gene_domain_fraction_top["Hugo_Symbol"].isin(top_genes)]

# Pivot for heatmap 
heatmap_df = gene_domain_fraction_top.pivot(
  index="Hugo_Symbol",
  columns="DOMAIN_NAME",
  values="Fraction_of_Domain"
).fillna(0) 


# Plot heatmap
plt.figure(figsize=(7,6))
sns.heatmap(heatmap_df, 
            cmap="Reds", 
            vmin=0, 
            vmax=1,
            linewidths=0.2)

plt.title("Top genes x Top domains (Oncogenic Variants)", fontsize=14, pad=12)
plt.xlabel("Domain Name", fontsize=12)
plt.ylabel("Gene (Hugo Symbol)", fontsize=12)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)

plt.tight_layout()
plt.savefig("plots/proteindomains/heatmap_oncogenic_top.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/proteindomains/heatmap_oncogenic_top.png'\n")
print("-"*50)


print("\nProtein domain analysis complete!🥳🥳\n")