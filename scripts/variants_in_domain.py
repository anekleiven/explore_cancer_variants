"""
====================================================================
Variant Protein Domain Analysis Script
====================================================================

Script: variants_in_domain.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across different protein domains and genes 

Major outputs: 
1. See how oncogenic and likely neutral variants are distributed (inside/outside protein domains)
2. Explore oncogenic vs neutral enrichment per domain 
3. Plot top domains by number of total variants 
4. Plot top domains enriched for oncogenic variants 
5. Combined heatmap of oncogenic and neutral variants across top domains and top genes 
6. Identify driver genes enriched in protein domains
7. Heatmap of oncogenic variants across top domains and top genes 


All plots are saved in:
    explore_cancer_variants/plots/

"""

print("\n========================================================")
print("VARIANT PROTEIN DOMAIN ANALYSIS")
print("========================================================")

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("LOAD VARIANT DATA")
print("------------------------------------------------------\n")

print("Loading variant data...\n")

variants = pd.read_csv(
    "annotation_pipeline/output/variants_with_func_sites.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.\n")

# ------------------------------------------------------------
# Count variants inside and outside protein domains
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("VARIANT COUNTS INSIDE/OUTSIDE DOMAINS")
print("------------------------------------------------------\n")

print("Counting variants inside/outside protein domains...\n")

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

print("ONCOGENIC\n")
print(count_in_out(oncogenic, "Oncogenic"))

print("\n","-"*60)

print("\nLIKELY NEUTRAL")
print(count_in_out(neutral, "Likely Neutral")) 


# ------------------------------------------------------------
# Expand DOMAIN_NAME for variants with multiple domains 
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("EXPLODE MULTI-DOMAIN VARIANTS")
print("------------------------------------------------------\n")

print("Exploding variants with multiple domains...\n")

variants_domains = (
  variants
  .dropna(subset=["DOMAIN_NAME"])
  .assign(DOMAIN_NAME = lambda df: df["DOMAIN_NAME"].str.split(";"))
  .explode("DOMAIN_NAME")
  )

variants_domains["DOMAIN_NAME"] = variants_domains["DOMAIN_NAME"].str.strip() 

print(f"After exploding: {len(variants_domains):,} domain-variant rows.\n")

# ------------------------------------------------------------
# Oncogenic vs neutral enrichment per domain 
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("ONCOGENIC VS NEUTRAL ENRICHMENT PER DOMAIN")
print("------------------------------------------------------\n")

print("Computing oncogenic vs neutral enrichment per domain...")

# extract oncogenic and likely neutral variants
variants_domains = variants_domains[variants_domains["ONCOGENIC"].isin(["Oncogenic", "Likely Neutral"])]

# count number of variants per class and domain 
domain_class_counts = (
  variants_domains
  .groupby(["DOMAIN_NAME","ONCOGENIC"])
  .size() 
  .reset_index(name="Count") 
)

# format into pivot 
domain_pivot = domain_class_counts.pivot(
    index="DOMAIN_NAME",
    columns="ONCOGENIC",
    values="Count"
  ).fillna(0) 

# calculate ratio 
domain_pivot["Onco_Neutral_Ratio"] = (
  (domain_pivot["Oncogenic"] + 1) / 
  (domain_pivot["Likely Neutral"] + 1) 
)

domain_enrichment = domain_pivot.sort_values("Onco_Neutral_Ratio", ascending=False)

print("\nPreview of the domain enrichment table:\n")
print(domain_enrichment.head(10), "\n")

# ------------------------------------------------------------
# Plot top domains by number of total variants 
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("PLOT TOP PROTEIN DOMAINS")
print("------------------------------------------------------\n")

# find total number of variants
domain_enrichment["Total"] = domain_enrichment["Oncogenic"] + domain_enrichment["Likely Neutral"]

# extract the domains with the highest number of variants 
high_count_domains = domain_enrichment.sort_values("Total", ascending=False).head(20) 

print("Plotting top protein domains by variant count...\n")

high_count_domains[["Oncogenic", "Likely Neutral"]].plot(
  kind="bar",
  stacked=True,
  figsize=(8,5),
  edgecolor="0.1",
  linewidth=0.3,
  color={"Oncogenic": "#C4473B", "Likely Neutral": "#7e8aa2"}
)

plt.title("Top Protein Domains by Variant Count\nOncogenic vs Neutral Distribution")
plt.xlabel("Domain Name")
plt.ylabel("Number of Variants")
plt.xticks(rotation=45, 
           ha="right")
plt.legend(title="Oncogenicity")
plt.tight_layout()
plt.savefig("explore_cancer_variants/plots/oncogenic_neutral_counts.png")

plt.show()

print("Plotting complete! Plot saved as 'explore_cancer_variants/plots/oncogenic_neutral_counts.png'")

# ------------------------------------------------------------
# Plot top domains enriched for oncogenic variants 
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("PLOT TOP PROTEIN DOMAINS ENRICHED FOR ONCOGENIC VARIANTS")
print("------------------------------------------------------\n")

print("Plotting top protein domains enriched for oncogenic variants...\n")

plt.figure(figsize=(8,5))
sns.barplot(
  data=domain_enrichment.head(15),
  x=domain_enrichment.head(15).index, 
  y="Onco_Neutral_Ratio",
  color="#C4473B",
  edgecolor="0.1",
  linewidth=0.3
)

plt.title("Domains Enriched for Oncogenic Variants \nOncogenic to Neutral Ratio", fontsize=14, pad=10)
plt.xlabel("Domain Name", fontsize=12)
plt.ylabel("Oncogenic-to-Neutral Ratio", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)


plt.tight_layout()
plt.savefig("explore_cancer_variants/plots/domain_oncogenic_enrichment.png")
plt.show()

print("Plotting complete! Plot saved as 'explore_cancer_variants/plots/domain_oncogenic_enrichment.png'")

# ------------------------------------------------------------
# Combined oncogenic + neutral heatmap 
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("COMBINED HEATMAP (TOP DOMAINS x TOP GENES)")
print("------------------------------------------------------\n")

# count oncogenic + neutral variants per domain and gene 
gene_domain_class_counts = (
  variants_domains
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

# compute oncogenic fraction 
gene_domain_matrix["Total"] = (
  gene_domain_matrix["Oncogenic"] + gene_domain_matrix["Likely Neutral"] 
)

gene_domain_matrix["Oncogenic_Fraction"] = (
  gene_domain_matrix["Oncogenic"] / gene_domain_matrix["Total"]
).fillna(0) 

# make sure top domain list is correct
top_domains = (
    domain_enrichment.sort_values("Total", ascending=False)
    .head(20)
    .index
)

# filter to top domains
combined_top = gene_domain_matrix[
    gene_domain_matrix["DOMAIN_NAME"].isin(top_domains)
]

# select top genes
top_genes_combined = (
    combined_top.groupby("Hugo_Symbol")["Total"]
    .sum()
    .sort_values(ascending=False)
    .head(20)
    .index
)

combined_top = combined_top[combined_top["Hugo_Symbol"].isin(top_genes_combined)]

# pivot for heatmap
heatmap_combined = combined_top.pivot(
    index="Hugo_Symbol",
    columns="DOMAIN_NAME",
    values="Oncogenic_Fraction"
).fillna(0)

print("Creating heatmap of top domains x top genes (neutral + oncogenic variants)...\n")

plt.figure(figsize=(7,6))
sns.heatmap(heatmap_combined, 
            cmap="Reds", 
            vmin=0, 
            vmax=1, 
            linewidths=0.2)
plt.title("Fraction of Oncogenic Variants per Gene × Domain", fontsize=14, pad=12)
plt.xlabel("Domain Name", fontsize=12)
plt.ylabel("Gene (Hugo Symbol)", fontsize=12)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)

plt.tight_layout()
plt.savefig("explore_cancer_variants/plots/heatmap_oncogenic_fraction.png",
            dpi=300, 
            bbox_inches="tight")
plt.show()

print("Heatmap complete! Saved as 'explore_cancer_variants/plots/heatmap_oncogenic_fraction.png'")

# ------------------------------------------------------------
# Identify driver genes enriched in protein domains
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("ONCOGENIC DRIVER GENES ENRICHED IN PROTEIN DOMAINS")
print("------------------------------------------------------\n")

print("Identifying oncogenic driver genes enriched in protein domains..\n")

# extract oncogenic variants from the exploded df 
oncogenic_variants = variants_domains[variants_domains["ONCOGENIC"] == "Oncogenic"].copy()

# count oncogenic variants per gene x domain 
# "How many oncogenic variants does each gene have in each domain?"
gene_domain_counts = (
    oncogenic_variants.groupby(["Hugo_Symbol", "DOMAIN_NAME"])
    .size()
    .reset_index(name="Variant_Count")
    .sort_values("Variant_Count", ascending=False)
)

# compute total oncogenic variants per domain 
domain_oncogenic_totals = (
    gene_domain_counts.groupby("DOMAIN_NAME")["Variant_Count"]
    .sum()
    .reset_index(name="Domain_Total")
)

# compute the fraction contributed per gene 
gene_domain_fraction = gene_domain_counts.merge(
  domain_oncogenic_totals, on="DOMAIN_NAME")

gene_domain_fraction["Fraction_of_Domain"] = (
    gene_domain_fraction["Variant_Count"] / gene_domain_fraction["Domain_Total"]
)

print("Preview of oncogenic driver genes:\n")
print(gene_domain_fraction.head(5), "\n")

# ------------------------------------------------------------
# Top domains x top genes 
# (Based on oncogenic variants only)
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("ONCOGENIC HEATMAP (TOP DOMAINS x TOP GENES)")
print("------------------------------------------------------\n")

print("Plotting heatmap of top domains x top genes (oncogenic variants)...\n")

n_genes = 20
n_domains = 20

# top domains by total oncogenic variants 
top_domains = (
  gene_domain_fraction.groupby("DOMAIN_NAME")["Variant_Count"]
  .sum()
  .sort_values(ascending=False) 
  .head(n_domains)
  .index
)

# extract only top domains
gene_domain_fraction_top = gene_domain_fraction[gene_domain_fraction["DOMAIN_NAME"].isin(top_domains)]

# pick top genes across top domains 
top_genes = (
  gene_domain_fraction_top.groupby("Hugo_Symbol")["Variant_Count"]
  .sum() 
  .sort_values(ascending=False) 
  .head(n_genes) 
  .index
)

gene_domain_fraction_top = gene_domain_fraction_top[gene_domain_fraction_top["Hugo_Symbol"].isin(top_genes)]

# pivot for heatmap 
heatmap_df = gene_domain_fraction_top.pivot(
  index="Hugo_Symbol",
  columns="DOMAIN_NAME",
  values="Fraction_of_Domain"
).fillna(0) 

plt.figure(figsize=(7,6))
sns.heatmap(heatmap_df, cmap="Reds", linewidths=0.2)

plt.title("Enrichment of Oncogenic Variants \n(Top genes x Top domains)", fontsize=14, pad=12)
plt.xlabel("Domain Name", fontsize=12)
plt.ylabel("Gene (Hugo Symbol)", fontsize=12)
plt.xticks(fontsize=9)
plt.yticks(fontsize=9)

plt.tight_layout()
plt.savefig("explore_cancer_variants/plots/heatmap_topgenes_topdomains.png",
            dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'explore_cancer_variants/plots/heatmap_topgenes_topdomains.png'\n")

print("========================================================")
print("Variant protein domain analysis complete!\n")
