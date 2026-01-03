"""
====================================================================
Variant Protein Domain Analysis Script
====================================================================
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
    visualize_variants/plots/

"""

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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

# ------------------------------------------------------------
# Count variants inside and outside protein domains
# ------------------------------------------------------------

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

print("\nONCOGENIC")
print(count_in_out(oncogenic, "Oncogenic"))

print("\nLIKELY NEUTRAL")
print(count_in_out(neutral, "Likely Neutral")) 


# ------------------------------------------------------------
# Expand DOMAIN_NAME for variants with multiple domains 
# ------------------------------------------------------------

print("\nExploding variants with multiple domains...\n")

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

print("\nComputing oncogenic vs neutral enrichment per domain...\n")

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
print(domain_enrichment.head(15), "\n")

# ------------------------------------------------------------
# Plot top domains by number of total variants 
# ------------------------------------------------------------

# find total number of variants
domain_enrichment["Total"] = domain_enrichment["Oncogenic"] + domain_enrichment["Likely Neutral"]

# extract the domains with the highest number of variants 
high_count_domains = domain_enrichment.sort_values("Total", ascending=False).head(20) 

high_count_domains[["Oncogenic", "Likely Neutral"]].plot(
  kind="bar",
  stacked=True,
  figsize=(10,6),
  edgecolor="black",
  linewidth=0.2,
  color={"Oncogenic": "#ef6f6c", "Likely Neutral": "#5b8fdc"}
)

plt.title("Top Protein Domains by Variant Count\nOncogenic vs Neutral Distribution",
          fontsize=16, 
          fontweight="bold")
plt.xlabel("Domain Name", 
           fontsize=13)
plt.ylabel("Number of Variants", 
           fontsize=13)
plt.xticks(rotation=45, 
           ha="right")
plt.legend(title="Variant Class")

plt.tight_layout()
plt.savefig("visualize_variants/plots/oncogenic_neutral_counts.png")

plt.show()

# ------------------------------------------------------------
# Plot top domains enriched for oncogenic variants 
# ------------------------------------------------------------

plt.figure(figsize=(10,6))
sns.barplot(
  data=domain_enrichment.head(15),
  x=domain_enrichment.head(15).index, 
  y="Onco_Neutral_Ratio",
  color="#ef6f6c"
)

plt.xticks(rotation=45, 
           ha="right")
plt.xlabel("Domain Name",
           fontsize=13)
plt.ylabel("Oncogenic-to-Neutral Ratio", 
           fontsize=13)
plt.title("Domains Enriched for Oncogenic Variants \nOncogenic to Neutral Ratio",
          fontsize=16, 
          fontweight="bold")

plt.tight_layout()
plt.savefig("visualize_variants/plots/domain_oncogenic_enrichment.png")

plt.show()


# ------------------------------------------------------------
# Combined oncogenic + neutral heatmap 
# ------------------------------------------------------------

# count oncogenic + neutral variants per domain and gene 
gene_domain_class_counts = (
  variants_domains
  .groupby(["Hugo_Symbol", "DOMAIN_NAME", "ONCOGENIC"])
  .size() 
  .reset_index(name="Count") 
)

print(gene_domain_class_counts.head())

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

plt.figure(figsize=(12,7))
sns.heatmap(heatmap_combined, 
            cmap="Reds", 
            vmin=0, 
            vmax=1, 
            linewidths=0.2)
plt.title("Fraction of Oncogenic Variants per Gene × Domain", 
          fontsize=14, 
          fontweight="bold")
plt.xlabel("Domain Name", 
           fontsize=14)
plt.ylabel("Gene (Hugo Symbol)",
           fontsize=14)
plt.tight_layout()
plt.savefig("visualize_variants/plots/heatmap_oncogenic_fraction.png",
            dpi=300, 
            bbox_inches="tight")
plt.show()


# ------------------------------------------------------------
# Identify driver genes enriched in protein domains
# ------------------------------------------------------------

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

print("Example output:\n")
print(gene_domain_fraction.head(5), "\n")

# ------------------------------------------------------------
# Top domains x top genes 
# (Based on oncogenic variants only)
# ------------------------------------------------------------

print("\nPlotting heatmap of top domains x top genes...\n")

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

plt.figure(figsize=(10,6))
sns.heatmap(heatmap_df, cmap="Reds", linewidths=0.2)
plt.title("Enrichment of Oncogenic Variants \n(Top genes x Top domains)",
          fontsize=14,
          fontweight="bold"
          )
plt.xlabel("Domain Name",
           fontsize=13)
plt.ylabel("Gene (Hugo Symbol)",
           fontsize=13)

plt.tight_layout()
plt.savefig("visualize_variants/plots/heatmap_topgenes_topdomains.png",
            dpi=300, bbox_inches="tight")
plt.show()