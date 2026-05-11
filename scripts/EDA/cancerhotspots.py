
# ====================================================================
# Cancer Hotspots Analysis
# ====================================================================

"""
Script: cancerhotspots.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants (oncogenic and likely neutral) distribute across different cancer hotspots 

Major outputs:
--------------
1. Overview of the cancer hotspot data 
2. Barplot of variant counts in cancer hotspots 
3. Barplot of variant fractions cancer hotspots
4. Genes with recurrent oncogenic variants in cancer hotspots
5. Gene-level summary of oncogenic hotspot enrichment
6. Visualization of frequently mutated genes with high hotspot fractions
7. Overview and plot of variants meeting ClinGen/CGC/VICC hotspot criteria (OS3, OM3)

All plots are saved in:
   plots/cancerhotspots 

"""

print("-"*50)
print("Variant Hotspot Analysis🤓")
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
        description="Explore cancer hotspots in variant data."
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

save_dir = "plots/cancerhotspots"
os.makedirs(save_dir, exist_ok=True) 


# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------

print(f"\nLoading variant file:\n{args.variants}")

variants = pd.read_csv(
    args.variants,
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.")

# ------------------------------------------------------------
# Overview of the data 
# ------------------------------------------------------------

# Total number of variants in data set
total_num = len(variants)

print("-"*50)
print("Variant Overview:")

# Number of variants in cancer hotspots
variants_in_hotspots = variants["In_Hotspot"].sum()
print(f"Found {variants_in_hotspots:,} variants in cancer hotspots.")

# Number of variants not in cancer hotspots 
variants_not_in_hotspots = (~variants["In_Hotspot"]).sum() 
print(f"Found {variants_not_in_hotspots:,} variants outside cancer hotspots.")

# Fraction of variants in cancer hotspots 
fraction_in_hotspots = variants_in_hotspots / total_num * 100 
print(f"{fraction_in_hotspots:.2f}% of the variants in the data is inside cancer hotspots.")
print("-"*50)


# ------------------------------------------------------------
# Fraction of Oncogenic Variants in Cancer Hotspots 
# ------------------------------------------------------------

# Extract oncogenic and likely neutral variants 
oncogenic = variants[variants["ONCOGENIC"] == "Oncogenic"]
likely_neutral = variants[variants["ONCOGENIC"] == "Likely Neutral"]

print("Oncogenic Variants:")

# Find the number of oncogenic variants in cancer hotspots 
oncogenic_in_hotspots = oncogenic["In_Hotspot"].sum() 
print(f"Found {oncogenic_in_hotspots:,} oncogenic variants in cancer hotspots.")

# Find the number of oncogenic variants outside cancer hotspots 
oncogenic_not_in_hotspots = (~oncogenic["In_Hotspot"]).sum() 
print(f"Found {oncogenic_not_in_hotspots:,} oncogenic variants outside cancer hotspots.")

# Fraction of oncogenic variants in cancer hotspots 
fraction_oncogenic_in_hotspots = oncogenic_in_hotspots / len(oncogenic) * 100 
print(f"{fraction_oncogenic_in_hotspots:.2f}% of oncogenic variants found inside cancer hotspots.")

# of all variants inside cancer hotspots, what fraction is oncogenic? 
fraction_oncogenic = oncogenic_in_hotspots / variants_in_hotspots * 100 
print(f"{fraction_oncogenic:.2f} % of all variants in cancer hotspots are oncogenic.\n")

# ------------------------------------------------------------
# Fraction of Likely Neutral Variants in Cancer Hotspots 
# ------------------------------------------------------------

print("-"*50)
print("Likely Neutral Variants:")

# Find the number of likely neutral variants in cancer hotspots 
neutral_in_hotspots = likely_neutral["In_Hotspot"].sum() 
print(f"Found {neutral_in_hotspots:,} likely neutral variants in cancer hotspots.")

# Find the number of likely neutral variants outside cancer hotspots 
neutral_not_in_hotspots = (~likely_neutral["In_Hotspot"]).sum() 
print(f"Found {neutral_not_in_hotspots:,} likely neutral variants outside cancer hotspots.")

# Fraction of likely neutral variants in cancer hotspots 
fraction_neutral_in_hotspots = neutral_in_hotspots / len(likely_neutral) * 100 
print(f"{fraction_neutral_in_hotspots:.2f} % of likely neutral variants found inside cancer hotspots.")

# Of all variants inside cancer hotspots, what fraction is likely neutral? 
fraction_neutral = neutral_in_hotspots / variants_in_hotspots * 100 
print(f"{fraction_neutral:.2f} % of all variants in cancer hotspots are likely neutral.\n")

# ------------------------------------------------------------
# Plot number of variants in cancer hotspots (oncogenic vs. neutral)
# ------------------------------------------------------------

print("-"*50)
print("Plotting number of variants in cancer hotspots...\n")

# Keep only neutral and oncogenic variants
classes = ["Oncogenic", "Likely Neutral"]
variants_onco_neutral = variants[variants["ONCOGENIC"].isin(classes)]

# Group variants by oncogenicity and cancer hotspots 
counts = (variants_onco_neutral
          .groupby(["In_Hotspot", "ONCOGENIC"])
          .size() 
          .reset_index(name="Variant_Count") 
)

palette = {"Oncogenic": "#c4314a", "Likely Neutral": "#88aed1"}

plt.figure(figsize=(8,5))

sns.set_style("white")

counts_plot = counts.copy() 
counts_plot['In_Hotspot'] = counts_plot['In_Hotspot'].map({True: 'In Hotspot', False: 'Outside Hotspot'})

sns.barplot(data=counts_plot, 
            x="In_Hotspot",
            y="Variant_Count",
            order=['In Hotspot', 'Outside Hotspot'],
            hue="ONCOGENIC", 
            palette=palette, 
            edgecolor="black",
            linewidth=0.8) 

sns.despine() 

plt.title("Variants in Cancer Mutation Hotspots", fontsize=14, pad=15)
plt.xlabel("Hotspot Status", fontsize=12, labelpad=10)
plt.ylabel("Counts", fontsize=12, labelpad=10) 
plt.xticks(rotation=0, fontsize=9)
plt.legend(title="Oncogenicity Class", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout() 
plt.savefig(f"{save_dir}/variants_in_hotspots.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/variants_in_hotspots.png'")

# ------------------------------------------------------------
# Plot fraction of oncogenic vs. neutral variants in cancer hotspots 
# ------------------------------------------------------------

print("-"*50)
print("Computing fraction of variants in cancer hotspots (oncogenic vs. neutral)..")

totals = variants_onco_neutral["ONCOGENIC"].value_counts().rename("Total")

counts = counts.merge(totals, left_on="ONCOGENIC", right_index=True)
counts["Fraction"] = counts["Variant_Count"] / counts["Total"] 

print("Fraction of variants in each feature type per class:")
print(counts)
print("-"*50)

print("Plotting fraction of variants in cancer hotspots (oncogenic vs. neutral)..")

plt.figure(figsize=(8,5)) 

sns.set_style("white")

counts_plot2 = counts.copy()
counts_plot2['In_Hotspot'] = counts_plot2['In_Hotspot'].map({True: 'In Hotspot', False: 'Outside Hotspot'})

sns.barplot(data=counts_plot2, 
            x="In_Hotspot",
            y="Fraction",
            order=['In Hotspot', 'Outside Hotspot'],
            hue="ONCOGENIC",
            palette=palette,
            edgecolor="black",
            linewidth=0.8)

sns.despine() 

plt.title("Fraction of Variants in Cancer Mutation Hotspots", fontsize=14, pad=15) 
plt.xlabel("Hotspot Status", fontsize=12, labelpad=10)
plt.ylabel("Fraction", fontsize=12, labelpad=10)
plt.xticks(rotation=0, fontsize=9)
plt.legend(title="Oncogenicity Class", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/fractions_in_hotspots.png", dpi=300, bbox_inches="tight")
plt.show() 

print(f"Plotting complete! Plot saved as '{save_dir}/fractions_in_hotspots.png'\n")

# ------------------------------------------------------------
# Identify Oncogenic Variants in Cancer Hotspots Across Genes 
# ------------------------------------------------------------

print("-"*50)
print("Identifying oncogenic driver genes in cancer hotspots..")

oncogenic = variants_onco_neutral[variants_onco_neutral["ONCOGENIC"] == "Oncogenic"] 

oncogenic_hotspots = oncogenic[oncogenic["In_Hotspot"] == True]

onco_genes = (oncogenic_hotspots
              .groupby("Hugo_Symbol")
              .size() 
              .reset_index(name="Hotspot_Variant_Count")
              .sort_values("Hotspot_Variant_Count", ascending=False)
)

print("Example of genes with a high number of oncogenic variants in cancer hotspots:")
print(onco_genes.head(), "\n")

# ------------------------------------------------------------
# Plot Oncogenic Variants Across Genes 
# ------------------------------------------------------------

print("-"*50)
print("Plotting Oncogenic Variants in Cancer Hotspots across Genes...\n")

top_oncogenes = onco_genes.head(20) 

plt.figure(figsize=(8,5))
sns.set_style("white")

sns.barplot(data=top_oncogenes,
            x="Hugo_Symbol",
            y="Hotspot_Variant_Count",
            color="#c4314a",
            edgecolor="black",
            linewidth=0.8) 

sns.despine() 

plt.title("Top Genes by Oncogenic Variant Burden in Cancer Mutation Hotspots", fontsize=14, pad=15) 
plt.xlabel("Gene (Hugo Symbol)", fontsize=12, labelpad=10) 
plt.ylabel("Number of Variants", fontsize=12, labelpad=10) 
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9) 

plt.tight_layout()
plt.savefig(f"{save_dir}/oncogenes_in_hotspots.png", dpi=300, bbox_inches="tight")
plt.show() 

print(f"Plotting complete! Plot saved as '{save_dir}/oncogenes_in_hotspots.png'\n")

# ------------------------------------------------------------
# Gene-level Hotspot Fraction 
# ------------------------------------------------------------

# Find out whether the genes are hotspot-driven or just frequently mutated 

print("-"*50)
print("Calculate the fraction of oncogenic variants in cancer hotspots for highly mutated genes..")

gene_totals = (
  oncogenic
  .groupby("Hugo_Symbol")
  .size() 
  .rename("Total_Oncogenic") 
)

gene_hotspots = (
  oncogenic_hotspots
  .groupby("Hugo_Symbol")
  .size() 
  .rename("Hotspot_Oncogenic") 
)

oncogenic_gene_summary = pd.concat([gene_totals, gene_hotspots], axis = 1).fillna(0)

oncogenic_gene_summary["Hotspot_Fraction"] = (
  oncogenic_gene_summary["Hotspot_Oncogenic"] / oncogenic_gene_summary["Total_Oncogenic"]
)

oncogenic_gene_summary = oncogenic_gene_summary.sort_values(
  "Total_Oncogenic", ascending=False
  )

print("Gene-level oncogenic hotspot summary:")
print(oncogenic_gene_summary.head(10))

# ------------------------------------------------------------
# Plot frequently mutated genes with hotspot enrichment 
# ------------------------------------------------------------

print("-"*50)
print("Plotting frequently mutated genes with hotspot enrichment..")

top_genes = oncogenic_gene_summary[
  (oncogenic_gene_summary["Total_Oncogenic"] >= 10) &
  (oncogenic_gene_summary["Hotspot_Fraction"] >= 0.3)
].sort_values("Total_Oncogenic", ascending=False) 

print("Example output frequently mutated genes:")
print(top_genes.head(5),"\n")

plt.figure(figsize=(8,5)) 

sns.set_style("white")

sns.barplot(
  data=top_genes.head(15), 
  x="Hugo_Symbol",
  y="Total_Oncogenic",
  hue="Hotspot_Fraction",
  palette="YlOrRd",
  edgecolor="black",
  linewidth=0.5
)

sns.despine()

plt.title("Genes with Cancer Mutation Hotspot Enrichment", fontsize=14, pad=15)
plt.xlabel("Gene", fontsize=12, labelpad=10)
plt.ylabel("Number of Variants (Oncogenic)", fontsize=12, labelpad=10)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Hotspot Fraction", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/oncogenes_hotspot_fraction.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Figure saved as '{save_dir}/oncogenes_hotspot_fraction.png'\n")

# ------------------------------------------------------------
# Cancer Hotspot Evidence (ClinGen / CGC / VICC Framework)
# Based on Horak et al.
# ------------------------------------------------------------

print("-"*50)
print("Identifying variants meeting ClinGen/CGC/VICC cancer hotspot criteria (OS3, OM3)..\n")

# DEFINE HOTSPOT EVIDENCE

# Count number of samples with a variant with the same amino acid position 
variants_onco_neutral["Pos_Total_Samples"] = (
    variants_onco_neutral
    .groupby(["Hugo_Symbol", "Protein_position"])["Samples"]
    .transform("sum")
)

# Count number of samples with a variant with the same amino acid change 
variants_onco_neutral["Exact_AA_Count"] = variants_onco_neutral["Samples"]

# Apply ClinGen/VICC cancer hotspot criteria

# OS3: Hotspot + >= 50 on protein position + >= 10 on exact change 
variants_onco_neutral["Meets_Hotspot_OS3"] = (
    variants_onco_neutral["In_Hotspot"] & 
    (variants_onco_neutral["Pos_Total_Samples"] >= 50) & 
    (variants_onco_neutral["Exact_AA_Count"] >= 10)
)

# OM3: Hotspot + < 50 on protein position + >= 10 on exact change 
variants_onco_neutral["Meets_Hotspot_OM3"] = (
    variants_onco_neutral["In_Hotspot"] & 
    (variants_onco_neutral["Pos_Total_Samples"] < 50) & 
    (variants_onco_neutral["Exact_AA_Count"] >= 10)
)

print("Preview of annotated variants:")
print(variants_onco_neutral.head(),"\n")


# ------------------------------------------------------------
# Summarize variant counts by oncogenicity
# ------------------------------------------------------------

summary_OS3 = (
    variants_onco_neutral
    .groupby(["ONCOGENIC", "Meets_Hotspot_OS3"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Summary of Cancer Hotspot Criterion OS3:")
print(summary_OS3, "\n")

summary_OM3 = (
    variants_onco_neutral
    .groupby(["ONCOGENIC", "Meets_Hotspot_OM3"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Summary of Cancer Hotspot Criterion OM3:")
print(summary_OM3)
print("-"*50)

# ------------------------------------------------------------
# Plot Cancer Hotspot Evidence
# ------------------------------------------------------------

print("Plotting cancer hotspot evidence..")

sns.set_style("white")

# ----------------
# OS3 PLOT
# ----------------

plt.figure(figsize=(8, 5))

ax = sns.barplot(
    data=summary_OS3,
    x="ONCOGENIC",
    y="Variant_Count",
    order=['Oncogenic', 'Likely Neutral'],
    hue="Meets_Hotspot_OS3",
    palette={True: "#f0483c", False: "#d9d7d7"},
    edgecolor="black",
    linewidth=0.8
)

sns.despine()

plt.title("ClinGen/CGC/VICC Cancer Hotspot Evidence (OS3)", fontsize=14, pad=15)
plt.xlabel("Oncogenicity Class", fontsize=12, labelpad=10)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10)

plt.legend(title="Hotspot Criterion Met", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/meets_hotspot_OS3.png", dpi=300, bbox_inches="tight")
plt.show()
print(f"OS3 plot saved as '{save_dir}/meets_hotspot_OS3.png'")


# ----------------
# OM3 PLOT
# ----------------

plt.figure(figsize=(8, 5))

ax = sns.barplot(
    data=summary_OM3,
    x="ONCOGENIC",
    y="Variant_Count",
    hue="Meets_Hotspot_OM3",
    order=['Oncogenic', 'Likely Neutral'],
    palette={True: "#ff8c45", False: "#d9d7d7"},
    edgecolor="black",
    linewidth=0.8
)

sns.despine()

plt.title("ClinGen/CGC/VICC Cancer Hotspot Evidence (OM3)", fontsize=14, pad=15)
plt.xlabel("Oncogenicity Class", fontsize=12, labelpad=10)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10)

plt.legend(title="Hotspot Criterion Met", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/meets_hotspot_OM3.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"OM3 plot saved as '{save_dir}/meets_hotspot_OM3.png'")
print("-"*50)

print("\nCancer hotspot analysis complete!🥳🥳\n")