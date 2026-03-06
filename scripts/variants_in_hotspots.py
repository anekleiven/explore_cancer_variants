"""
====================================================================
Variant Hotspots Analysis Script
====================================================================

Script: variants_in_hotspots.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across different cancer hotspots 

Major outputs:
--------------
1. Overview of variant distribution inside and outside cancer hotspots
2. Number of oncogenic and likely neutral variants in cancer hotspots
3. Fraction of oncogenic and likely neutral variants in cancer hotspots
4. Bar plot of variant counts in cancer hotspots (oncogenic vs. neutral)
5. Bar plot of fractions of variants in cancer hotspots (oncogenic vs. neutral)
6. Identification of genes with recurrent oncogenic variants in cancer hotspots
7. Gene-level summary of oncogenic hotspot enrichment
8. Visualization of frequently mutated genes with high hotspot fractions
9. Overview and plot of variants meeting  ClinGen/CGC/VICC hotspot criteria

All plots are saved in:
   plots/

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

# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------

print("\nLoading variant data..")

variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.")

# ------------------------------------------------------------
# Define hotspot membership for each variant
# ------------------------------------------------------------

# create boolean, check if variants is in hotspot true/false
variants["In_Hotspot"] = variants["Hotspot_Type"].notna()

# ------------------------------------------------------------
# Overview of the data 
# ------------------------------------------------------------

# total number of variants in data set
total_num = len(variants)

print("-"*30)
print("Variant Overview:")

# number of variants in cancer hotspots
variants_in_hotspots = variants["In_Hotspot"].sum()
print(f"Found {variants_in_hotspots:,} variants in cancer hotspots.")

# number of variants not in cancer hotspots 
variants_not_in_hotspots = (~variants["In_Hotspot"]).sum() 
print(f"Found {variants_not_in_hotspots:,} variants outside cancer hotspots.")

# fraction of variants in cancer hotspots 
fraction_in_hotspots = variants_in_hotspots / total_num * 100 
print(f"{fraction_in_hotspots:.2f}% of the variants in the data is inside cancer hotspots.")
print("-"*30)


# ------------------------------------------------------------
# Fraction of Oncogenic Variants in Cancer Hotspots 
# ------------------------------------------------------------

# extract oncogenic and likely neutral variants 
oncogenic = variants[variants["ONCOGENIC"] == "Oncogenic"]
likely_neutral = variants[variants["ONCOGENIC"] == "Likely Neutral"]

print("Oncogenic Variants:")

# find the number of oncogenic variants in cancer hotspots 
oncogenic_in_hotspots = oncogenic["In_Hotspot"].sum() 
print(f"Found {oncogenic_in_hotspots:,} oncogenic variants in cancer hotspots.")

# find the number of oncogenic variants outside cancer hotspots 
oncogenic_not_in_hotspots = (~oncogenic["In_Hotspot"]).sum() 
print(f"Found {oncogenic_not_in_hotspots:,} oncogenic variants outside cancer hotspots.")

# fraction of oncogenic variants in cancer hotspots 
fraction_oncogenic_in_hotspots = oncogenic_in_hotspots / len(oncogenic) * 100 
print(f"{fraction_oncogenic_in_hotspots:.2f}% of oncogenic variants found inside cancer hotspots.")

# of all variants inside cancer hotspots, what fraction is oncogenic? 
fraction_oncogenic = oncogenic_in_hotspots / variants_in_hotspots * 100 
print(f"{fraction_oncogenic:.2f} % of all variants in cancer hotspots are oncogenic.\n")

# ------------------------------------------------------------
# Fraction of Likely Neutral Variants in Cancer Hotspots 
# ------------------------------------------------------------

print("-"*30)
print("Likely Neutral Variants:")

# find the number of likely neutral variants in cancer hotspots 
neutral_in_hotspots = likely_neutral["In_Hotspot"].sum() 
print(f"Found {neutral_in_hotspots:,} likely neutral variants in cancer hotspots.")

# find the number of likely neutral variants outside cancer hotspots 
neutral_not_in_hotspots = (~likely_neutral["In_Hotspot"]).sum() 
print(f"Found {neutral_not_in_hotspots:,} likely neutral variants outside cancer hotspots.")

# fraction of likely neutral variants in cancer hotspots 
fraction_neutral_in_hotspots = neutral_in_hotspots / len(likely_neutral) * 100 
print(f"{fraction_neutral_in_hotspots:.2f} % of likely neutral variants found inside cancer hotspots.")

# of all variants inside cancer hotspots, what fraction is likely neutral? 
fraction_neutral = neutral_in_hotspots / variants_in_hotspots * 100 
print(f"{fraction_neutral:.2f} % of all variants in cancer hotspots are likely neutral.\n")

# ------------------------------------------------------------
# Plot number of variants in cancer hotspots (oncogenic vs neutral)
# ------------------------------------------------------------

print("-"*30)
print("Plotting number of variants in cancer hotspots...\n")

# keep only neutral and oncogenic variants
classes = ["Oncogenic", "Likely Neutral"]
variants_onco_neutral = variants[variants["ONCOGENIC"].isin(classes)]

# group variants by oncogenicity and cancer hotspots 
counts = (variants_onco_neutral
          .groupby(["In_Hotspot", "ONCOGENIC"])
          .size() 
          .reset_index(name="Variant_Count") 
)

palette = {"Oncogenic": "#C4473B", "Likely Neutral": "#7e8aa2"}

plt.figure(figsize=(8,5))
sns.barplot(data=counts, 
            x="In_Hotspot",
            y="Variant_Count",
            hue="ONCOGENIC", 
            palette=palette, 
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Number of Variants in Cancer Hotspots", fontsize=14, pad=10)
plt.xlabel("Variant in Hotspot", fontsize=12)
plt.ylabel("Counts", fontsize=12) 
plt.tight_layout() 
plt.savefig("plots/variants_in_hotspots.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'plots/variants_in_hotspots.png'")

# ------------------------------------------------------------
# Plot fraction of oncogenic vs neutral variants in cancer hotspots 
# ------------------------------------------------------------

print("-"*30)
print("Computing fraction of variants in cancer hotspots (oncogenic vs. neutral)..")

totals = variants_onco_neutral["ONCOGENIC"].value_counts().rename("Total")

counts = counts.merge(totals, left_on="ONCOGENIC", right_index=True)
counts["Fraction"] = counts["Variant_Count"] / counts["Total"] 

print("Fraction of variants in each feature type per class:")
print(counts)
print("-"*30)

print("Plotting fraction of variants in cancer hotspots (oncogenic vs. neutral)..")

plt.figure(figsize=(8,5)) 
sns.barplot(data=counts, 
            x="In_Hotspot",
            y="Fraction",
            hue="ONCOGENIC",
            palette=palette,
            edgecolor="0.1",
            linewidth=0.3)

plt.title("Fraction of Variants in Cancer Hotspots", fontsize=14, pad=10) 
plt.xlabel("Variant in Hotspot", fontsize=12)
plt.ylabel("Fraction", fontsize=12)
plt.tight_layout()
plt.savefig("plots/fractions_in_hotspots.png", dpi=300)
plt.show() 

print("Plotting complete! Plot saved as 'plots/fractions_in_hotspots.png'\n")

# ------------------------------------------------------------
# Identify Oncogenic Variants in Cancer Hotspots Across Genes 
# ------------------------------------------------------------

print("-"*30)
print("Identifying oncogenic driver genes in cancer hotspots..")

oncogenic = variants_onco_neutral[variants_onco_neutral["ONCOGENIC"] == "Oncogenic"] 

oncogenic_hotspots = oncogenic[oncogenic["In_Hotspot"] == True]

onco_genes = (oncogenic_hotspots
              .groupby("Hugo_Symbol")
              .size() 
              .reset_index(name="Hotspot_Variant_Count")
              .sort_values("Hotspot_Variant_Count", ascending=False)
)

print("Example of genes with a high number of variants in cancer hotspots:")
print(onco_genes.head(), "\n")

# ------------------------------------------------------------
# Plot Oncogenic Variants Across Genes 
# ------------------------------------------------------------

print("-"*30)
print("Plotting Oncogenic Variants in Cancer Hotspots across Genes...\n")

top_oncogenes = onco_genes.head(20) 

plt.figure(figsize=(8,5))
sns.barplot(data=top_oncogenes,
            x="Hugo_Symbol",
            y="Hotspot_Variant_Count",
            color="#C4473B",
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Top Oncogenic Genes in Cancer Hotspots", fontsize=14, pad=10) 
plt.xlabel("Hugo Symbol", fontsize=12) 
plt.ylabel("Number of Variants", fontsize=12) 
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9) 

plt.tight_layout()
plt.savefig("plots/oncogenes_in_hotspots.png", dpi=300)
plt.show() 

print("Plotting complete! Plot saved as 'plots/oncogenes_in_hotspots.png'\n")

# ------------------------------------------------------------
# Gene-level Hotspot Fraction 
# ------------------------------------------------------------

# Find out whether the genes are hotspot-driven or just frequently mutated 

print("-"*30)
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

print("-"*30)
print("Plotting frequently mutated genes with hotspot enrichment..")

top_genes = oncogenic_gene_summary[
  (oncogenic_gene_summary["Total_Oncogenic"] >= 10) &
  (oncogenic_gene_summary["Hotspot_Fraction"] >= 0.3)
].sort_values("Total_Oncogenic", ascending=False) 

print("Example output frequently mutated genes:")
print(top_genes.head(5),"\n")

plt.figure(figsize=(10,6)) 
sns.barplot(
  data=top_genes.head(15), 
  x="Hugo_Symbol",
  y="Total_Oncogenic",
  hue="Hotspot_Fraction",
  palette="YlOrRd",
  dodge=False,
  edgecolor="0.1",
  linewidth=0.3
)

plt.title("Frequently Mutated Genes with Hotspot Enrichment", fontsize=14, pad=10)
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Total Oncogenic Mutations", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Hotspot Fraction", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("plots/oncogenes_hotspot_fraction.png", dpi=300)
plt.show()

print("Plotting complete! Figure saved as 'plots/oncogenes_hotspot_fraction.png'\n")

# ------------------------------------------------------------
# Cancer Hotspot Evidence (ClinGen / CGC / VICC Framework)
# Based on Horak et al.
# ------------------------------------------------------------

print("-"*30)
print("Identifying variants meeting ClinGen/CGC/VICC cancer hotspot criteria (OS3, OM3)..\n")


# Define hotspot-related evidence components

# Variant occurs at a hotspot position with ≥50 somatic samples
variants_onco_neutral["Hotspot_Pos_50"] = variants_onco_neutral["Samples"] >= 50

# Count how many hotspot samples share the exact amino acid change
variants_onco_neutral["Exact_AA_Hotspot_Count"] = (
    variants_onco_neutral
    .groupby("HGVSp")["In_Hotspot"]
    .transform("sum")
)

# Variant has ≥10 hotspot samples with the same amino acid change
variants_onco_neutral["Exact_AA_10"] = (
    variants_onco_neutral["Exact_AA_Hotspot_Count"] >= 10
)

# Apply ClinGen/VICC cancer hotspot criteria

# OS3: Hotspot + ≥50 samples at position + ≥10 with exact AA change
variants_onco_neutral["Meets_Hotspot_OS3"] = (
    variants_onco_neutral["In_Hotspot"] &
    variants_onco_neutral["Hotspot_Pos_50"] &
    variants_onco_neutral["Exact_AA_10"]
)

# OM3: Hotspot + ≥10 samples with exact AA change
variants_onco_neutral["Meets_Hotspot_OM3"] = (
    variants_onco_neutral["In_Hotspot"] &
    variants_onco_neutral["Exact_AA_10"]
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
print("-"*30)

# ------------------------------------------------------------
# Plot Cancer Hotspot Evidence
# ------------------------------------------------------------

print("Plotting cancer hotspot evidence..")

hotspot_palette = {True: "#B53226", False: "#A5BAE4"}


# OS3 PLOT
plt.figure(figsize=(8, 5))

sns.barplot(
    data=summary_OS3,
    x="ONCOGENIC",
    y="Variant_Count",
    hue="Meets_Hotspot_OS3",
    palette=hotspot_palette,
    edgecolor="0.1",
    linewidth=0.3
)

plt.title("ClinGen/CGC/VICC Cancer Hotspot Evidence (OS3)", fontsize=14, pad=10)
plt.xlabel("Oncogenicity", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12)
plt.legend(title="Meets OS3 Criterion", bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()
plt.savefig("plots/meets_hotspot_OS3.png", dpi=300)
plt.show()

print("OS3 plot saved as 'plots/meets_hotspot_OS3.png'")


# OM3 PLOT
plt.figure(figsize=(8, 5))

sns.barplot(
    data=summary_OM3,
    x="ONCOGENIC",
    y="Variant_Count",
    hue="Meets_Hotspot_OM3",
    palette=hotspot_palette,
    edgecolor="0.1",
    linewidth=0.3
)

plt.title("ClinGen/CGC/VICC Cancer Hotspot Evidence (OM3)", fontsize=14, pad=10)
plt.xlabel("Oncogenicity", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12)
plt.legend(title="Meets OM3 Criterion", bbox_to_anchor=(1.05, 1), loc="upper left")

plt.tight_layout()
plt.savefig("plots/meets_hotspot_OM3.png", dpi=300)
plt.show()

print("OM3 plot saved as 'plots/meets_hotspot_OM3.png'")

print("\nCancer hotspot analysis complete!🥳🥳\n")