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

All plots are saved in:
    visualize_variants/plots/

"""

print("\n========================================================")
print("VARIANT HOTSPOT ANALYSIS")
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

print("\nLoading variant data...\n")

variants = pd.read_csv(
    "annotation_pipeline/output/variants_with_func_sites.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.\n")


# ------------------------------------------------------------
# Define hotspot membership for each variant
# ------------------------------------------------------------

variants["In_Hotspot"] = variants["Hotspot_Type"].notna()

# ------------------------------------------------------------
# Overview of the data 
# ------------------------------------------------------------

# total number of variants in data set
total_num = len(variants)

print("------------------------------------------------------")
print("Variant Overview:")
print("------------------------------------------------------\n")

# number of variants in cancer hotspots
variants_in_hotspots = variants["In_Hotspot"].sum()
print(f"Found {variants_in_hotspots:,} variants in cancer hotspots.")

# number of variants not in cancer hotspots 
variants_not_in_hotspots = (~variants["In_Hotspot"]).sum() 
print(f"Found {variants_not_in_hotspots:,} variants outside cancer hotspots.")

# fraction of variants in cancer hotspots 
fraction_in_hotspots = variants_in_hotspots / total_num * 100 
print(f"{fraction_in_hotspots:.2f}% of the variants in the data is inside cancer hotspots.\n")


# ------------------------------------------------------------
# Fraction of Oncogenic Variants in Cancer Hotspots 
# ------------------------------------------------------------

# extract oncogenic and likely neutral variants 
oncogenic = variants[variants["ONCOGENIC"] == "Oncogenic"]
likely_neutral = variants[variants["ONCOGENIC"] == "Likely Neutral"]

print("------------------------------------------------------")
print("Oncogenic Variants:")
print("------------------------------------------------------\n")

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

print("------------------------------------------------------")
print("Likely Neutral Variants:")
print("------------------------------------------------------\n")

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

print("------------------------------------------------------")
print("VARIANTS IN CANCER HOTSPOTS PLOT")
print("------------------------------------------------------\n")

# keep only neutral and oncogenic variants
classes = ["Oncogenic", "Likely Neutral"]
variants_onco_neutral = variants[variants["ONCOGENIC"].isin(classes)]

# group variants by oncogenicity and cancer hotspots 
counts = (variants_onco_neutral
          .groupby(["In_Hotspot", "ONCOGENIC"])
          .size() 
          .reset_index(name="Variant_Count") 
)

print("Plotting number of variants in cancer hotspots...\n")

palette = palette = {"Oncogenic": "#ef6f6c", "Likely Neutral": "#5b8fdc"}

plt.figure(figsize=(6,4))
sns.barplot(data=counts, 
            x="In_Hotspot",
            y="Variant_Count",
            hue="ONCOGENIC", 
            palette=palette, 
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Number of Variants in Cancer Hotspots")
plt.ylabel("Counts") 
plt.tight_layout() 
plt.savefig("visualize_variants/plots/variants_in_hotspots.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'visualize_variants/plots/variants_in_hotspots.png'")


# ------------------------------------------------------------
# Plot fraction of oncogenic vs neutral variants in cancer hotspots 
# ------------------------------------------------------------

print("\n------------------------------------------------------")
print("FRACTION OF VARIANTS IN CANCER HOTSPOTS PLOT")
print("------------------------------------------------------\n")

print("Computing fraction of variants in cancer hotspots (oncogenic vs. neutral)...\n")

totals = variants_onco_neutral["ONCOGENIC"].value_counts().rename("Total")

counts = counts.merge(totals, left_on="ONCOGENIC", right_index=True)
counts["Fraction"] = counts["Variant_Count"] / counts["Total"] 

print("Fraction of variants in each feature type per class:\n")
print(counts, "\n")


print("Plotting fraction of variants in cancer hotspots (oncogenic vs. neutral)...\n")

plt.figure(figsize=(6,4)) 
sns.barplot(data=counts, 
            x="In_Hotspot",
            y="Fraction",
            hue="ONCOGENIC",
            palette=palette,
            edgecolor="0.1",
            linewidth=0.3)

plt.title("Fraction of Variants in Cancer Hotspots") 
plt.ylabel("Fraction")
plt.tight_layout()
plt.savefig("visualize_variants/plots/fractions_in_hotspots.png", dpi=300)
plt.show() 

print("\nPlotting complete! Plot saved as 'visualize_variants/plots/fractions_in_hotspots.png'\n")

# ------------------------------------------------------------
# Identify Oncogenic Variants in Cancer Hotspots Across Genes 
# ------------------------------------------------------------

print("------------------------------------------------------")
print("ONCOGENIC VARIANTS IN CANCER HOTSPOTS ACROSS GENES") 
print("------------------------------------------------------\n")

print("Identifying oncogenic driver genes in cancer hotspots...\n")

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

print("------------------------------------------------------")
print("COUNT OF ONCOGENIC VARIANTS PER GENE IN CANCER HOTSPOTS") 
print("------------------------------------------------------\n")
top_oncogenes = onco_genes.head(20) 

print("Plotting Oncogenic Variants in Cancer Hotspots across Genes...\n")

plt.figure(figsize=(6,4))
sns.barplot(data=top_oncogenes,
            x="Hugo_Symbol",
            y="Hotspot_Variant_Count",
            color="#ef6f6c",
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Top Oncogenic Genes in Cancer Hotspots") 
plt.xlabel("Hugo Symbol") 
plt.ylabel("Number of Variants") 
plt.xticks(rotation=45, ha="right")

plt.tight_layout()
plt.savefig("visualize_variants/plots/oncogenes_in_hotspots.png", dpi=300)
plt.show() 

print("Plotting complete! Plot saved as 'visualize_variants/plots/oncogenes_in_hotspots.png'\n")

# ------------------------------------------------------------
# Gene-level Hotspot Fraction 
# ------------------------------------------------------------

print("------------------------------------------------------")
print("GENE-LEVEL HOTSPOT FRACTION") 
print("------------------------------------------------------\n")

# Find out whether the genes are hotspot-driven or just frequently mutated 

print("Finding the fraction of oncogenic variants in cancer hotspots for highly mutated genes...\n")

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
print(oncogenic_gene_summary.head(10),"\n")

# ------------------------------------------------------------
# Plot frequently mutated genes with hotspot enrichment 
# ------------------------------------------------------------

print("------------------------------------------------------")
print("GENE-LEVEL HOTSPOT ENRICHMENT PLOT") 
print("------------------------------------------------------\n")

print("Plotting frequently mutated genes with hotspot enrichment...\n")

top_genes = oncogenic_gene_summary[
  (oncogenic_gene_summary["Total_Oncogenic"] >= 10) &
  (oncogenic_gene_summary["Hotspot_Fraction"] >= 0.3)
].sort_values("Total_Oncogenic", ascending=False) 

print("Example output:")
print(top_genes.head(5),"\n")

plt.figure(figsize=(10, 6)) 
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

plt.xlabel("Gene", fontsize=12)
plt.ylabel("Total Oncogenic Mutations", fontsize=12)
plt.title("Frequently Mutated Genes with Hotspot Enrichment", fontsize=14)
plt.xticks(rotation=45, ha="right")
plt.legend(title="Hotspot Fraction", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("visualize_variants/plots/oncogenes_hotspot_fraction.png", dpi=300)
plt.show()

print("Plotting complete! Figure saved as 'visualize_variants/plots/oncogenes_hotspot_fraction.png'\n")

print("========================================================")
print("Variant Hotspot Analysis Complete!")
print("========================================================\n")
