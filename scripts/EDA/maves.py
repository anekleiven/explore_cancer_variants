
# ====================================================================
# Variants MAVEs Analysis Script 
# ====================================================================

"""
Script: maves.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants with different oncogenicity distribute across 
Multiplexed Assays of Variant Effect (MAVEs)

Major outputs:
--------------
    1. Number of variants with MAVE scores
    2. Number of variants with MAVE scores per oncogenicity class (oncogenic vs neutral)
    3. Top genes with MAVE scores
    4. Filter out wanted MAVE experiments (by number of variants with MAVE score)
    5. Box plot of MaveDB score distributions within top MAVE experiments (single genes)
    6. KDE plot of MaveDB score distributions within top MAVE experiments (single genes)

All plots are saved in:
   plots/maves

"""

print("-"*50)
print("MAVE score analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse
from pathlib import Path

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Explore MaveDB scores in variant data."
    ) 

    parser.add_argument(
        "--variants", 
        type=Path, 
        required=True, 
        help="Path to the input file with variant data."
    )

    parser.add_argument(
        "--output", 
        type=Path, 
        required=True, 
        help="Path to the filtered MAVE output file."
    )

    return parser.parse_args() 


args = getargs() 

# ------------------------------------------------------------
# Create output directory 
# ------------------------------------------------------------

save_dir = "plots/maves"
os.makedirs(save_dir, exist_ok=True) 

# ------------------------------------------------------------
# Load variant data 
# ------------------------------------------------------------

print(f"Loading variant file:\n{args.variants}")

variants = pd.read_csv(
    args.variants,
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} somatic variants (exploded version).")
print("-"*30)

# ------------------------------------------------------------
# Explore number of variants with MAVE scores
# ------------------------------------------------------------

# remove rows with no Mave score 
variants_with_mave = variants[variants["MaveDB_score"].notna()]

# make sure mave column is numeric 
variants_with_mave["MaveDB_score"] = pd.to_numeric(variants_with_mave["MaveDB_score"], errors='coerce')

# display number of mave rows 
print(f"Number of rows with MAVE scores: {len(variants_with_mave):,} rows.") 
print(f"Percentage of rows with MAVE scores: {(len(variants_with_mave)/len(variants)*100):.2f}%.\n")

# find the number of unique Mave-variants
variant_id_cols = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
wanted_classes = ["Oncogenic", "Likely Neutral"]

mave_subset = variants_with_mave[variants_with_mave["ONCOGENIC"].isin(wanted_classes)]
unique_mave_var = mave_subset.drop_duplicates(variant_id_cols) 

unique_counts = unique_mave_var["ONCOGENIC"].value_counts() 

# print number of unique variants per class 
print("Unique variants with MAVE scores:")
for c in wanted_classes:
    count = unique_counts.get(c,0)
    print(f"{c}: {count}")
print("-"*30)

# ------------------------------------------------------------
# Explore number of variants with MAVE scores for
# the two oncogenicity classes 
# ------------------------------------------------------------

print("Plotting MAVE variants per oncogenicity class...")

unique_counts_df = (
    unique_mave_var["ONCOGENIC"]
    .value_counts()
    .reset_index(name="count")
    .rename(columns={"index": "ONCOGENIC"})
)

oncogenicity_palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1"
}

plt.figure(figsize=(8,5))
sns.set_style("white")

sns.barplot(
    data=unique_counts_df, 
    x="ONCOGENIC",
    y="count",
    hue="ONCOGENIC", 
    palette=oncogenicity_palette, 
    edgecolor="black",
    linewidth=0.5) 

sns.despine() 

plt.title("Variants with MAVE Scores per Oncogenicity Class", fontsize=14, pad=15)
plt.xlabel("Oncogenicity Class", fontsize=12, labelpad=10)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10) 
plt.xticks(rotation=0, ha='right', fontsize=9)

plt.tight_layout() 
plt.savefig(f"{save_dir}/oncogenicity_classes_maves.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/oncogenicity_classes_maves.png'\n")
print("-"*30)

# ------------------------------------------------------------
# Explore genes with MaveDB_score
# ------------------------------------------------------------

# Summary per gene 
gene_summary_unique = (
  unique_mave_var.groupby("Hugo_Symbol")
  .size()
  .reset_index(name="Variant_Count")
  .rename(columns={"Hugo_Symbol": "Gene"})
  .sort_values("Variant_Count", ascending=False)
)
# Print summary 
print("Number of MaveDB_scores per gene:")
print(gene_summary_unique.head(10), "\n") 

# Plot summary 
print("Plotting genes with MaveDB score (by variant count)..")

sns.set_style("white")
plt.figure(figsize=(8,5))

sns.barplot(data=gene_summary_unique.head(10), 
            x="Gene",
            y="Variant_Count",
            color="#c4314a",
            edgecolor="black",
            linewidth=0.5) 

sns.despine()

plt.title("Genes with MAVE Scores", fontsize=14, pad=15)
plt.xlabel("Gene", fontsize=12, labelpad=10)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10) 

plt.xticks(rotation=0, ha='right', fontsize=9)
plt.tight_layout() 
plt.savefig(f"{save_dir}/top_genes_with_maves.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/top_genes_with_maves.png'")
print("-"*30)

# ------------------------------------------------------------
# Find oncogenicity distribution per gene and MaveDB score set: 
# ------------------------------------------------------------

variants_with_mave_filtered = variants_with_mave[variants_with_mave["ONCOGENIC"].isin(wanted_classes)]

mave_summary = (
  variants_with_mave_filtered.groupby(["Hugo_Symbol", "MaveDB_urn", "ONCOGENIC"])
  .size()) 

print("Oncogenic and likely neutral variant counts per gene and MaveDB score set:")
print(mave_summary, "\n")

# -------------------------------------------------------------
# Explore MAVE genes 
# -------------------------------------------------------------

print("Filtering out wanted MaveDB score sets..\n")

# Wanted MAVE score sets
wanted_urns = ["urn:mavedb:00000068-a-1", "urn:mavedb:00000013-a-1", "urn:mavedb:00000115-a-7", "urn:mavedb:00000081-a-2"] 
mave_filtered = variants_with_mave_filtered[variants_with_mave_filtered["MaveDB_urn"].isin(wanted_urns)]

# -------------------------------------------------------------
# Plot score distributions for wanted MaveDB score sets
# -------------------------------------------------------------

print("Plotting score distributions for wanted MaveDB score sets..")

palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

#--------------
# KDE plot
#--------------

unique_genes = mave_filtered["Hugo_Symbol"].unique() 

for gene in unique_genes: 
    gene_data = mave_filtered[mave_filtered["Hugo_Symbol"] == gene]

    plt.figure(figsize=(8,5))

    sns.kdeplot(
        data=gene_data,
        x="MaveDB_score",
        hue="ONCOGENIC",
        hue_order=wanted_classes,
        palette=palette,
        fill=True,
        alpha=0.5
    )

    n_neutral = (gene_data["ONCOGENIC"] == "Likely Neutral").sum()
    n_oncogenic = (gene_data["ONCOGENIC"] == "Oncogenic").sum()

    
    plt.gca().annotate(
        f"n(Likely Neutral) = {n_neutral}\nn(Oncogenic) = {n_oncogenic}",
        xy=(0.05, 0.95),
        xycoords="axes fraction",
        fontsize=8,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.5)
    )
    
    sns.despine()

    plt.title(f"MAVE Score Distribution: {gene}", fontsize=14, pad=15)
    plt.xlabel("MAVE Functional Score (MaveDB)", fontsize=12, labelpad=10)
    plt.ylabel("Density", fontsize=12, labelpad=10)
    
    plt.tight_layout()
    file_name = f"mave_density_{gene}.png"
    plt.savefig(f"{save_dir}/{file_name}", dpi=300, bbox_inches="tight")
    
    print(f"Saved plot for {gene} to {file_name}")

#--------------
# Box plot
#--------------

unique_genes = mave_filtered["Hugo_Symbol"].unique() 

for gene in unique_genes: 
    gene_data = mave_filtered[mave_filtered["Hugo_Symbol"] == gene]

    plt.figure(figsize=(8,5))
    sns.set_style("white")

    sns.boxplot(
        data=gene_data,
        x="ONCOGENIC",
        y="MaveDB_score",
        order=wanted_classes,
        hue="ONCOGENIC",
        hue_order=wanted_classes,
        palette=palette,
        legend=False, 
        width=0.6,
        linewidth=1.2
    )

    sns.stripplot(
        data=gene_data,
        x="ONCOGENIC",
        y="MaveDB_score",
        order=wanted_classes,
        color="black",
        alpha=0.3,
        size=4,
        jitter=True
    )

    n_neutral = (gene_data["ONCOGENIC"] == "Likely Neutral").sum()
    n_oncogenic = (gene_data["ONCOGENIC"] == "Oncogenic").sum()

    
    plt.gca().annotate(
        f"n(Likely Neutral) = {n_neutral}\nn(Oncogenic) = {n_oncogenic}",
        xy=(0.05, 0.95),
        xycoords="axes fraction",
        fontsize=9,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.5)
    )
    
    sns.despine()

    plt.title(f"MAVE Score Distribution: {gene}", fontsize=14, pad=15)
    plt.xlabel("MAVE Functional Score (MaveDB)", fontsize=12, labelpad=10)
    plt.ylabel("Density", fontsize=12, labelpad=10)
    
    plt.tight_layout()

    file_name = f"mave_boxplot_{gene}.png"
    plt.savefig(f"{save_dir}/{file_name}", dpi=300, bbox_inches="tight")
    
    print(f"Saved plot for {gene} to {file_name}")


#----------------------------------------------------------
# Save filtered MAVE data as .tsv for statistics
#----------------------------------------------------------

mave_filtered.to_csv(args.output, sep="\t", index=False)

print("-"*50)
print(f"Saved filtered MAVE data to: \n {args.output}")
print("-"*50 + "\n")

print("MAVE exploratory analyses complete!🥳🥳\n")