
# ============================================================
# Top Genes Analysis
# ============================================================

""" 
Script: top_genes.py
Author: Ane Kleiven

This script performs exploratory analyses to identify genes that harbor the
highest numbers of somatic variants across different oncogenicity classes.
Specifically, it quantifies and compares the distribution of Oncogenic,
Likely Oncogenic, and Likely Neutral variants across genes, with a focus on
genes enriched for oncogenic variants.

The analysis aims to highlight recurrently altered cancer genes and to
examine whether genes with many oncogenic variants also accumulate a
substantial fraction of likely neutral variation.

Script content:
--------------
1. Load annotated somatic variant data
2. Extract variants by oncogenicity class (Oncogenic, Likely Oncogenic,
   Likely Neutral)
3. Identify and visualize the top genes by number of variants per class
4. Compare oncogenic and likely neutral variant counts within the genes with the most oncogenic variants 
5. Visualize both absolute counts and relative (percentage) distributions
   of oncogenicity classes per gene

All plots are saved in:
    plots/top_genes
"""

print("-"*50)
print("Top genes analysis🤓")
print("-"*50)

#--------------------------------------------------------------------
# Import libraries
#--------------------------------------------------------------------

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os 
import argparse
from pathlib import Path

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Explore top genes within oncogenicity classes"
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

save_dir = "plots/top_genes"
os.makedirs(save_dir, exist_ok=True) 

#--------------------------------------------------------------------
# Load variant data
#--------------------------------------------------------------------

print(f"Loading variant file:\n{args.variants}")

variants = pd.read_csv(args.variants, sep="\t", low_memory=False)

print(f"Loaded {len(variants)} variants.")
print("-"*30)

#--------------------------------------------------------------------
# Check number of unique genes
#--------------------------------------------------------------------

# Check number of unique genes 
n_unique_genes = variants["Hugo_Symbol"].nunique() 

print(f"The variant dataset consists of {n_unique_genes:,} unique genes.")

#--------------------------------------------------------------------
# Extract 'Oncogenic' variants
#--------------------------------------------------------------------

print("Extracting oncogenic variants..")

oncogenic_variants = variants[variants['ONCOGENIC'] == 'Oncogenic']

oncogenic_genes = ( 
  oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

n_unique_onco = oncogenic_variants["Hugo_Symbol"].nunique() 
print(f"Total number of genes with oncogenic variants: {n_unique_onco}")

#--------------------------------------------------------------------
# Extract 'Likely Oncogenic' variants
#--------------------------------------------------------------------

print("Extracting likely oncogenic variants..")

likely_oncogenic_variants = variants[variants['ONCOGENIC'] == 'Likely Oncogenic']

likely_oncogenic_genes = ( 
  likely_oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

n_unique_likely_onco = likely_oncogenic_variants["Hugo_Symbol"].nunique() 
print(f"Total number of genes with likely oncogenic variants: {n_unique_likely_onco}")

#--------------------------------------------------------------------
# Extract 'Likely Neutral' variants
#--------------------------------------------------------------------

print("Extracting likely neutral variants..")

neutral_variants = variants[variants['ONCOGENIC'] == 'Likely Neutral']

neutral_genes = ( 
  neutral_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
  )

n_unique_neutral = neutral_variants["Hugo_Symbol"].nunique() 
print(f"Total number of genes with neutral variants: {n_unique_neutral}")

#--------------------------------------------------------------------
# Check number of overlapping genes
#--------------------------------------------------------------------

from matplotlib_venn import venn2

oncogenic_genes_set = set(oncogenic_variants["Hugo_Symbol"])
neutral_genes_set = set(neutral_variants["Hugo_Symbol"])

overlap = oncogenic_genes_set & neutral_genes_set

print(f"Number of overlapping genes: {len(overlap)}")
print(f"Overlapping genes: {overlap}")

venn2([oncogenic_genes_set, neutral_genes_set], 
      set_labels=('Oncogenic', 'Neutral'),
      set_colors=('#c4314a', '#88aed1'),
      alpha=0.6)

plt.savefig(f"{save_dir}/venn.png", dpi=300, bbox_inches="tight")

#-------------------------------------------------------------------
# Visualize top 30 genes for all oncogenicity classes 
#--------------------------------------------------------------------

print("-"*30)
print("Plot top genes for all \noncogenicity class")
print("-"*30)

sns.set_theme(style="white") 

# Function to make consistent plots for each oncogenicity class

def plot_top_genes(df, title, color, plotname):
    """
    Create barplot for top genes based on oncogenicity class
    """

    plt.figure(figsize=(8,5))
    
    # sort count values in descending order 
    df_sorted = df.sort_values("Variant_Count", ascending=False)

    # create barplot
    sns.barplot(
        data=df_sorted,
        x="Gene",
        y="Variant_Count",
        color=color,
        edgecolor="black", 
        linewidth=0.8
    )

    sns.despine() 

    # Style titles and labels
    plt.title(title, fontsize=14, pad=15)
    plt.xlabel("Gene", fontsize=12, labelpad=10)
    plt.ylabel("Number of Variants", fontsize=12, labelpad=10)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    plt.tight_layout()
    plt.savefig(f"{save_dir}/{plotname}.png", dpi=300, bbox_inches="tight")
    plt.show()
     
# ============================================================
# Plot each oncogenicity class 
# ============================================================

plot_top_genes(
    oncogenic_genes.head(30),
    "Top Genes (Oncogenic Variants)",
    color="#c4314a",
    plotname="top_oncogenic" 
)

plot_top_genes(
    likely_oncogenic_genes.head(30),
    "Top Genes (Likely Oncogenic Variants)",
    color="#FF834E",
    plotname="top_likely_oncogenic" 
)

plot_top_genes(
    neutral_genes.head(30),
    "Top Genes (Likely Neutral Variants)",
    color="#88aed1",
    plotname="top_likely_neutral"
)

print(f"Plotting complete! Plot saved in folder '{save_dir}'\n")


# ============================================================
# Class distribution in the top oncogenic genes 
# (Oncogenic and Likely Neutral)
# ============================================================

print("-"*30)
print("Count class distribution \nin top oncogenic genes")
print("-"*30)

print("Exploring oncogenicity distribution within the gene with the highest number of oncogenic variants..")

top_onco_genes = oncogenic_genes.head(20)["Gene"] 

top_genes_variants = variants[variants["Hugo_Symbol"].isin(top_onco_genes)]

distribution = ( 
    top_genes_variants
    .groupby(["Hugo_Symbol", "ONCOGENIC"])
    .size() 
    .reset_index(name="Count")
)

print("Example output oncogenicity distribution:")
print(distribution.head())

# extract only oncogenic and likely netrual variants
wanted_classes = ["Oncogenic", "Likely Neutral"]

distribution_filtered = distribution[
    distribution["ONCOGENIC"].isin(wanted_classes)]

plt.figure(figsize=(10,6)) 

sns.set_style(style="white") 

sns.barplot(
    data=distribution_filtered,
    x="Hugo_Symbol",
    y="Count", 
    hue="ONCOGENIC", 
    palette=["#88aed1","#c4314a"],
    edgecolor="black", 
    linewidth=0.5
)

sns.despine() 

print("\nPlotting distribution..")

plt.title("Class Distribution in Genes with Most Oncogenic Variants", fontsize=14, pad=15)
plt.xlabel("Gene", fontsize=12, labelpad=10)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10) 
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity Class", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, title_fontsize=11, markerscale=0.9, frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/distribution_top_onco.png", dpi=300, bbox_inches="tight")
plt.show() 

print(f"Plotting complete! Plot saved in folder '{save_dir}'\n")

# ============================================================
# Class distribution in genes with most oncogenic variants (PIVOT)
# (Oncogenic and Likely Neutral) 
# ============================================================

print("-"*30)
print("Percentage class distribution \nin genes with most oncogenic variants")
print("-"*30)

print("Plotting the percentage class distribution in genes with most oncogenic variants..")

pivot = distribution_filtered.pivot(
    index="Hugo_Symbol",
    columns="ONCOGENIC",
    values="Count"
).fillna(0)

pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100 
pivot_pct = pivot_pct.sort_values(by="Oncogenic", ascending=False)

sns.set_style(style="white") 

pivot_pct.plot(
    kind="bar",
    stacked=True, 
    figsize=(10,6),
    color=["#88aed1","#c4314a"],
    edgecolor="black",
    linewidth=0.6
)

sns.despine() 

plt.title("Oncogenicity Distribution in Genes with Most Oncogenic Variants (%)", fontsize=14, pad=15) 
plt.xlabel("Gene", fontsize=12, labelpad=10)
plt.ylabel("Percentage of Variants (%)", fontsize=12, labelpad=10) 
plt.xticks(rotation=45, ha="right", fontsize=9) 
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity Class", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, title_fontsize=11, markerscale=0.9, frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/percentage_top_onco.png", dpi=300, bbox_inches="tight")
plt.show() 

print(f"Plotting complete! Plot saved in folder '{save_dir}'\n")

