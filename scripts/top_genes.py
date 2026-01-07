"""
# ============================================================
# ANALYSIS: Top Genes by Variant Oncogenicity
# ============================================================

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
4. Compare oncogenic and likely neutral variant counts within the top
   oncogenic genes
5. Visualize both absolute counts and relative (percentage) distributions
   of oncogenicity classes per gene

All plots are saved in:
    explore_cancer_variants/plots/
"""

print("\n========================================================")
print("TOP GENES ANALYSIS")
print("========================================================")

#--------------------------------------------------------------------
# Import libraries
#--------------------------------------------------------------------

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os 

# directory to save plots 
save_dir = "explore_cancer_variants/plots"
os.makedirs(save_dir, exist_ok=True) 

#--------------------------------------------------------------------
# Load variant data
#--------------------------------------------------------------------

print("\n------------------------------------------------------")
print("LOAD VARIANT DATA")
print("------------------------------------------------------\n")

print("Loading variant data...\n")

variants = pd.read_csv(
  "annotation_pipeline/output/variants_with_func_sites.tsv", 
  sep="\t", 
  low_memory=False
  )

#--------------------------------------------------------------------
# Extract 'Oncogenic' variants
#--------------------------------------------------------------------

print("\n------------------------------------------------------")
print("EXTRACT ONCOGENIC VARIANTS")
print("------------------------------------------------------\n")

print("Extracting oncogenic variants...\n")

oncogenic_variants = variants[variants['ONCOGENIC'] == 'Oncogenic']

oncogenic_genes = ( 
  oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

#--------------------------------------------------------------------
# Extract 'Likely Oncogenic' variants
#--------------------------------------------------------------------

print("\n------------------------------------------------------")
print("EXTRACT LIKELY ONCOGENIC VARIANTS")
print("------------------------------------------------------\n")

print("Extracting likely oncogenic variants...\n")

likely_oncogenic_variants = variants[variants['ONCOGENIC'] == 'Likely Oncogenic']

likely_oncogenic_genes = ( 
  likely_oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

#--------------------------------------------------------------------
# Extract 'Likely Neutral' variants
#--------------------------------------------------------------------

print("\n------------------------------------------------------")
print("EXTRACT LIKELY NEUTRAL VARIANTS")
print("------------------------------------------------------\n")

print("Extracting likely neutral variants...\n")

neutral_variants = variants[variants['ONCOGENIC'] == 'Likely Neutral']

neutral_genes = ( 
  neutral_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

#--------------------------------------------------------------------
# Visualize top 30 genes for all oncogenicity classes 
#--------------------------------------------------------------------

print("\n------------------------------------------------------")
print("VISUALIZATION OF TOP GENES PER ONCOGENICITY CLASS")
print("------------------------------------------------------\n")

sns.set_theme(style="whitegrid", context="talk") 

# Function to make consistent plots for each oncogenicity class

def plot_top_genes(df, title, color, plotname):
    """
    Create a consistent barplot for top genes based on oncogenicity class
    """
    plt.figure(figsize=(8,5))
    
    # sort count values in descending order 
    df_sorted = df.sort_values("Variant_Count", ascending=False)

    # create barplot
    sns.barplot(
        data=df_sorted,
        x="Gene",
        y="Variant_Count",
        color=color
    )

    # Add value labels above bars
    for i, v in enumerate(df_sorted["Variant_Count"]):
        plt.text(i, v + 0.5, str(v), ha="center", va="bottom", fontsize=7, color="black")

    # Style titles and labels
    plt.title(title, fontsize=14, pad=10)
    plt.xlabel("Gene", fontsize=12)
    plt.ylabel("Number of Variants", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    # Clean style
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{save_dir}/{plotname}.png", dpi=300, bbox_inches="tight")
    plt.show()
     

# ============================================================
# Plot each oncogenicity class 
# ============================================================

print("Plotting top genes per oncogenicity class...\n")

plot_top_genes(
    oncogenic_genes.head(30),
    "Top Genes by Number of Oncogenic Variants",
    color="#C4473B",
    plotname="top_oncogenic" 
)

plot_top_genes(
    likely_oncogenic_genes.head(30),
    "Top Genes by Number of Likely Oncogenic Variants",
    color="#D98C6A",
    plotname="top_likely_oncogenic" 
)

plot_top_genes(
    neutral_genes.head(30),
    "Top Genes by Number of Likely Neutral Variants",
    color="#7e8aa2",
    plotname="top_likely_neutral"
)

print("Plotting complete. All plots saved in 'explore_cancer_variants/plots'\n")

# ============================================================
# Class distribution in the top oncogenic genes 
# (Oncogenic and Likely Neutral)
# ============================================================

print("\n------------------------------------------------------")
print("COUNTS CLASS DISTRIBUTION IN THE TOP ONCOGENIC GENES")
print("------------------------------------------------------\n")

print("Exploring oncogenicity distribution within the top oncogenic genes...\n")

top_onco_genes = oncogenic_genes.head(30)["Gene"] 

top_genes_variants = variants[variants["Hugo_Symbol"].isin(top_onco_genes)]

distribution = ( 
    top_genes_variants
    .groupby(["Hugo_Symbol", "ONCOGENIC"])
    .size() 
    .reset_index(name="Count")
)

print("Example output oncogenicity distribution:\n")
print(distribution.head())

# extract only oncogenic and likely netrual variants
wanted_classes = ["Oncogenic", "Likely Neutral"]

distribution_filtered = distribution[
    distribution["ONCOGENIC"].isin(wanted_classes)]

sns.set_style(style="whitegrid") 

plt.figure(figsize=(10,6)) 

sns.barplot(
    data=distribution_filtered,
    x="Hugo_Symbol",
    y="Count", 
    hue="ONCOGENIC", 
    palette=["#7e8aa2","#C4473B"]
)

print("\nPlotting distribution...\n")

plt.title("Distribution of Variant Oncogenicity in Top Oncogenic Genes", fontsize=14, pad=10)
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12) 
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(
    title="Oncogenicity Class",
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    fontsize=10,
    title_fontsize=11,
    markerscale=0.9
)

plt.tight_layout()

plt.savefig(f"{save_dir}/distribution_top_onco.png", dpi=300, bbox_inches="tight")

plt.show() 

print("Plotting complete! Plot saved in folder 'explore_cancer_variants/plots'\n")

# ============================================================
# Class distribution in the oncogenic genes (PIVOT)
# (Oncogenic and Likely Neutral) 
# ============================================================

print("\n------------------------------------------------------")
print("PERCENTAGE CLASS DISTRIBUTION IN THE TOP ONCOGENIC GENES")
print("------------------------------------------------------\n")

print("Plotting the percentage class distribution in the top oncogenic genes...\n")

pivot = distribution_filtered.pivot(
    index="Hugo_Symbol",
    columns="ONCOGENIC",
    values="Count"
).fillna(0)

pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100 

sns.set_style(style="whitegrid") 

pivot_pct.plot(
    kind="bar",
    stacked=True, 
    figsize=(10,6),
    color=["#7e8aa2","#C4473B"],
    edgecolor="0.1",
    linewidth=0.3
)

plt.title("Percentage Distribution of Variant Oncogenicity in Top Oncogenic Genes", fontsize=14, pad=10) 
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Percentage of Variants (%)", fontsize=12) 
plt.xticks(rotation=45, ha="right", fontsize=9) 
plt.yticks(fontsize=9)
plt.legend(
    title="Oncogenicity Class",
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    fontsize=10,
    title_fontsize=11,
    markerscale=0.9
)

plt.tight_layout()

plt.savefig(f"{save_dir}/percentage_top_onco.png", dpi=300, bbox_inches="tight")

plt.show() 

print("Plotting complete! Plot saved in folder 'explore_cancer_variants/plots'\n")


print("========================================================")
print("Top genes analysis complete!\n")
