"""
====================================================================
Variants MAVEs Analysis Script 
====================================================================

Script: variants_with_maves.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants with different oncogenicity distribute across 
Multiplexed Assays of Variant Effect (MAVEs)

Major outputs:
--------------
1. Number of variants with MAVEdb scores
2. Number of MaveDB_scores for each oncogenicity class 
3. Top genes with MaveDB_score 
4. Oncogenicity distribution within top genes with MaveDB_score 
5. Oncogenic vs neutral counts within top genes with MaveDB_score
6. Descriptive statistics of MAVE scores 
7. Boxplot of MAVE score distributions per class 
8. Density plot of MAVE score distributions per class 


All plots are saved in:
   plots/

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

# ------------------------------------------------------------
# Load variant data 
# ------------------------------------------------------------

print("\nLoading variant data..")

variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} somatic variants.")
print("-"*30)

# ------------------------------------------------------------
# Explore number of variants with MAVE scores
# ------------------------------------------------------------

variants_with_mave = variants[variants["MaveDB_score"].notna()]

print(f"Number of variants with MAVE scores: {len(variants_with_mave):,} variants.") 
print(f"Percentage of variants with MAVE scores: {(len(variants_with_mave)/len(variants)*100):.2f}%.")
print("-"*30)

# ------------------------------------------------------------
# Explore number of variants with MAVE scores for
# different oncogenicity classes 
# ------------------------------------------------------------

# create summary table 
oncogenicity_summary = variants_with_mave.groupby("ONCOGENIC").size().reset_index(name="Variant_Count")

# rename column
oncogenicity_summary.columns = ["Oncogenicity_Class", "Variant_Count"]

# sort by count 
oncogenicity_summary = oncogenicity_summary.sort_values("Variant_Count", ascending = False) 

print("Number of MAVE scores per oncogenicity class:")
print(oncogenicity_summary,"\n") 

# plot summary 
print("Plotting MAVE score per oncogenicity class...")

oncogenicity_palette = {
    "Oncogenic": "#c4314a",
    "Likely Oncogenic": "#D98C6A",
    "Likely Neutral": "#88aed1",
    "Inconclusive": "#f9c74f",
    "Unknown": "#848a8e",
    "Resistance": "#ba7ad4"
}

plt.figure(figsize=(8,5))
sns.barplot(data=oncogenicity_summary, 
            x="Oncogenicity_Class",
            y="Variant_Count",
            hue="Oncogenicity_Class", 
            palette=oncogenicity_palette, 
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Variants with MAVE Scores per Oncogenicity Class", fontsize=14, pad=10)
plt.xlabel("Oncogenicity Class", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12) 
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.tight_layout() 
plt.savefig("plots/oncogenicity_classes_maves.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'plots/oncogenicity_classes_maves.png'\n")
print("-"*30)

# ------------------------------------------------------------
# Explore top genes with MAVEdb_score 
# ------------------------------------------------------------

# summary per gene 
gene_summary = variants_with_mave.groupby("Hugo_Symbol").size().reset_index(name="Variant_Count")

# rename columne
gene_summary.columns = ["Gene", "Variant_Count"]

# sort by count 
gene_summary = gene_summary.sort_values("Variant_Count", ascending = False) 

print("Number of MaveDB_scores per gene:")
print(gene_summary) 

# plot summary 
print("Plotting MaveDB score per gene...")


plt.figure(figsize=(8,5))
sns.barplot(data=gene_summary.head(10), 
            x="Gene",
            y="Variant_Count",
            color="#9bc7de",
            edgecolor="0.1",
            linewidth=0.3) 

plt.title("Top Genes with MAVE Scores", fontsize=14, pad=10)
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12) 
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.tight_layout() 
plt.savefig("plots/top_genes_with_maves.png", dpi=300)
plt.show()

print("Plotting complete! Plot saved as 'plots/top_genes_with_maves.png'")
print("-"*30)

# ------------------------------------------------------------
# Explore oncogenicity within top genes with MAVEdb_score 
# ------------------------------------------------------------

# filter for specific oncogenicity classes
oncogenicity_classes_of_interest = ['Oncogenic','Likely Oncogenic', 'Likely Neutral']
variants_filtered = variants_with_mave[variants_with_mave["ONCOGENIC"].isin(oncogenicity_classes_of_interest)]

# summary per gene and oncogenicity class 
gene_oncogenicity_summary = variants_filtered.groupby(["Hugo_Symbol", "ONCOGENIC"]).size().reset_index(name="Variant_Count")

# sort by count 
gene_oncogenicity_summary = gene_oncogenicity_summary.sort_values("Variant_Count", ascending = False) 

print("Oncogenicity distribution within genes with MaveDB_score:")
print(gene_oncogenicity_summary,"\n") 

# plot stacked bar chart 
# find top genes (total MAVE score count) 
top_genes = gene_oncogenicity_summary.groupby("Hugo_Symbol")["Variant_Count"].sum().nlargest(6).index

# filter to top genes
top_genes_data = gene_oncogenicity_summary[gene_oncogenicity_summary['Hugo_Symbol'].isin(top_genes)]

# pivot for stacked bar plot
top_genes_pivot = top_genes_data.pivot(index="Hugo_Symbol", columns ="ONCOGENIC", values="Variant_Count").fillna(0) 

# sort by total counts 
top_genes_pivot = top_genes_pivot.loc[top_genes_pivot.sum(axis=1).sort_values(ascending=False).index]

three_classes_palette = {
    "Oncogenic": "#c4314a",
    "Likely Oncogenic": "#D98C6A",
    "Likely Neutral": "#88aed1",
}

print("Plotting oncogenicity distribution in top genes with MAVE scores..")

plt.figure(figsize=(8,5)) 
top_genes_pivot.plot(kind='bar', stacked=True, color=[three_classes_palette.get(col, "#cccccc") for col in top_genes_pivot.columns], edgecolor='white', linewidth=0.5, ax=plt.gca()) 

plt.title("Oncogenicity Distribution in Top Genes with MAVE Scores", fontsize=14, pad=10)
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1, 1), loc='upper right')
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.tight_layout() 
plt.savefig("plots/top_genes_mave_oncogenicity.png", dpi=300, bbox_inches='tight')
plt.show() 

print("Stacked bar plot saved as 'plots/top_genes_mave_oncogenicity.png'.")
print("-"*30)

# ------------------------------------------------------------
# Oncogenic vs neutral within top genes with MaveDB_score
# ------------------------------------------------------------

# keep only variants that are oncogenic or likely neutral 
oncogenic_neutral_classes =  ["Oncogenic", "Likely Neutral"]
variants_onco_vs_neutral = variants_with_mave[variants_with_mave["ONCOGENIC"].isin(oncogenic_neutral_classes)]

# create summary
onco_vs_neutral_summary = variants_onco_vs_neutral.groupby(["Hugo_Symbol", "ONCOGENIC"]).size().reset_index(name="Variant_Count") 

# sort by count 
onco_vs_neutral_summary = onco_vs_neutral_summary.sort_values("Variant_Count", ascending = False) 

print("Oncogenic vs neutral counts within genes with MaveDB_score:")
print(onco_vs_neutral_summary,"\n") 

# plot stacked bar chart 
# find top genes (total MAVE score count) 
top_genes_comparison = onco_vs_neutral_summary.groupby("Hugo_Symbol")["Variant_Count"].sum().nlargest(6).index

# filter to top genes
top_genes_comparison_data = onco_vs_neutral_summary[onco_vs_neutral_summary['Hugo_Symbol'].isin(top_genes_comparison)]

# pivot for stacked bar plot
top_genes_comparison_pivot = top_genes_comparison_data.pivot(index="Hugo_Symbol", columns = "ONCOGENIC", values="Variant_Count").fillna(0) 

# sort by total counts 
top_genes_comparison_pivot = top_genes_comparison_pivot.loc[top_genes_comparison_pivot.sum(axis=1).sort_values(ascending=False).index]

two_classes_palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

print("Plotting variant counts in top genes with mave scores..")

plt.figure(figsize=(8,5)) 
top_genes_comparison_pivot.plot(kind='bar', stacked=True, color=[two_classes_palette.get(col, "#cccccc") for col in top_genes_comparison_pivot.columns], edgecolor='white', linewidth=0.5, ax=plt.gca()) 

plt.title("Oncogenic vs Neutral Counts in Top Genes with MAVE Scores", fontsize=14, pad=10)
plt.xlabel("Gene", fontsize=12)
plt.ylabel("Number of Variants", fontsize=12)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1, 1), loc='upper right')
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.tight_layout() 
plt.savefig("plots/onco_neutral_mave_counts.png", dpi=300, bbox_inches='tight')
plt.show() 

print("Stacked bar plot saved as 'plots/onco_neutral_mave_counts.png'.")
print("-"*30)

# ------------------------------------------------------------
# Descriptive statistics of MAVE scores per oncogenicity class 
# ------------------------------------------------------------

maves_score_summary = variants_onco_vs_neutral.groupby("ONCOGENIC")["MaveDB_score"].agg(
    count='count',
    min ='min',
    max='max',
    median='median',
    mean='mean',
    std='std'
).reset_index()

print("Descriptive statistics of MAVE score distributions:")
print(maves_score_summary)
print("-"*30)


# ------------------------------------------------------------
# Box plot og MAVE scores per oncogenicity class 
# ------------------------------------------------------------

# remove rows with missing MAVE data 
mavescore_plot = variants_onco_vs_neutral.dropna(subset=["MaveDB_score"])

# standardize colors 
palette={
  "Oncogenic": "#c4314a",
  "Likely Neutral": "#88aed1"
  }

print("Plotting boxplot of MAVE score data...\n")

plt.figure(figsize=(8,5))

sns.boxplot(
  data=mavescore_plot, 
  x="ONCOGENIC", 
  y="MaveDB_score", 
  palette=palette)

plt.title("MAVE Score by Oncogenicity Class")
plt.xlabel("Oncogenicity")
plt.ylabel("MaveDB Score")

plt.tight_layout() 
plt.savefig("plots/boxplot_maves.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Boxplot saved as 'plots/boxplot_maves.png'.")
print("-"*30)


# ------------------------------------------------------------
# Density plot og MAVE scores per oncogenicity class 
# ------------------------------------------------------------

# Plot settings 
sns.set_theme(style="whitegrid")

# Density plot with oncogenic and neutral variants 
print("Plotting distribution of mave scores per oncogenicity class...\n")

plt.figure(figsize=(8,5)) 
sns.kdeplot(
  data=mavescore_plot,
  x="MaveDB_score",
  hue="ONCOGENIC",
  fill=True,
  common_norm=False,
  palette=palette, 
  alpha=0.5,
  linewidth=0.2, 
  legend=True
)

plt.title("Distribution of MAVE Scores: Oncogenicity", fontsize=14)
plt.xlabel("MaveDB Score (Functional Impact)", fontsize=12)
plt.ylabel("Density", fontsize=12)

plt.tight_layout() 
plt.savefig("plots/mave_scores_comparison.png", bbox_inches='tight')
plt.show()

print("Plotting complete! Plot saved as 'plots/mave_scores_comparison.png'.")
print("-"*30)



# ------------------------------------------------------------
# Density plot og MAVE scores per mutation effect (LOF/GOF) 
# ------------------------------------------------------------

mechanism_map = {
    "Likely Loss-of-function": "Loss-of-function",
    "Loss-of-function": "Loss-of-function",
    "Likely Gain-of-function": "Gain-of-function",
    "Gain-of-function": "Gain-of-function"
}

variants['Simplified_Mechanism'] = variants['MUTATION_EFFECT'].map(mechanism_map)

# Filter data
onco_set = variants[
    (variants["ONCOGENIC"] == "Oncogenic") & 
    (variants["Simplified_Mechanism"])
].copy()

palette = {"Loss-of-function": "#c4314a", "Gain-of-function": "#88aed1"}

# Plot settings 
sns.set_theme(style="whitegrid")
plt.figure(figsize=(8,5)) 

# Density plot with oncogenic and neutral variants 
print("Plotting distribution of mave scores between LoF and GoF variants..")

sns.kdeplot(
  data=onco_set,
  x="MaveDB_score",
  hue="Simplified_Mechanism",
  fill=True,
  common_norm=False,
  palette=palette, 
  alpha=0.5, 
  linewidth=0.2, 
  legend=True
)

plt.gca().get_legend().set_title("Variant Mechanism")
plt.title('Distribution of MAVE Scores: GoF vs. LoF (Oncogenic Variants)', fontsize=14)
plt.xlabel('MaveDB Score (Functional Impact)', fontsize=12)
plt.ylabel('Density', fontsize=12)

plt.tight_layout() 
plt.savefig("plots/mave_scores_lof_gof.png", bbox_inches='tight')
plt.show()

print("Plotting complete! Plot saved as 'plots/mave_scores_lof_gof.png'.")
print("-"*30)


# Filter out variants without MAVE-score 
mave_variants = variants[variants['MaveDB_score'].notna()]

# Make cross table
mave_effect_counts = pd.crosstab(
    mave_variants['MUTATION_EFFECT'], 
    mave_variants['ONCOGENIC'], 
    margins=True, 
    margins_name="Total"
)

print("Variants with MAVE-score per mutation effect:")
print(mave_effect_counts)

print("\nMAVEs analysis complete!🥳🥳\n")