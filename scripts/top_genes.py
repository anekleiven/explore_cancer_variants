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
6. Perform statistics on all features for the top 10 oncogenic genes 

All plots are saved in:
    plots/
"""

print("\nTop genes analysis🤓")
print("-"*50)

#--------------------------------------------------------------------
# Import libraries
#--------------------------------------------------------------------

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os 

# directory to save plots 
save_dir = "plots"
os.makedirs(save_dir, exist_ok=True) 

#--------------------------------------------------------------------
# Load variant data
#--------------------------------------------------------------------

print("Loading variant data..")

variants = pd.read_csv(
  "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv", 
  sep="\t", 
  low_memory=False
  )

print(f"Loaded {len(variants)} variants.")
print("-"*30)

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

#--------------------------------------------------------------------
# Visualize top 30 genes for all oncogenicity classes 
#--------------------------------------------------------------------

print("-"*30)
print("Plot top genes for all oncogenicity class")
print("-"*30)

sns.set_theme(style="whitegrid", context="talk") 

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

print("Plotting complete. All plots saved in 'plots'\n")

# ============================================================
# Class distribution in the top oncogenic genes 
# (Oncogenic and Likely Neutral)
# ============================================================

print("-"*30)
print("Count class distribution in top oncogenic genes")
print("-"*30)

print("Exploring oncogenicity distribution within the top oncogenic genes..")

top_onco_genes = oncogenic_genes.head(30)["Gene"] 

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

sns.set_style(style="whitegrid") 

plt.figure(figsize=(10,6)) 

sns.barplot(
    data=distribution_filtered,
    x="Hugo_Symbol",
    y="Count", 
    hue="ONCOGENIC", 
    palette=["#7e8aa2","#C4473B"]
)

print("\nPlotting distribution..")

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

print("Plotting complete! Plot saved in folder 'plots'\n")

# ============================================================
# Class distribution in the oncogenic genes (PIVOT)
# (Oncogenic and Likely Neutral) 
# ============================================================

print("-"*30)
print("Percentage class distribution in top oncogenic genes")
print("-"*30)

print("Plotting the percentage class distribution in the top oncogenic genes..")

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

print("Plotting complete! Plot saved in folder 'plots'\n")


# ============================================================
# Statistics: Top genes
# ============================================================

print("-"*30)
print("Perform statistics on top 10 oncogenic genes")
print("-"*30)

# import libraries
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio 
from statsmodels.stats.multitest import multipletests
import numpy as np

# extract top 10 oncogenic genes
top_10_onco = oncogenic_genes.head(10)["Gene"].tolist()  

# extract variants in top 10 oncogenic genes 
top_10_variants = variants[variants["Hugo_Symbol"].isin(top_10_onco)]

# add feature 
top_10_variants["has_gnomAD_AF"] = (
    (top_10_variants["gnomAD_AF"].notna()) &
    (top_10_variants["gnomAD_AF"] != "NA") & 
    (top_10_variants["gnomAD_AF"] != "")
)

variants["has_gnomAD_AF"] = ( 
    (variants["gnomAD_AF"].notna()) & 
    (variants["gnomAD_AF"] != "NA") & 
    (variants["gnomAD_AF"] != "") 
)

# define the features 
features = ["gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE", "Germline_Proximity", "MaveDB_score"]

def analyze_top_genes(df, label="Dataset"):
    print("\n" + "-"*50)
    print(f"Statistics: {label}")
    print("-"*50 + "\n")

    results = []

    for f in features:

        if f in ["gnomAD_AF", "Germline_Proximity", "MaveDB_score"]: 
            oncogenic = df[df["ONCOGENIC"] == "Oncogenic"][f].dropna()
            neutral = df[df["ONCOGENIC"] == "Likely Neutral"][f].dropna()

            if len(oncogenic) == 0 or len(neutral) == 0:
                print(f"[{f}] Skipped (not enough data)\n")
                continue

            # perform Mann-Whitney U test 
            stat, p = mannwhitneyu(oncogenic, neutral, alternative="two-sided")
            # calculate rank-biserial correlation 
            # (effect size for mann-whitney u) 
            n1 = len(oncogenic)
            n2 = len(neutral)
            r = (2 * stat) / (n1 * n2) - 1

            # calculate probability 
            probability = (1+r)/2 

            print(f"[{f}]")
            print(f"Mann-Whitney U: {stat:.3f}, p-value: {p:.4f}")
            print(f"{'Reject H₀: distributions differ.' if p < 0.05 else 'Failed to reject H₀.'}")
            print(f"Rank-biserial r: {r:.3f} | P(oncogenic > neutral): {probability*100:.2f}%\n")
            results.append({"feature": f, "test": "Mann-Whitney", "p_value": p})


        else: 
            onc_in  = len(df[(df["ONCOGENIC"] == "Oncogenic") & (df[f] == True)])
            onc_out = len(df[(df["ONCOGENIC"] == "Oncogenic") & (df[f] == False)])
            neu_in  = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (df[f] == True)])
            neu_out = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (df[f] == False)])

            observed_table = [[onc_in, neu_in], [onc_out, neu_out]]

            # Odds ratio and 95% CI
            or_result = scipy_odds_ratio(observed_table)
            ci = or_result.confidence_interval(confidence_level=0.95)

            # select test based on number of observations in each cell
            # < 5 variants in one cell: Fisher, else: Chi-Square
            if min(onc_in, onc_out, neu_in, neu_out) < 5: 
                _, p = fisher_exact(observed_table)           
                test_used = "Fisher"
            else: 
                chi2, p, _, _ = chi2_contingency(observed_table)     
                test_used = "Chi-Square"
                n = sum(sum(row) for row in observed_table)
                k = 1  # only for 2x2 table
                cramers_v = np.sqrt(chi2 / (n * k))

            if test_used == "Fisher":
                print(f"[{f}]")
                print(f"Test: {test_used}")
                print(f"OR: {or_result.statistic:.3f} (95% CI: {ci.low:.3f}–{ci.high:.3f}) | p-value: {p:.4f}\n")
            else:
                print(f"[{f}]")
                print(f"Test: {test_used}")
                print(f"OR: {or_result.statistic:.3f} (95% CI: {ci.low:.3f}–{ci.high:.3f})")
                print(f"p-value: {p:.4f} | Cramer's V: {cramers_v:.3f}\n")

            results.append({"feature": f, "test": test_used, "p_value": p})
                

    results_df = pd.DataFrame(results)
    _, q_values, _, _ = multipletests(results_df["p_value"], method="fdr_bh")
    results_df["q_value"] = q_values.round(4)
    results_df["p_value"] = results_df["p_value"].round(4)


    print(f"\n{'-'*60}")
    print("FDR-corrected results (Benjamini-Hochberg)")
    print(f"{'-'*60}")
    print(results_df.to_string(index=False))


analyze_top_genes(top_10_variants, label="Top 10 Oncogenic genes")
analyze_top_genes(variants, label="All variants") 





