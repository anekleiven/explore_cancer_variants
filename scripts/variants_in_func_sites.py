"""
====================================================================
Variant Functional-Site Analysis Script
====================================================================

Script: variants_in_func_sites.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variants distribute across different functional protein features.

Major outputs:
--------------
1. Counts of Oncogenic vs Likely Neutral variants inside/outside functional sites
2. Expansion of multi functional site annotations (e.g., "Binding site; Region")
3. Variant counts per functional site type per class
4. Fractions of variants in each functional site type
5. Identification of genes enriched in functional sites
6. Comparison of Oncogenic vs Neutral variants in the same genes
7. Oncogenic-to-neutral ratio plots for each functional site type
8. Heatmap to see which genes dominate in each functional site type
9. Statistical analysis: Chi-Square, Cramer's V and odds-ratio 

Key functional site handled:
---------------------
- Binding site
- Modified residue
- Region
- Topological domain

All plots are saved in:
    plots/

"""

print("-"*50)
print("Variant functional site analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# ------------------------------------------------------------
# Setup for analysis 
# ------------------------------------------------------------

# Create output folder if it doesn't exist
os.makedirs("plots", exist_ok=True)

print("Loading variant data..")

variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.")

# ------------------------------------------------------------
# Keep only Oncogenic + Likely Neutral variants
# ------------------------------------------------------------

print("-"*30)
print("Filtering data to only contain oncogenic & likely oncogenic variants..")

classes = ["Oncogenic", "Likely Neutral"]
variants = variants[variants["ONCOGENIC"].isin(classes)]
print(f"Remaining variants after filtering: {len(variants):,}")

# Ensure boolean data type
variants["IN_FUNC_SITE"] = variants["IN_FUNC_SITE"].astype(bool)

# ------------------------------------------------------------
# Summary: variants inside vs outside functional sites
# ------------------------------------------------------------

print("-"*30)
print("Creating summary of variants inside/outside functional sites..")
summary = []
for c in classes:
    subset = variants[variants["ONCOGENIC"] == c]
    total = len(subset)
    inside = subset["IN_FUNC_SITE"].sum()
    fraction_in = inside / total if total > 0 else 0

    summary.append({
        "Class": c,
        "Total_variants": total,
        "Inside_func_site": inside,
        "Outside_func_site": total - inside,
        "Fraction_inside": fraction_in
    })

summary_df = pd.DataFrame(summary)
print("Distribution of variants inside and outside functional sites:")
print(summary_df)

# ------------------------------------------------------------
# Expand FEATURE_TYPE 
# ------------------------------------------------------------

print("-"*30)
print("Expanding FEATURE_TYPE so each functional site type is one row..")

expanded = (
    variants
    .dropna(subset=["FEATURE_TYPE"])
    .assign(FEATURE_TYPE=lambda df: df["FEATURE_TYPE"].str.split(";"))
    .explode("FEATURE_TYPE")
)


expanded["FEATURE_TYPE"] = expanded["FEATURE_TYPE"].str.strip()
print(f"\nExpanded to {len(expanded):,} feature-variant rows.")

# ------------------------------------------------------------
# Count variants per functional site type per class
# ------------------------------------------------------------

print("-"*30)
print("Counting variants per functional site per class..")

counts = (
    expanded
    .groupby(["FEATURE_TYPE", "ONCOGENIC"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Number of variants in each functional site for both classes:")
print(counts)

# ------------------------------------------------------------
# Plot raw counts
# ------------------------------------------------------------

print("-"*30)
print("Plotting raw counts of variants in functional sites..")

palette = {"Oncogenic": "#c4314a", "Likely Neutral": "#88aed1"}

plt.figure(figsize=(8,5))
sns.barplot(
    data=counts,
    x="FEATURE_TYPE",
    y="Variant_Count",
    hue="ONCOGENIC",
    palette=palette, 
    edgecolor="0.1",
    linewidth=0.3
)
plt.title("Number of Variants per Functional Site Type", fontsize=14, pad=10)
plt.xlabel("Functional Site", fontsize=12)
plt.ylabel("Variant Count", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("plots/counts_per_func_site.png", dpi=300)
plt.show()

print(f"Plotting complete. Saved as 'plots/counts_per_func_site.png'")

# ------------------------------------------------------------
# Compute fractions
# ------------------------------------------------------------

print("-"*30)
print("Computing fractions of variants in the different functional sites..")

totals = variants["ONCOGENIC"].value_counts().rename("Total")

counts = counts.merge(totals, left_on="ONCOGENIC", right_index=True)
counts["Fraction"] = counts["Variant_Count"] / counts["Total"]

print("\nFraction of variants in each functional site type per class:")
print(counts)

# ------------------------------------------------------------
# Plot fractions
# ------------------------------------------------------------

print("-"*30)
print("Plotting fractions of variants for each functional site type..\n")

plt.figure(figsize=(8,5))
sns.barplot(
    data=counts,
    x="FEATURE_TYPE",
    y="Fraction",
    hue="ONCOGENIC",
    palette=palette, 
    edgecolor="0.1",
    linewidth=0.3
)
plt.title("Fraction of Variants per Functional Site Type", fontsize=14, pad=10)
plt.xlabel("Functional Site Type", fontsize=12)
plt.ylabel("Fraction", fontsize=12)
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(fontsize=9)
plt.legend(title="Oncogenicity", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("plots/fraction_per_func_site.png", dpi=300)
plt.show()

print("Plotting complete. Saved as 'plots/fraction_per_func_site.png'")

# ------------------------------------------------------------
# Identify genes enriched in functional sites
# ------------------------------------------------------------

print("-"*30)
print("Identifying oncogenic driver genes enriched in functional sites..")

filtered_sites = ["Binding site", "Modified residue", "Region", "Topological domain"]

expanded_filtered = expanded[expanded["FEATURE_TYPE"].isin(filtered_sites)].copy()

# Oncogenic counts
onco = expanded_filtered[expanded_filtered["ONCOGENIC"] == "Oncogenic"].copy()

gene_feature_counts = (
    onco.groupby(["Hugo_Symbol", "FEATURE_TYPE"])
    .size()
    .reset_index(name="Variant_Count")
    .sort_values("Variant_Count", ascending=False)
)

feature_totals = (
    gene_feature_counts.groupby("FEATURE_TYPE")["Variant_Count"]
    .sum()
    .reset_index(name="Feature_Total")
)

gene_feature_fraction = gene_feature_counts.merge(feature_totals, on="FEATURE_TYPE")
gene_feature_fraction["Fraction_of_Feature"] = (
    gene_feature_fraction["Variant_Count"] / gene_feature_fraction["Feature_Total"]
)

print("Example output:")
print(gene_feature_fraction.head(5))

# ------------------------------------------------------------
# Plot top genes per functional site
# ------------------------------------------------------------

print("-"*30)
print("Plotting top genes per functional site..")

for ft in filtered_sites:
    subset = (
        gene_feature_fraction[gene_feature_fraction["FEATURE_TYPE"] == ft]
        .sort_values("Fraction_of_Feature", ascending=False)
        .head(20)
    )
    
    plt.figure(figsize=(8,5))

    plt.bar(subset["Hugo_Symbol"], 
            subset["Fraction_of_Feature"], 
            color="#C4473B", 
            edgecolor="0.1",
            linewidth=0.3)

    plt.title(f"Top Driver Genes in '{ft}'", fontsize=14, pad=10)
    plt.xlabel("Gene", fontsize=12)
    plt.ylabel("Fraction of variants in feature", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    plt.tight_layout()
    plt.savefig(f"plots/topgenes_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved in 'plots'")

# ------------------------------------------------------------
# Compare likely neutral variants to oncogenic variants in top oncogenic driver genes
# ------------------------------------------------------------

print("-"*30)
print("COMPARE VARIANTS IN TOP ONCOGENIC DRIVER GENES")
print("-"*30)

print("Extracting likely neutral variants..")

likely_neutral = expanded_filtered[expanded_filtered["ONCOGENIC"] == "Likely Neutral"].copy()

print("Counting likely neutral variants per functional site type..")
likely_neutral_counts = (
    likely_neutral.groupby(["Hugo_Symbol", "FEATURE_TYPE"])
    .size()
    .reset_index(name="Variant_Count")
)

print("Merging oncogenic variants to likely neutral variants..")
comparison = gene_feature_counts.merge(
    likely_neutral_counts,
    on=["Hugo_Symbol", "FEATURE_TYPE"],
    how="outer",
    suffixes=("_oncogenic", "_likely_neutral")
).fillna(0)

top_onco_genes = gene_feature_counts["Hugo_Symbol"].unique()
comparison_top = comparison[comparison["Hugo_Symbol"].isin(top_onco_genes)].copy()

print("Calculating the oncogenic-neutral ratio for the top driver genes..")
comparison_top["ratio_onco_neutral"] = (
    (comparison_top["Variant_Count_oncogenic"] + 1) /
    (comparison_top["Variant_Count_likely_neutral"] + 1)
)

print("Example output oncogenic-neutral ratio:")
print(comparison_top.head())

# ------------------------------------------------------------
# Plot oncogenic-neutral ratios
# ------------------------------------------------------------

print("-"*30)
print("Plotting oncogenic-neutral ratio for all functional site types..")

for ft in filtered_sites:
    subset = (
        comparison_top[comparison_top["FEATURE_TYPE"] == ft]
        .sort_values("ratio_onco_neutral", ascending=False)
        .head(20)
    )

    plt.figure(figsize=(8,5))
    plt.bar(subset["Hugo_Symbol"], 
            subset["ratio_onco_neutral"], 
            color="#C4473B",
            edgecolor="0.1",
            linewidth=0.3)

    plt.title(f"Oncogenic vs Neutral Variant Ratio in '{ft}'", fontsize=14, pad=10)
    plt.xlabel("Gene", fontsize=12)
    plt.ylabel("Ratio", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=9)
    plt.yticks(fontsize=9)

    plt.tight_layout()
    plt.savefig(f"plots/onco-neutral-ratio_in_{ft}.png", dpi=300)
    plt.show()

print("Plotting complete! Saved in 'plots'")

# ------------------------------------------------------------
# Which genes dominate in each functional site? 
# ------------------------------------------------------------

print("-"*30)
print("Extracting top genes per functional site type..")

top_genes_per_feature = (
    gene_feature_fraction
    .sort_values(['FEATURE_TYPE', 'Fraction_of_Feature'], ascending=[True,False])
    .groupby('FEATURE_TYPE')
    .head()
    .reset_index(drop=True)
)

print("Example output:")
print(top_genes_per_feature.head(10),"\n")

print("Plotting top genes per functional site..")

pivot = top_genes_per_feature.pivot(
    index='Hugo_Symbol', 
    columns='FEATURE_TYPE', 
    values='Fraction_of_Feature'
).fillna(0)

plt.figure(figsize=(10,6))
ax = sns.heatmap(pivot, 
            cmap='Reds', 
            annot=True, 
            fmt='.2f',
            cbar_kws={"label":"Fraction of variants in functional site"},
            )

ax.set_xlabel("Feature Type", fontsize=12)
ax.set_ylabel("Gene (Hugo Symbol)", fontsize=12)

plt.title('Gene Distribution Across Functional Site Types', fontsize=14, pad=10)

plt.tight_layout()
plt.savefig("plots/top_genes_per_functional_site.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/top_genes_per_functional_site.png'")

# ------------------------------------------------------------
# Statistical Test: Chi-Square 
# ------------------------------------------------------------

# Check whether there is an association between oncogenicity and variants
# inside / outside functional sites 

# Model assumptions: 
# Both variables are categorical.
# All observations are independent.
# Each observation must only contribute to one cell.
# The expected frequency in each cell should be at least five.

# Hypotheses: 
# H0:   There is no association between variant classification (oncogenic vs. likely neutral) 
#       and their localization relative to functional sites 
# H1:   There is a statistically significant association between variant classification and 
#       localization relative to functional sites. 

# import libraries
from scipy.stats import chi2_contingency 
import numpy as np 
from scipy.stats.contingency import odds_ratio

print("-"*30)
print("Performing statistics on functional sites data..\n")
print("Running Chi-square test..\n")

# prepare the data 
oncogenic_in = variants.query("ONCOGENIC == 'Oncogenic' and IN_FUNC_SITE")
oncogenic_out = variants.query("ONCOGENIC == 'Oncogenic' and IN_FUNC_SITE == False")
neutral_in = variants.query("ONCOGENIC == 'Likely Neutral' and IN_FUNC_SITE")
neutral_out = variants.query("ONCOGENIC == 'Likely Neutral' and IN_FUNC_SITE == False")

# create contingency table 
observed_table = pd.DataFrame([
  [len(oncogenic_in), len(oncogenic_out)],
  [len(neutral_in), len(neutral_out)]],
  index= ["Oncogenic", "Likely Neutral"],
  columns=["Inside func_site", "Outside func_site"])

print("Contingency table (observed values):")
print(observed_table)

# run Chi-square 
chi2, p, dof, expected = chi2_contingency(observed_table)

# contingency table expected values (under H0) 
expected_table = pd.DataFrame(
    expected,
    columns=["Inside func_site", "Outside func_site"],
    index=["Oncogenic", "Likely Neutral"]
)

print("\nContingency table (expected values):")
print(expected_table.round(2))

# print results from Chi-square 
print("\nResults:")
print(f"Chi-square: {chi2:.3f}, p-value: {p:.4f}")

# cramers v (effect size for Chi-square) 
n = observed_table.values.sum()
k = min(observed_table.shape) - 1 
cramers_v = np.sqrt(chi2 / (n*k))
print(f"Cramér's V (effect size Chi-square): {cramers_v:.3f}")

# calculate the odds-ratio 
result = odds_ratio(observed_table)
ci = result.confidence_interval(confidence_level=0.95)

print(f"Odds-ratio: {result.statistic:.2f}")
print(f"95% CI: [{ci.low:.2f}, {ci.high:.2f}]")
print("-"*30)


# ------------------------------------------------------------
# Statistical Test: Feature by feature 
# ------------------------------------------------------------
# check each feature statistically using Chi-Square or Fisher 

print("Performing Chi-Square and Fisher test on feature types (one-by-one)..")

from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio   # ← added
from statsmodels.stats.multitest import multipletests

def analyze_func_sites(df):
    features = df["FEATURE_TYPE"].unique() 
    results = []

    for f in features: 

        if pd.isna(f): continue

        # create contingency table
        onc_in  = len(df[(df["ONCOGENIC"] == "Oncogenic") & (df["FEATURE_TYPE"] == f)])
        onc_out = len(df[(df["ONCOGENIC"] == "Oncogenic") & (df["FEATURE_TYPE"] != f)])
        neu_in  = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (df["FEATURE_TYPE"] == f)])
        neu_out = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (df["FEATURE_TYPE"] != f)])

        observed_table = [[onc_in, neu_in], [onc_out, neu_out]]

        # Odds ratio and 95% CI
        or_result = scipy_odds_ratio(observed_table)
        or_value = or_result.statistic
        ci = or_result.confidence_interval(confidence_level=0.95)
        ci_low, ci_high = ci.low, ci.high

        # select test based on number of observations in each cell
        # < 5 variants in one cell: Fisher, else: Chi-Square
        if min(onc_in, onc_out, neu_in, neu_out) < 5: 
            _, p = fisher_exact(observed_table)           
            test_used = "Fisher"
        else: 
            _, p, _, _ = chi2_contingency(observed_table)      
            test_used = "Chi-Square"
        
        # calculate effect size for chi-square: Cramer's V
        chi2, _, _, _ = chi2_contingency(observed_table)
        n = sum(sum(row) for row in observed_table)
        k = 1  # only for 2x2 table
        cramers_v = np.sqrt(chi2 / (n * k))

        results.append({
            "Feature":    f, 
            "p-value":    p, 
            "Odds_Ratio": or_value,
            "CI_95_low":  ci_low,
            "CI_95_high": ci_high,
            "Cramer's V": cramers_v if test_used == "Chi-Square" else np.nan,  
            "Test":       test_used
        })
    
    return pd.DataFrame(results) 

statistics_results = analyze_func_sites(expanded) 

# adjust for multiple testing (Benjamini-Hochberg)
statistics_results["p_adj"] = multipletests(statistics_results['p-value'], method='fdr_bh')[1]

# add information about significancy 
statistics_results["Significant"] = statistics_results["p_adj"].apply(lambda x: "Yes" if x < 0.05 else "No")

# sort and print results 
statistics_results = statistics_results.sort_values("p_adj")
print("Results one-by-one statistics functional sites:")
print(statistics_results.to_string(index=False))
print("-"*30)

print("\nFunctional site analysis complete!🥳\n")