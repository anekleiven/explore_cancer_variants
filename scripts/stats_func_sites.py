# ============================================================
# Statistics: Functional sites (one-by-one)
# ============================================================

print("-"*50)
print("Statistics Functional Sites (One-by-one)")
print("-"*50)

# import libraries
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio 
from statsmodels.stats.multitest import multipletests

# define statistics function 
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
        if test_used == "Chi-Square":
            chi2, _, _, _ = chi2_contingency(observed_table)
            n = sum(sum(row) for row in observed_table)
            k = 1  # only for 2x2 table
            cramers_v = np.sqrt(chi2 / (n * k))
        else:
            cramers_v = np.nan

        results.append({
            "Feature":    f, 
            "p_value":    p, 
            "Odds_Ratio": or_value,
            "CI_95_low":  ci_low,
            "CI_95_high": ci_high,
            "Cramer's V": cramers_v, 
            "Test":       test_used
        })
    
    # build dataframe and adjust for multiple testing (Benjamini-Hochberg)
    results_df = pd.DataFrame(results)
    results_df["p_adj"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
    results_df["Significant"] = results_df["p_adj"].apply(lambda x: "Yes" if x < 0.05 else "No")
    results_df = results_df.sort_values("p_adj").reset_index(drop=True)

    return results_df

# ------------------------------------------------
# Functional site analysis: all variants
# ------------------------------------------------

print("-"*40) 
print("All variants") 
print("-"*40)

# import variant data 
variants = pd.read_csv(
  "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
  sep="\t",
  low_memory=False)

print(f"Loaded {len(variants)} variants.")

# expand data (one row per functional site) 
expanded = (
    variants
    .dropna(subset=["FEATURE_TYPE"])
    .assign(FEATURE_TYPE=lambda df: df["FEATURE_TYPE"].str.split(";"))
    .explode("FEATURE_TYPE")
)

expanded["FEATURE_TYPE"] = expanded["FEATURE_TYPE"].str.strip()
print(f"\nExpanded to {len(expanded):,} feature-variant rows.\n")

print("Performing statistics on feature types (one-by-one)..\n")

statistics_results = analyze_func_sites(expanded)
print("Results all variants:")
print(statistics_results.to_string(index=False))
print("-"*30)

# ------------------------------------------------
# Functional site analysis: top 10 oncogenic genes 
# ------------------------------------------------

print("-"*40) 
print("Top 10 oncogenic genes") 
print("-"*40)

oncogenic_variants = expanded[expanded['ONCOGENIC'] == 'Oncogenic']

oncogenic_genes = ( 
  oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
)

top_10_onco_genes = oncogenic_genes["Hugo_Symbol"].head(10).tolist() 
print("Top 10 oncogenic genes (by variant count):")
print(top_10_onco_genes, "\n")

# extract top oncogenic variants
top_10_gene_variants = expanded[expanded["Hugo_Symbol"].isin(top_10_onco_genes)]

statistics_results = analyze_func_sites(top_10_gene_variants)
print("Results top 10 oncogenic genes:")
print(statistics_results.to_string(index=False))
print("-"*30)


print("Statistics on functional sites complete!🥳")