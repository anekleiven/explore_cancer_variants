# ============================================================
# Statistics: Functional sites (one-by-one)
# ============================================================

"""
Script purpose: 

Perform statistics on different functional site types from UniProt.

The discriminatory power of different functional sites are tested using 
Chi-Square test with Cramer's V/OR or Fisher's Exact test with OR.

p-values are adjusted for multiple testing using the Benjamini-Hochberg procedure. 

"""

print("-"*50)
print("Statistics Functional Sites (One-by-one)")
print("-"*50)

# -------------------------------------------
# import libraries 
# -------------------------------------------

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio 
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import argparse 

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Perform statistics on variant data."
    ) 

    parser.add_argument(
        "--variants", 
        type=Path, 
        required=False, 
        default="/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_tsg_og.tsv",
        help="Path to the input file with variant data."
    )

    parser.add_argument(
        "--alpha",
        type=float,
        required=False,
        default=0.05,
        help="Significance level for hypothesis testing (default: 0.05)."
    )

    return parser.parse_args() 

args = getargs() 

# -------------------------------------------
# define statistics function 
# -------------------------------------------

def analyse_func_sites(df, alpha=0.05):
    # Find all unique functional site types from the semicolon-separated column 
    all_features = set()
    df['FEATURE_TYPE'].dropna().str.split(';').apply(
        lambda x: [all_features.add(f.strip()) for f in x]
        )

    print(f"Significance level: α = {alpha}\n")

    results = []
    for f in sorted(list(all_features)):
        if not f: continue 

    for f in all_features:
        # Check original df 
        # Does this variant have feature x 
        # Split by semicolon, strip whitespace, and check if 'f' is in the resulting list
        has_f = df['FEATURE_TYPE'].str.split(';').apply(
        lambda x: f in [i.strip() for i in x] if isinstance(x, list) else False)
        
        onc_in  = len(df[(df["ONCOGENIC"] == "Oncogenic") & (has_f)])
        onc_out = len(df[(df["ONCOGENIC"] == "Oncogenic") & (~has_f)])
        neu_in  = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (has_f)])
        neu_out = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (~has_f)])

        observed_table = np.array([[onc_in, onc_out], [neu_in, neu_out]])
        total = onc_in + onc_out + neu_in + neu_out

        min_total = 10

        # skip functional sites with less than ten samples in total 
        if total < min_total:
            print(f"[{f}] Skipped (insufficient sample size: n={total})\n")
            print(f"{'Group':<12} | {'In Feature':<10} | {'Out Feature':<10} | {'Total':<6}")
            print("-" * 45)
            print(f"{'Oncogenic':<12} | {onc_in:<10} | {onc_out:<10} | {onc_in + onc_out:<6}")
            print(f"{'Neutral':<12} | {neu_in:<10} | {neu_out:<10} | {neu_in + neu_out:<6}")
            continue

        # skip functional sites with 0 values in a whole row or column 
        if np.any(observed_table.sum(axis=0) == 0) or np.any(observed_table.sum(axis=1) == 0):
            print(f"[{f}] Skipped: One or more categories have zero total observations.\n")
            continue

        # Select test based on expected values in each cell.
        # Expected < 5 in any cell: use Fisher's exact, else: Chi-Square.
        # Expected = 0 in any cell: use Fisher's exact, else: Chi-Square
        chi2, _, _, expected = chi2_contingency(observed_table)

        if expected.min() < 5 or any(0 in row for row in observed_table):
            _, p = fisher_exact(observed_table)
            test_used = "Fisher"
            cramers_v = np.nan
        else:
            _, p, _, _ = chi2_contingency(observed_table)
            test_used = "Chi-Square"
            n = onc_in + onc_out + neu_in + neu_out
            k = 1  # only for 2x2 table
            cramers_v = np.sqrt(chi2 / (n * k)) if n > 0 else 0

        # Odds ratio and 95% CI
        or_result = scipy_odds_ratio(observed_table)
        or_value = or_result.statistic
        ci = or_result.confidence_interval(confidence_level=0.95)

        results.append({
            "Feature":     f, 
            "p_value":    p, 
            "Odds_Ratio": or_value,
            "CI_95_low":  ci.low,
            "CI_95_high": ci.high,
            "Count_Onco_in": onc_in,
            "Count_Neu_in":  neu_in,
            "Cramer's V": cramers_v, 
            "Test":       test_used
        })
            
    # build dataframe and adjust for multiple testing (Benjamini-Hochberg)
    results_df = pd.DataFrame(results)
    results_df["p_adj"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
    results_df["Significant"] = results_df["p_adj"].apply(lambda x: "Yes" if x < alpha else "No")
    results_df = results_df.sort_values("p_adj").reset_index(drop=True)

    return results_df

# ------------------------------------------------
# Functional site analysis: all variants
# ------------------------------------------------

print("-"*40) 
print("All variants") 
print("-"*40)

# -------------------------------------------
# Load and clean variant data 
# -------------------------------------------

print(f"Loading variant file:\n{args.variants}")
variants = pd.read_csv(
    args.variants, 
    sep = "\t", 
    low_memory=False)

print(f"Loaded {len(variants)} variants.")

# clean column
variants['FEATURE_TYPE'] = variants['FEATURE_TYPE'].str.strip()

# -------------------------------------------
# Perform statistics 
# -------------------------------------------

print("Performing statistics on feature types (one-by-one)..\n")

# filter to rows with a functional site annotation 
df_with_func_sites = variants.dropna(subset=["FEATURE_TYPE"])
stats_all = analyse_func_sites(df_with_func_sites, alpha=args.alpha)

print("Results all variants:")
print(stats_all.to_string(index=False))
print("-"*30)

# ------------------------------------------------
# Functional site analysis: top 10 oncogenic genes 
# ------------------------------------------------

print("-"*40) 
print("Top 10 oncogenic genes") 
print("-"*40)

top_10_onco_genes = (
    variants[variants['ONCOGENIC'] == 'Oncogenic']['Hugo_Symbol']
    .value_counts()
    .head(10)
    .index.tolist()
)

print(f"Top 10 genes: {top_10_onco_genes}")

# extract top oncogenic variants
df_top_genes = df_with_func_sites[df_with_func_sites["Hugo_Symbol"].isin(top_10_onco_genes)]
stats_top_10 = analyse_func_sites(df_top_genes, alpha=args.alpha) 

print("Results top 10 oncogenic genes:")
print(stats_top_10.to_string(index=False))
print("-"*30)


print("Statistics on functional sites complete!🥳")