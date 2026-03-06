
#========================================================================
# Statistics Cancer Variants
#========================================================================

"""
Script: stats_function.py
Author: Ane Kleiven

Description: 
Statistics function for cancer variants. 

Purpose: 
Perform statistics on variant features such as: 
  gnomAD_AF (population frequency data)
  has_gnomAD_AF (presence of population frequency data)
  cancer hotspots
  protein domains
  functional sites 
  germline proximity 
  MAVE scores 

Statistical tests: 
*   Mann Whitney U for continuous data (Rank-biserion correlation for effect size) 
*   Chi-Square or Fisher test for categorical data (Cramer's V and Odds-ratio for effect size). 
*   Benjamini Hochberg to adjust for multiple testing 

"""

# import libraries 
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio
from statsmodels.stats.multitest import multipletests

def stats_func(df, features, label="Dataset"):
    print("\n" + "-"*50)
    print(f"Statistics: {label}")
    print("-"*50 + "\n")

    results = []

    for f in features:

        if f in ["gnomAD_AF", "Germline_Proximity", "MaveDB_score"]: 
            oncogenic = df[df["ONCOGENIC"] == "Oncogenic"][f].dropna()
            neutral = df[df["ONCOGENIC"] == "Likely Neutral"][f].dropna()

            if f == "gnomAD_AF":
                oncogenic = oncogenic[oncogenic > 0]
                neutral = neutral[neutral > 0]

            if len(oncogenic) == 0 or len(neutral) == 0:
                print(f"[{f}] Skipped (not enough data)\n")
                continue

            # perform Mann-Whitney U test 
            stat, p = mannwhitneyu(oncogenic, neutral, alternative="two-sided")
            test_used = "Mann-Whitney U"
            # calculate rank-biserial correlation 
            # (effect size for mann-whitney u) 
            n1 = len(oncogenic)
            n2 = len(neutral)
            r = (2 * stat) / (n1 * n2) - 1

            # calculate probability 
            probability = (1+r)/2 

            print(f"[{f}]")
            print(f"Test: {test_used}")
            print(f"Mann-Whitney U: {stat:.3f}, p-value: {p:.4f}")
            print(f"Rank-biserial r: {r:.3f} | P(oncogenic > neutral): {probability*100:.2f}%")
            results.append({"feature": f, "test": "Mann-Whitney", "p_value": p})
            print(f"{'Reject H₀: distributions differ.' if p < 0.05 else 'Failed to reject H₀.'}\n")


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
            cramers_v = None 

            if min(onc_in, onc_out, neu_in, neu_out) < 5: 
                _, p = fisher_exact(observed_table)           
                test_used = "Fisher"
            else: 
                chi2, p, _, _ = chi2_contingency(observed_table)     
                test_used = "Chi-Square"
                n = sum(sum(row) for row in observed_table)
                k = 1 # only for 2x2 tables
                cramers_v = np.sqrt(chi2 / (n * k))

            if test_used == "Fisher":
                print(f"[{f}]")
                print(f"Test: {test_used}")
                print(f"OR: {or_result.statistic:.3f} (95% CI: {ci.low:.3f}–{ci.high:.3f}) | p-value: {p:.4f}")
                print(f"{'Reject H₀: association between oncogenicity and ' + f + '.' if p < 0.05 else 'Failed to reject H₀.'}\n")
            else:
                print(f"[{f}]")
                print(f"Test: {test_used}")
                print(f"OR: {or_result.statistic:.3f} (95% CI: {ci.low:.3f}–{ci.high:.3f})")
                print(f"p-value: {p:.4f} | Cramer's V: {cramers_v:.3f}")
                print(f"{'Reject H₀: association between oncogenicity and ' + f + '.' if p < 0.05 else 'Failed to reject H₀.'}\n")

            results.append({"feature": f, "test": test_used, "p_value": p})
                
     
    results_df = pd.DataFrame(results)
    results_df["p_value"] = results_df["p_value"].round(4)

    if len(features) > 1:
        _, q_values, _, _ = multipletests(results_df["p_value"], method="fdr_bh")
        results_df["q_value"] = q_values.round(4)


        print(f"{'-'*50}")
        print("FDR-corrected results (Benjamini-Hochberg)")
        print(f"{'-'*50}")
        print(results_df.to_string(index=False))