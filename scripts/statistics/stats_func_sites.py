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
    # Find all unique functional site types from the semicolon-separated column 
    all_features = set()
    df['FEATURE_TYPE'].dropna().str.split(';').apply(
        lambda x: [all_features.add(f.strip()) for f in x]
        )
    
    results = []
    for f in sorted(list(all_features)):
        if not f: continue 

    for f in all_features:
        # Check original df 
        # Does this variant have feature x 
        pattern = rf"(^|;)\s*{f}\s*($|;)"
        has_f = df['FEATURE_TYPE'].str.contains(pattern, na=False, regex=True)
        
        onc_in  = len(df[(df["ONCOGENIC"] == "Oncogenic") & (has_f)])
        onc_out = len(df[(df["ONCOGENIC"] == "Oncogenic") & (~has_f)])
        neu_in  = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (has_f)])
        neu_out = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (~has_f)])

        observed_table = [[onc_in, neu_in], [onc_out, neu_out]]

        # skip feature if there is not enough samples
        if (onc_in + neu_in) == 0: continue
        
        # Odds ratio and 95% CI
        or_result = scipy_odds_ratio(observed_table)
        or_value = or_result.statistic
        ci = or_result.confidence_interval(confidence_level=0.95)

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
            n = onc_in + onc_out + neu_in + neu_out
            k = 1  # only for 2x2 table
            cramers_v = np.sqrt(chi2 / (n * k)) if n > 0 else 0 
        else:
            cramers_v = np.nan

        results.append({
            "Feature":     f, 
            "p_value":    p, 
            "Odds_Ratio": or_value,
            "CI_95_low":  ci.low,
            "CI_95_high": ci.high,
            "Count_Onco": onc_in,
            "Count_Neu":  neu_in,
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

# clean column for regex 
variants['FEATURE_TYPE'] = variants['FEATURE_TYPE'].str.strip()

print(f"Loaded {len(variants)} variants.")

# filter to rows with a functional site annotation 
df_with_func_sites = variants.dropna(subset=["FEATURE_TYPE"])
stats_all = analyze_func_sites(df_with_func_sites)

print("Performing statistics on feature types (one-by-one)..\n")

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
stats_top_10 = analyze_func_sites(df_top_genes) 

print("Results top 10 oncogenic genes:")
print(stats_top_10.to_string(index=False))
print("-"*30)


print("Statistics on functional sites complete!🥳")