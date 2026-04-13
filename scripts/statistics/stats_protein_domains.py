# ============================================================
# Statistics: Protein Domains (Top 10) - One-by-one
# ============================================================

print("-" * 50)
print("Statistics Protein Domains (Top 10)")
print("-" * 50)

# -------------------------------------------
# Import libraries 
# -------------------------------------------

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio 
from statsmodels.stats.multitest import multipletests
import argparse
from pathlib import Path

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

    return parser.parse_args() 

args = getargs() 

# -------------------------------------------
# Create stats function protein domains 
# -------------------------------------------

def analyze_top_domains(df, n_top=10):
    """
    Analyzing n top protein domains 
    Using regex to match domain names in semicolon lists 
    """
    # Find top domains (by variant counts)
    all_counts = df['DOMAIN_NAME'].dropna().str.split(';').explode().str.strip().value_counts()
    top_domains = all_counts.head(n_top).index.tolist()
    
    print(f"Analyzing the top {len(top_domains)} most frequent domains...")
    
    results = []
    for d in top_domains:
        if not d: continue
        
        # Regex-matching: Make sure 'Pkinase' do not match 'Pkinase_Tyr' by mistake
        pattern = rf"(^|;)\s*{d}\s*($|;)"
        has_d = df['DOMAIN_NAME'].str.contains(pattern, na=False, regex=True)
        
        # Contingency table
        onc_in  = len(df[(df["ONCOGENIC"] == "Oncogenic") & (has_d)])
        onc_out = len(df[(df["ONCOGENIC"] == "Oncogenic") & (~has_d)])
        neu_in  = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (has_d)])
        neu_out = len(df[(df["ONCOGENIC"] == "Likely Neutral") & (~has_d)])

        observed_table = [[onc_in, neu_in], [onc_out, neu_out]]
        
        # Odds ratio 
        or_result = scipy_odds_ratio(observed_table)
        or_value = or_result.statistic
        ci = or_result.confidence_interval()

        # Statistical testing
        if min(onc_in, onc_out, neu_in, neu_out) < 5: 
            _, p = fisher_exact(observed_table)           
            test_used = "Fisher"
        else: 
            _, p, _, _ = chi2_contingency(observed_table)      
            test_used = "Chi-Square"
        
        results.append({
            "Domain":      d, 
            "p_value":     p, 
            "Odds_Ratio":  or_value,
            "CI_95_low":   ci.low,
            "CI_95_high":  ci.high,
            "Count_Onco":  onc_in,
            "Count_Neu":   neu_in,
            "Test":        test_used
        })
    
    # Adjust for multiple testing 
    results_df = pd.DataFrame(results)
    if not results_df.empty:
        results_df["p_adj"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
        results_df["Significant"] = results_df["p_adj"].apply(lambda x: "Yes" if x < 0.05 else "No")
        results_df = results_df.sort_values("p_adj").reset_index(drop=True)

    return results_df

# -------------------------------------------
# Perform statistics
# -------------------------------------------

# Import variant data 
variants = pd.read_csv(
  args.variants,
  sep="\t",
  low_memory=False)

df_domains = variants.dropna(subset=["DOMAIN_NAME"]).copy()
df_domains['DOMAIN_NAME'] = df_domains['DOMAIN_NAME'].str.strip()

# Analyze for all variants
print("\nResults: Top 10 Domains (All Variants)")
stats_domains_all = analyze_top_domains(df_domains, n_top=10)
print(stats_domains_all.to_string(index=False))

# Analyze for top 10 oncogenes

top_10_onco_genes = (
    variants[variants['ONCOGENIC'] == 'Oncogenic']['Hugo_Symbol']
    .value_counts()
    .head(10)
    .index.tolist()
)

print("\nResults: Top 10 Domains (Top 10 Oncogenic Genes)")
df_top_genes_domains = df_domains[df_domains["Hugo_Symbol"].isin(top_10_onco_genes)]
stats_domains_top = analyze_top_domains(df_top_genes_domains, n_top=10)
print(stats_domains_top.to_string(index=False))


print("\nProtein domain statistics complete!🥳\n")