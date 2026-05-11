
# ============================================================
# Statistics: Germline proximity in top oncogenic genes 
# ============================================================

"""
Script purpose: 

Perform statistics on germline proximity in top oncogenic genes. 

The discriminatory power of germline proximity in selected genes is tested
using Mann-Whitney U with rank-biserial correlation. 

p-values are adjusted for multiple testing using the Benjamini-Hochberg procedure. 

"""

print("-"*50)
print("Statistics Germline Proximity🔎🔢")
print("-"*50)

# -------------------------------------------
# Import libraries 
# -------------------------------------------

import pandas as pd
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import argparse

# -------------------------------------------
# Import statistics function 
# -------------------------------------------

from scripts.statistics.stats_function import stats_func

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
        default="/home/anekl/git/master/explore_cancer_variants/output/germline_dist_filtered.tsv",
        help="Path to the input file with variant data."
    )

    return parser.parse_args() 

args = getargs() 

# -------------------------------------------
# Load variant data
# -------------------------------------------

print(f"Loading variants with germline distances in top oncogenic genes:\n{args.variants}")
variants = pd.read_csv(
    args.variants, 
    sep = "\t", 
    low_memory=False)

print(f"Loaded {len(variants)} variants from selected oncogenic genes.")
print("-"*30)

# Select feature 
features = ['Germline_Proximity'] 

# -------------------------------------------
# Loop through each gene and run statistics 
# -------------------------------------------

# Fetch gene names 
top_genes = variants['Hugo_Symbol'].unique()

p_values_list = []
gene_names = []

# Loop through each gene 
for gene in top_genes:
    # Create a subset for the given gene
    gene_df = variants[variants['Hugo_Symbol'] == gene]
    
    onc = gene_df[gene_df["ONCOGENIC"] == "Oncogenic"]['Germline_Proximity']
    neu = gene_df[gene_df["ONCOGENIC"] == "Likely Neutral"]['Germline_Proximity']
    
    res = stats_func(gene_df, features, label=f"{gene}")

    if res is not None:
        for _, row in res.iterrows():
            p_values_list.append(row['p_value'])
            gene_names.append(gene)

# -------------------------------------------
# Benjamini Hochberg (multiple testing) 
# -------------------------------------------

if p_values_list:
    reject, q_values, _, _ = multipletests(p_values_list, method='fdr_bh')
    
    summary = pd.DataFrame({
        'Gene': gene_names,
        'P_value': p_values_list,
        'P_adj': q_values,
        'Significant_FDR': reject 
    }).sort_values("P_adj")

    # round values
    summary["P_value"] = summary["P_value"].round(4)
    summary["P_adj"] = summary["P_adj"].round(4)

    print("\n" + "="*60)
    print("OVERALL GENE SUMMARY (FDR-CORRECTED)")
    print("="*60)
    print(summary.to_string(index=False))
    print("="*60)

print("\nGermline proximity statistics complete!🧬\n")