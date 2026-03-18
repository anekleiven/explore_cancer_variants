
# ============================================================
# Statistics: Germline proximity in top oncogenic genes 
# ============================================================

"""
Script: stats_germline_proximity.py
Author: Ane Kleiven

Major output: 
  1. Run Mann-Whitney U and Rank-biserion correlation on variants with germline distance in top oncogenic genes. 
     Purpose: See if there is a difference in distances between oncogenic and neutral variants
     for the given genes. 
     Using the stats_func() 

""" 

print("-"*50)
print("Statistics Germline Proximity🔎🔢")
print("-"*50)

# Import libraries 
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio
from statsmodels.stats.multitest import multipletests

# Import statistics function 
from scripts.statistics.stats_function import stats_func

#-----------------------------------------------------------------------------------
# Load variant data
#-----------------------------------------------------------------------------------

# Load filtered germline data
print("Loading filtered variants with germline distances in top oncogenic genes (filtered and cleaned)..")
variants = pd.read_csv(
    "/home/anekl/git/master/explore_cancer_variants/output/germline_dist_filtered.tsv", 
    sep = "\t", 
    low_memory=False)

print(f"Loaded {len(variants)} variants from selected oncogenic genes.")
print("-"*30)

# Select feature 
features = ['Germline_Proximity'] 

#-----------------------------------------------------------------------------------
# Loop through each gene and run statistics 
#-----------------------------------------------------------------------------------

# Fetch gene names 
top_genes = variants['Hugo_Symbol'].unique()

print(f"Starting statistics for {len(top_genes)} genes...")

p_values_list = []
gene_names = []


# Loop through each gene 
for gene in top_genes:
    # Create a subset for the given gene
    gene_df = variants[variants['Hugo_Symbol'] == gene]
    
    onc = gene_df[gene_df["ONCOGENIC"] == "Oncogenic"]['Germline_Proximity']
    neu = gene_df[gene_df["ONCOGENIC"] == "Likely Neutral"]['Germline_Proximity']
    
    if len(onc) > 0 and len(neu) > 0:
        # Run statistics function 
        stats_func(gene_df, features, label=f"{gene}")
        print(f"{len(onc)} oncogenic variants for the {gene} gene.")
        print(f"{len(neu)} likely neutral variants for the {gene} gene.")

        _, p = mannwhitneyu(onc, neu, alternative="two-sided")
        p_values_list.append(p)
        gene_names.append(gene) 

    else:
        print("\n" + "-"*50)
        print(f"Statistics: {gene}")
        print("-"*50)
        print(f"\nGene: {gene} skipped (Onc: {len(onc)}, Neu: {len(neu)}).")
        

#-----------------------------------------------------------------------------------
# Benjamini Hochberg (multiple testing) 
#-----------------------------------------------------------------------------------

if p_values_list:
    # Run Benjamini Hochberg on all p-values
    reject, q_values, _, _ = multipletests(p_values_list, method='fdr_bh')

    # Create summary
    summary = pd.DataFrame({
        'Gene': gene_names,
        'P_value': p_values_list,
        'Q_value': q_values,
        'Significant_FDR': reject 
    }).sort_values("Q_value")

    print("\n" + "-"*50)
    print("Summary FDR-corrected testing:")
    print("-"*50)
    print(summary.to_string(index=False))
    print("-"*50)


print("\nGermline proximity statistics complete! 🧬")