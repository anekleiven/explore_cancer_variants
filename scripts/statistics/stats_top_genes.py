
# ============================================================
# Statistics: Top 10 oncogenic genes (by variant counts)
# ============================================================

print("-"*50)
print("Statistics Top 10 Oncogenic Genes🔎🔢")
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

# Load variant data 
print("Loading variants..")
variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_tsg_og.tsv", 
    sep = "\t", 
    low_memory=False)

print(f"Loaded {len(variants)} variants.")
print("-"*30)

# Extract top oncogenic genes 
oncogenic_variants = variants[variants['ONCOGENIC'] == 'Oncogenic']

oncogenic_genes = ( 
  oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
)

top_10_onco_genes = oncogenic_genes["Hugo_Symbol"].head(10).tolist() 
print("Top 10 oncogenic genes (by variant count):")
print(top_10_onco_genes)

# Extract top oncogenic variants
top_10_gene_variants = variants[variants["Hugo_Symbol"].isin(top_10_onco_genes)]

# Add feature about gnomAD presence 
top_10_gene_variants["has_gnomAD_AF"] = ( 
    (top_10_gene_variants["gnomAD_AF"].notna()) & 
    (top_10_gene_variants["gnomAD_AF"] != "NA") & 
    (top_10_gene_variants["gnomAD_AF"] != "") 
)

# Define the features 
features = ["gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE"]

# Run statistics function 
stats_func(top_10_gene_variants, features, "Top 10 oncogenic genes")