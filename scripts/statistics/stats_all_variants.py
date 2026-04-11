# ============================================================
# Statistics: All variants 
# ============================================================

print("-"*50)
print("Statistics all variants🔎🔢")
print("-"*50)

# import libraries 
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

# Add feature about gnomAD presence 
variants["has_gnomAD_AF"] = ( 
    (variants["gnomAD_AF"].notna()) & 
    (variants["gnomAD_AF"] != "NA") & 
    (variants["gnomAD_AF"] != "") 
)

# Define the features 
features = ["gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE", "is_null_variant", "is_null_var_tsg"]

# Run statistics function 
stats_func(variants, features, "All variants") 

