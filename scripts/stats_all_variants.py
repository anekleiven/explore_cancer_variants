# ============================================================
# Statistics: All variants 
# ============================================================

print("-"*50)
print("Statistics Top 10 Oncogenic Genes🔎🔢")
print("-"*50)

# import libraries 
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio
from statsmodels.stats.multitest import multipletests

# import statistics function 
from stats_function import stats_func

# load variant data 
print("Loading variants..")
variants = pd.read_csv(
    "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv", 
    sep = "\t", 
    low_memory=False)

print(f"Loaded {len(variants)} variants.")
print("-"*30)

# add feature about gnomAD presence 
variants["has_gnomAD_AF"] = ( 
    (variants["gnomAD_AF"].notna()) & 
    (variants["gnomAD_AF"] != "NA") & 
    (variants["gnomAD_AF"] != "") 
)

# define the features 
features = ["gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE", "Germline_Proximity", "MaveDB_score"]

if __name__ == "__main__":
  stats_func(variants, features, "Top 10 oncogenic genes")