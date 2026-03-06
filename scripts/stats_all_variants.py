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

stats_func(variants, features, "All variants") 

##-------------------------------------------
# Statistics MAVEs LOF vs GOF
#-------------------------------------------

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

# check distribution of lof and gof 
effect_summary = variants.groupby("MUTATION_EFFECT").size()
print(effect_summary, "\n")

# create lof and gof labels
lof_labels = ["Likely Loss-of-function", "Loss-of-function"]
gof_labels = ["Likely Gain-of-function", "Gain-of-function"]

# define the data
lof_variants = variants[variants["MUTATION_EFFECT"].isin(lof_labels)]
gof_variants = variants[variants["MUTATION_EFFECT"].isin(gof_labels)]

lof_oncogenic = lof_variants[lof_variants["ONCOGENIC"] == "Oncogenic"]["MaveDB_score"].dropna()
gof_oncogenic = gof_variants[gof_variants["ONCOGENIC"] == "Oncogenic"]["MaveDB_score"].dropna()

# check number of variants 
print(f"LoF oncogenic variants with MAVE score: {len(lof_oncogenic)}")
print(f"GoF oncogenic variants with MAVE score: {len(gof_oncogenic)}")

# perform statistics 
stat, p = mannwhitneyu(lof_oncogenic, gof_oncogenic, alternative="two-sided")
n1, n2 = len(lof_oncogenic), len(gof_oncogenic)
r = (2 * stat) / (n1 * n2) - 1
probability = (1 + r) / 2

print(f"\nMann-Whitney U: {stat:.3f}, p-value: {p:.4f}")
print(f"Rank-biserial r: {r:.3f} | P(LoF > GoF): {probability*100:.2f}%")
print(f"{'Reject H₀: MAVE scores differ between GoF and LoF.' if p < 0.05 else 'Failed to reject H₀.'}")