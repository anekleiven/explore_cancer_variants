
# ============================================================
# Statistics: Top 10 oncogenic genes 
# ============================================================

"""
Script purpose: 

Perform statistics on top 10 oncogenic genes (by oncogenic variant counts)
using the stats_func() 

Features included: 
"gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE", "is_null_var_tsg", "is_null_variant"

The statistics function performs Mann-Whitney U test with rank-biserial correlation on continuous features, 
Chi-Square test with Cramer's V/OR or Fisher's Exact test with OR on categorical features. 

p-values are adjusted for multiple testing using the Benjamini-Hochberg procedure. 

"""

# -------------------------------------------
# import libraries 
# -------------------------------------------

import pandas as pd
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
        default="/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_tsg_og.tsv",
        help="Path to the input file with variant data."
    )

    return parser.parse_args() 

args = getargs() 

# -------------------------------------------
# Load variant data 
# -------------------------------------------

print(f"Loading variant file:\n{args.variants}")
variants = pd.read_csv(
    args.variants, 
    sep = "\t", 
    low_memory=False)

print(f"Loaded {len(variants)} variants.")
print("-"*30)

# -------------------------------------------
# Add feature about gnomAD presence 
# -------------------------------------------

variants["has_gnomAD_AF"] = ( 
    (variants["gnomAD_AF"].notna()) & 
    (variants["gnomAD_AF"] != "NA") & 
    (variants["gnomAD_AF"] != "") & 
    (variants["gnomAD_AF"] > 0)
)

# -------------------------------------------
# Extract top oncogenic genes 
# -------------------------------------------

oncogenic_variants = variants[variants['ONCOGENIC'] == 'Oncogenic']

oncogenic_genes = ( 
  oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
)

top_10_onco_genes = oncogenic_genes["Hugo_Symbol"].head(10).tolist() 
print("Top 10 oncogenic genes (by variant count):")
print(top_10_onco_genes)

# -------------------------------------------
# Extract top oncogenic variants
# -------------------------------------------

top_10_gene_variants = variants[variants["Hugo_Symbol"].isin(top_10_onco_genes)]


# -------------------------------------------
# Define the features 
# -------------------------------------------

features = ["gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE", "is_null_var_tsg", "is_null_variant"]

# -------------------------------------------
# Run statistics function 
# -------------------------------------------

stats_func(top_10_gene_variants, features, "Top 10 oncogenic genes")
print("\n")