
# ============================================================
# Statistics: Top 10 oncogenic genes (by variant counts)
# ============================================================

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
# Add feature about gnomAD presence 
# -------------------------------------------

top_10_gene_variants["has_gnomAD_AF"] = ( 
    (top_10_gene_variants["gnomAD_AF"].notna()) & 
    (top_10_gene_variants["gnomAD_AF"] != "NA") & 
    (top_10_gene_variants["gnomAD_AF"] != "") 
)

# -------------------------------------------
# Define the features 
# -------------------------------------------

features = ["gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE", "is_null_variant", "is_null_var_tsg"]

# -------------------------------------------
# Run statistics function 
# -------------------------------------------

stats_func(top_10_gene_variants, features, "Top 10 oncogenic genes")
print("\n")