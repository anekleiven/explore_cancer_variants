# ============================================================
# Statistics: All variants 
# ============================================================

"""
Script: stats_all_variants.py
Author: Ane Kleiven

Main purpose: 

    Perform statistics on all variants
    using the stats_func() (see stats_function.py)

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
        required=True, 
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
# Define the features 
# -------------------------------------------

features = ["gnomAD_AF", "has_gnomAD_AF", "In_Hotspot", "IN_DOMAIN", "IN_FUNC_SITE", "is_null_variant", "is_null_var_tsg"]

# -------------------------------------------
# Run statistics 
# -------------------------------------------

stats_func(variants, features, "All variants") 
print("\n")