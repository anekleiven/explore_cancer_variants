
# ============================================================
# Statistics: MAVE experiment for TP53 gene
# ============================================================

print("-"*50)
print("Statistics TP53 MAVE EXPERIMENT🔎🔢")
print("-"*50)

# -------------------------------------------
# Import libraries 
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
        default="/home/anekl/git/master/explore_cancer_variants/output/mave_urn_filtered.tsv",
        help="Path to the input file with variant data."
    )

    return parser.parse_args() 

args = getargs() 

# -------------------------------------------
# Load filtered MAVE data 
# -------------------------------------------

print(f"Loading MAVE variants from selected experiments:\n{args.variants}")

variants = pd.read_csv(
    args.variants, 
    sep = "\t", 
    low_memory=False)

print(f"Loaded {len(variants)} MAVE variants from selected experiments.")
print("-"*30)

# Select feature 
features = ['MaveDB_score']

#-----------------------------------------------------------------------------------

# ---TP53 gene---

# Filter to wanted urn 
urn = "urn:mavedb:00000068-a-1"
print("---TP53---")
print(f"\nFilter to variants in {urn} ")
variants_tp53 = variants[variants["MaveDB_urn"] == urn]
print(f"After filtering, the data contains {len(variants_tp53)} variants.")

print(f"Oncogenic TP53 variants: {len(variants_tp53[variants_tp53['ONCOGENIC'] == 'Oncogenic'])}")
print(f"Neutral TP53 variants: {len(variants_tp53[variants_tp53['ONCOGENIC'] == 'Likely Neutral'])}\n")

# Run statistics 
stats_func(variants_tp53, features, "urn:mavedb:00000068-a-1 TP53 variants") 
print("-"*30)

print("MAVE TP53 statistics complete🥳!\n") 