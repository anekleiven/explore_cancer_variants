# ============================================================
# Statistics: MAVE experiments
# ============================================================

"""
Script: stats_maves.py
Author: Ane Kleiven

Major output: 
  1. Run Mann-Whitney U and Rank-biserion correlation on the urn:mavedb:00000068-a-1 experiment
     Purpose: See if there is a difference in MAVE scores between oncogenic and neutral variants
     for the given MAVE experiment in gene TP53 
     Using the stats_func() 
  2. Run statistics on PTEN, KRAS and BRCA1 (note: insufficient data)

""" 

print("-"*50)
print("Statistics MAVE EXPERIMENT🔎🔢")
print("-"*50)

# Import libraries 
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency, fisher_exact
from scipy.stats.contingency import odds_ratio as scipy_odds_ratio
from statsmodels.stats.multitest import multipletests

# Import statistics function 
from stats_function import stats_func

#-----------------------------------------------------------------------------------

# Load filtered MAVE data 
print("Loading MAVE variants from selected experiments (filtered and cleaned)..")
variants = pd.read_csv(
    "/home/anekl/git/master/explore_cancer_variants/output/mave_urn_filtered.tsv", 
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


#-----------------------------------------------------------------
# Extra: Statistics on other MAVE experiments 
# (NB: insufficient data)
#-----------------------------------------------------------------

# ---PTEN gene--- 

# Filter to wanted urn 
urn = "urn:mavedb:00000013-a-1"
print("---PTEN---")
print(f"Filter to data to variants in {urn} ")
variants_PTEN = variants[variants["MaveDB_urn"] == urn]
print(f"After filtering, the data contains {len(variants_PTEN)} variants.")

print(f"Oncogenic PTEN variants: {len(variants_PTEN[variants_PTEN['ONCOGENIC'] == 'Oncogenic'])}")
print(f"Neutral PTEN variants: {len(variants_PTEN[variants_PTEN['ONCOGENIC'] == 'Likely Neutral'])}\n")

# Run statistics 
stats_func(variants_PTEN, features, "urn:mavedb:00000013-a-1 PTEN variants") 
print("-"*30)

#-----------------------------------------------------------------------------------

# ---KRAS gene---

# Filter to wanted urn 
urn = "urn:mavedb:00000115-a-7"
print("---KRAS---")
print(f"Filter to data to variants in {urn} ")
variants_KRAS = variants[variants["MaveDB_urn"] == urn]
print(f"After filtering, the data contains {len(variants_KRAS)} variants.")

print(f"Oncogenic KRAS variants: {len(variants_KRAS[variants_KRAS['ONCOGENIC'] == 'Oncogenic'])}")
print(f"Neutral KRAS variants: {len(variants_KRAS[variants_KRAS['ONCOGENIC'] == 'Likely Neutral'])}\n")

# Run statistics 
stats_func(variants_KRAS, features, "urn:mavedb:00000115-a-7 KRAS variants") 
print("-"*30)


#-----------------------------------------------------------------------------------

# ---BRCA1 gene---

# Filter to wanted urn 
urn = "urn:mavedb:00000081-a-2"
print("---BRCA1---")
print(f"Filter to data to variants in {urn} ")
variants_BRCA1 = variants[variants["MaveDB_urn"] == urn]
print(f"After filtering, the data contains {len(variants_BRCA1)} variants.")

print(f"Oncogenic BRCA1 variants: {len(variants_BRCA1[variants_BRCA1['ONCOGENIC'] == 'Oncogenic'])}")
print(f"Neutral BRCA1 variants: {len(variants_BRCA1[variants_BRCA1['ONCOGENIC'] == 'Likely Neutral'])}\n")

# Run statistics 
stats_func(variants_BRCA1, features, "urn:mavedb:00000081-a-2 BRCA1 variants") 
print("-"*30)

print("MAVE statistics complete🥳!")