
#========================================================================
# Neutral variants ClinVar
#========================================================================

"""
Script: neutral_dataset_clinvar.py
Author: Ane Kleiven

This script extracts neutral variants from the ClinVar data set. 
Only includes Gene IDs from the GENIE variant file. 

Filtered on the following: 
  Assembly = GRCh37
  Origin = Germline
  ClinicalSignificance = Benign, Likely Benign 
  Review status >= 1 review star 
  
"""

# Import libraries 
import pandas as pd 

# -----------------------------------------------
# Read variant files 
# -----------------------------------------------

# Load somatic variant file 
print("\nLoading somatic variants..")
somatic_var = pd.read_csv("/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_tsg_og.tsv.tsv",
    sep="\t",
    low_memory=False)
print(f"\nLoaded {len(somatic_var)} somatic variants.\n")

# Extract unique gene IDs 
somatic_gene_ids = set(somatic_var["Entrez_Gene_Id"].dropna().unique())

# Load ClinVar variant file 
print("Loading ClinVar variant file..\n")

clinvar_path = "/home/anekl/git/master/cancer_variants_annotation_pipeline/data/variant_summary.txt/variant_summary.txt"

# Columns to read: 
germline_cols = [
  "Assembly",
  "Origin",
  "Name",
  "GeneID",
  "Chromosome", 
  "PositionVCF",
  "ReferenceAlleleVCF",
  "AlternateAlleleVCF",
  "GeneSymbol", 
  "ClinicalSignificance",
  "ReviewStatus"
]

# Read file in chunks to save memory 
chunks =  []

for chunk in pd.read_csv(
  clinvar_path, 
  sep="\t", 
  usecols=germline_cols, 
  chunksize=100000,
  low_memory=False
):
  chunks.append(chunk[chunk["GeneID"].isin(somatic_gene_ids)])

clinvar_variants = pd.concat(chunks, ignore_index=True)

print(f"Loaded {len(clinvar_variants):,} variants from ClinVar.\n")

# -----------------------------------------------
# Filter ClinVar variants 
# -----------------------------------------------

# Extract variants using the following filters: 
#   Assembly == GRCh37
#   Oncogenicity == Benign and Likely Benign 
#   Origin == germline
#   ReviewStatus >= 1 review star

review_statuses = [
    "practice guideline",
    "reviewed by expert panel",
    "criteria provided, multiple submitters, no conflicts",
    "criteria provided, single submitter"
]

# Value check before filtering 
print(clinvar_variants["ReviewStatus"].value_counts()) 
print(clinvar_variants["ClinicalSignificance"].value_counts()) 
print(clinvar_variants["Origin"].value_counts())

neutral_clinvar = clinvar_variants[
    (clinvar_variants["Assembly"] == "GRCh37") &
    (clinvar_variants["Origin"] == "germline") & 
    (clinvar_variants["ClinicalSignificance"].isin(["Benign", "Likely Benign"])) & 
    (clinvar_variants["ReviewStatus"].isin(review_statuses))
].copy()

print(f"After filtering, the dataset consists of {len(neutral_clinvar)} neutral variants from ClinVar")
print(neutral_clinvar.head())

# Extract HGVSp from 'Name' column 
neutral_clinvar["HGVSp"] = neutral_clinvar["Name"].str.extract(r'\((p\.\w+)\)')

# Remove variants with NA HGVSp 
neutral_clinvar = neutral_clinvar[neutral_clinvar["HGVSp"].notna()].copy()

print(f"Number of neutral variants with protein mapping: {len(neutral_clinvar)}")
print(neutral_clinvar[['Name', 'HGVSp']].head())


# save file 
output_path = "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/neutral_clinvar_filtered.tsv"
print(f"\nSaving {len(neutral_clinvar)} variants to {output_path}...")
neutral_clinvar.to_csv(output_path, sep="\t", index=False)


print("\nNeutral ClinVar data complete!🎉")
