# ====================================================================
# Oncogene and tumor suppressor gene analysis 
# ====================================================================

"""
Script: tsg_og.py
Author: Ane Kleiven

This script performs a multi-step analysis to explore how somatic cancer
variant classes distribute oncogenes and tumor suppressor genes. 

Major outputs:
    1. Distribution of variants in tumor suppressor genes and oncogenes (Oncogenic and Likely Neutral variants) 
    2. Prevalence of null variants by oncogenicity class 
    3. Prevelance of null variants in tsg by oncogenicity class 
    4. Genes with oncogenic null variants in tsg 


All plots are saved in:
    plots/tsg_og

"""

print("-"*50)
print("TSG/OG Analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import os

# -------------------------------------------
# Argparse function for user input file paths
# -------------------------------------------

def getargs(): 
    parser = argparse.ArgumentParser(
        description="Explore oncogenes and tumor suppressor genes in variant data."
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

# ------------------------------------------------------------
# Create output directory 
# ------------------------------------------------------------

save_dir = "plots/tsg_og"
os.makedirs(save_dir, exist_ok=True) 

# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------

print(f"Loading variant file:\n{args.variants}")

variants = pd.read_csv(
    args.variants,
    sep="\t",
    low_memory=False
)

print(f"Loaded {len(variants):,} variants.")
print("-"*50)

# ------------------------------------------------------------
# Filter to oncogenic and likely neutral variants
# ------------------------------------------------------------

print("Filtering to wanted oncogenic classes..")

wanted_classes = ["Oncogenic", "Likely Neutral"]

variants_filtered = variants[variants["ONCOGENIC"].isin(wanted_classes)]

print("-"*50)

# ------------------------------------------------------------
# Create gene type column 
# ------------------------------------------------------------

print("Creating gene type column (oncogene/tumor suppressor gene)..")

def get_gene_type(row):
    if row["ncg_oncogene"] and row["ncg_tsg"]: return "Oncogene and TSG"
    if row["ncg_oncogene"]: return "Oncogene"
    if row["ncg_tsg"]: return "TSG"
    return "Other"

variants_filtered["gene_category"] = variants_filtered.apply(get_gene_type, axis=1)

print("Example output of new column:")
print(variants_filtered["gene_category"].head(10))


og_tsg_summary = variants_filtered.groupby(["gene_category", "ONCOGENIC"]).size().unstack(fill_value=0)
print("\nSummary OG/TSG:")
print(og_tsg_summary)
print("-"*50)

# ------------------------------------------------------------
# Plot gene category distributions for 
# oncogenic and likely neutral variants 
# ------------------------------------------------------------

# plotting summary 

print("Plotting gene category distributions for oncogenic and likely neutral variants..")

og_tsg_summary.plot(kind='bar', figsize=(8, 5), color=['#88aed1', '#c4314a'], edgecolor="0.1", linewidth=0.3)

plt.title("Distribution of Variant Classes per Gene Category", fontsize=14)
plt.ylabel("Number of Variants", fontsize=12)
plt.xlabel("Gene Category", fontsize=12)
plt.xticks(rotation=45, fontsize=9)
plt.legend(title="Oncogenicity Class", fontsize=12)

plt.tight_layout()
plt.savefig(f"{save_dir}/gene_categories.png", dpi=300)
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/gene_categories.png'")
print("-"*50)

# ------------------------------------------------------------
# Check null variants distribution 
# for oncogenic and likely neutral classes 
# ------------------------------------------------------------

print("Checking null variant distribution among oncogenic and neutral variants..")

null_var_summary = variants_filtered.groupby(["is_null_variant", "ONCOGENIC"]).size().unstack(fill_value=0)
print("\nSummary:")
print(null_var_summary)
print("-"*50)

# ------------------------------------------------------------
# Plot null variant distribution 
# for oncogenic and likely neutral classes 
# ------------------------------------------------------------

print("Plotting null variant distributions for oncogenic and likely neutral variants..")

null_var_summary.plot(kind='bar', figsize=(8, 5), color=['#88aed1', '#c4314a'], edgecolor="0.1", linewidth=0.3)

plt.title("Prevalence of Null Variants by Oncogenicity Class", fontsize=14)
plt.ylabel("Number of Variants", fontsize=12)
plt.xlabel("Null Variant (True/False)", fontsize=12)
plt.xticks(rotation=45, fontsize=9)
plt.legend(title="Oncogenicity Class", fontsize=12)

plt.tight_layout()
plt.savefig(f"{save_dir}/null_var.png", dpi=300)
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/null_var.png'")
print("-"*50)

# ------------------------------------------------------------
# Check null variants in tsg distribution 
# for oncogenic and likely neutral classes 
# ------------------------------------------------------------

print("Checking null variant in tsg distribution among oncogenic and neutral variants..\n")
null_var_tsg_summary = variants_filtered.groupby(["is_null_var_tsg", "ONCOGENIC"]).size().unstack(fill_value=0)
print("Summary:")
print(null_var_tsg_summary)
print("-"*50)

# ------------------------------------------------------------
# Plot null variants in tsg distribution 
# for oncogenic and likely neutral classes 
# ------------------------------------------------------------

print("Plotting null variant in tsg distribution for oncogenic and likely neutral variants..")

null_var_tsg_summary.plot(kind='bar', figsize=(8, 5), color=['#88aed1', '#c4314a'], edgecolor="0.1", linewidth=0.3)

plt.title("Prevalence of Null Variants in TSGs by Oncogenicity Class", fontsize=14)
plt.ylabel("Number of Variants", fontsize=12)
plt.xlabel("Null Variant in TSG (True/False)", fontsize=12)
plt.xticks(rotation=45, fontsize=9)
plt.legend(title="Oncogenicity Class", fontsize=12)

plt.tight_layout()
plt.savefig(f"{save_dir}/null_var_in_tsg.png", dpi=300)
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/null_var_in_tsg.png'")
print("-"*50)


# ------------------------------------------------------------
# Check gene distribution among positive null variants in TSG 
# ------------------------------------------------------------

# Filter out rows where is_null_var_tsg is True and Oncogenic == "Oncogenic"
pos_null_var_tsg = variants_filtered[
    (variants_filtered["is_null_var_tsg"] == True) & 
    (variants_filtered["ONCOGENIC"] == "Oncogenic")
]

# Which genes do these variants appear in? 
gene_counts = pos_null_var_tsg["Hugo_Symbol"].value_counts()

print("Genes with oncogenic null variants in tumor suppressor genes:")
print(gene_counts)

# Plot results
plt.figure(figsize=(8, 5))
gene_counts.plot(kind="bar", color="#c4314a", edgecolor="0.1", linewidth=0.3)

plt.title("Genes with Oncogenic Null Variants in TSGs", fontsize=14)
plt.ylabel("Number of Variants", fontsize=12)
plt.xlabel("Gene Symbol", fontsize=12)
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig(f"{save_dir}/genes_null_var_in_tsg.png", dpi=300)
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/genes_null_var_in_tsg.png'")
print("-"*50)


print("Oncogene and tumor suppressor gene analysis complete!🥳")