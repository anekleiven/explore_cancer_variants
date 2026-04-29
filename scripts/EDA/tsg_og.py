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
import seaborn as sns
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

# Melt data from wide to long format 
plot_data = og_tsg_summary.reset_index().melt(id_vars='gene_category')
plot_data.columns = ['Gene Category', 'Oncogenicity', 'Variant Count']

# Define colors for each class 
palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

# Plot 
plt.figure(figsize=(8,5))
sns.set_style("white")

ax = sns.barplot(
    data=plot_data,
    x='Gene Category',
    y='Variant Count',
    hue='Oncogenicity',
    palette=palette,
    edgecolor="black",
    linewidth=0.8
)

sns.despine()

plt.title("Distribution of Variants Across Gene Categories", fontsize=14, pad=15)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10)
plt.xlabel("Gene Category", fontsize=12, labelpad=10)
plt.legend(title="Oncogenicity Class", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/og_tsg_distribution.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/og_tsg_distribution.png'")

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

# set style
plt.figure(figsize=(8, 5))
sns.set_style("white")
palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

# Melt data from wide to long format 
plot_data_null = null_var_summary.reset_index().melt(id_vars='is_null_variant')
plot_data_null.columns = ['Is_Null', 'Oncogenicity', 'Count']

# Plot
ax = sns.barplot(
    data=plot_data_null,
    x='Is_Null',
    y='Count',
    hue='Oncogenicity',
    palette=palette,
    edgecolor="black",
    linewidth=0.8
)

sns.despine()

plt.title("Prevalence of Null Variants by Oncogenicity Class", fontsize=14, pad=15)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10)
plt.xlabel("Null Variant Status", fontsize=12, labelpad=10)
plt.legend(title="Oncogenicity Class", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/null_var.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/null_var.png'")

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

# set style
plt.figure(figsize=(8, 5))
sns.set_style("white")
palette = {
    "Oncogenic": "#c4314a",
    "Likely Neutral": "#88aed1",
}

# Melt data from wide to long format 
plot_data_null_tsg = null_var_tsg_summary.reset_index().melt(id_vars='is_null_var_tsg')
plot_data_null_tsg.columns = ['Is_Null_TSG', 'Oncogenicity', 'Count']

plot_data_null_tsg['Is_Null_TSG'] = plot_data_null_tsg['Is_Null_TSG'].map({True: 'Null Variant in TSG', False: 'Other Variants'})

# Plot
ax = sns.barplot(
    data=plot_data_null_tsg,
    x='Is_Null_TSG',
    y='Count',
    order=['Null Variant in TSG', 'Other Variants'],
    hue='Oncogenicity',
    palette=palette,
    edgecolor="black",
    linewidth=0.8
)

sns.despine()

plt.title("Prevalence of Null Variants in TSGs by Oncogenicity Class", fontsize=14, pad=15)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10)
plt.xlabel("Variant Type", fontsize=12, labelpad=10)
plt.xticks(rotation=0, fontsize=9)
plt.legend(title="Oncogenicity Class", bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)

plt.tight_layout()
plt.savefig(f"{save_dir}/null_var_in_tsg.png", dpi=300, bbox_inches="tight")
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
sns.set_style("white") 

sns.barplot(gene_counts,
            color="#c4314a", 
            edgecolor="black",
            linewidth=0.8)

sns.despine() 

plt.title("Genes with Oncogenic Null Variants in TSGs", fontsize=14, pad=15)
plt.ylabel("Number of Variants", fontsize=12, labelpad=10)
plt.xlabel("Gene", fontsize=12, labelpad=10)
plt.xticks(rotation=45, fontsize=9)

plt.tight_layout()
plt.savefig(f"{save_dir}/genes_null_var_in_tsg.png", dpi=300, bbox_inches="tight")
plt.show()

print(f"Plotting complete! Plot saved as '{save_dir}/genes_null_var_in_tsg.png'")
print("-"*50)


print("Oncogene and tumor suppressor gene analysis complete!🥳")