"""
# ============================================================
# ANALYSIS: Oncogenicity Distribution
# ============================================================

Script: oncogenicity.py
Author: Ane Kleiven

This script performs explorative analyses regarding oncogenicity distribution
in cancer variant data. 

Script content:
--------------
1. Group data by oncogenicity and sort descending. 
2. Plot distribution of oncogenicity classes 

All plots are saved in:
    plots/

"""
print("-"*50)
print("Oncogenicity Distribution Analysis🤓")
print("-"*50)

# ------------------------------------------------------------
# Import libraries 
# ------------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 

# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------

print("Loading variant data..\n")

variants = pd.read_csv("/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_tsg_og.tsv", sep="\t", low_memory=False)

print(f"Successfully loaded {len(variants):,} somatic cancer variants!")

# ------------------------------------------------------------
# Group by oncogenicity and sort descending
# ------------------------------------------------------------

print("-"*30)
print("GROUP AND SORT VARIANT DATA") 

oncogenicity_df = (
variants["ONCOGENIC"]
.value_counts()
.reset_index()
)

oncogenicity_df.columns = ["Oncogenicity", "Count"]
oncogenicity_df = oncogenicity_df.sort_values("Count", ascending=False) 

print("\nSummary of Oncogenicity classes:")
print(oncogenicity_df, "\n") 

print("Table output for the latex report:")
print(oncogenicity_df.to_latex(index=False))


# ------------------------------------------------------------
# Plot distribution of Oncogenicity classes 
# ------------------------------------------------------------

print("-"*30)
print("Plotting distribution of oncogenicity classes..")

palette = {
    "Oncogenic": "#C4473B",
    "Likely Oncogenic": "#D98C6A",
    "Likely Neutral": "#7e8aa2",
    "Inconclusive": "#f9c74f",
    "Unknown": "#848a8e",
    "Resistance": "#ba7ad4"
}

plt.figure(figsize=(8,5))
sns.barplot(
    data=oncogenicity_df,
    x="Oncogenicity",
    y="Count",
    dodge=False, 
    palette=palette,
    edgecolor="0.1",
    linewidth=0.3)

plt.yscale("log")
plt.title("Distribution of Oncogenicity Classes", fontsize=14, pad=10)
plt.xlabel("Oncogenicity Class", fontsize=12)
plt.ylabel("Number of Variants (log-scaled)", fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=9)
plt.tight_layout()

plt.savefig("plots/oncogenicity.png", dpi=300, bbox_inches="tight")
plt.show()

print("Plotting complete! Plot saved as 'plots/oncogenicity.png'\n")


print("-"*30)
print("Oncogenicity distribution analysis complete!🥳\n")
