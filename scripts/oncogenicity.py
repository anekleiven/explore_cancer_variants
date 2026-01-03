import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 

# ------------------------------------------------------------
# Load variant data
# ------------------------------------------------------------
variants = pd.read_csv("annotation_pipeline/output/variants_with_func_sites.tsv", sep="\t", low_memory=False)
print(variants.columns)

# ------------------------------------------------------------
# Group by oncogenicity and sort descending
# ------------------------------------------------------------

oncogenicity_df = (
variants["ONCOGENIC"]
.value_counts()
.reset_index()
)

oncogenicity_df.columns = ["Oncogenicity", "Count"]
oncogenicity_df = oncogenicity_df.sort_values("Count", ascending=False) 

print("\nSummary of Oncogenicity classes:\n")
print(oncogenicity_df, "\n") 

print("\nTable output for the latex report:\n")
print(oncogenicity_df.to_latex(index=False))

# ------------------------------------------------------------
# Plot distribution of Oncogenicity classes 
# ------------------------------------------------------------

palette = {
    "Oncogenic": "#ef6f6c",
    "Likely Oncogenic": "#f5b971",
    "Likely Neutral": "#5b8fdc",
    "Inconclusive": "#bdbab2",
    "Unknown": "#848a8e",
    "Resistance": "#ba7ad4"
}

plt.figure(figsize=(8,4))
sns.barplot(
    data=oncogenicity_df,
    x="Oncogenicity",
    y="Count",
    dodge=False, 
    palette=palette)

plt.yscale("log")
plt.title("Distribution of Oncogenicity Classes", fontsize=14, fontweight="bold")
plt.xlabel("Oncogenicity class", fontsize=12)
plt.ylabel("Number of Variants (log-scale)", fontsize=12)
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()

plt.savefig("visualize_variants/plots/oncogenicity.png", dpi=300, bbox_inches="tight")
plt.show()

