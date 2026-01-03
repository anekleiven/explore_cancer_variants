import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os 

# directory to save plots 
save_dir = "visualize_variants/plots"
os.makedirs(save_dir, exist_ok=True) 

#--------------------------------------------------------------------
# Read variants file 
#--------------------------------------------------------------------

variants = pd.read_csv(
  "annotation_pipeline/output/variants_with_func_sites.tsv", 
  sep="\t", 
  low_memory=False
  )

#--------------------------------------------------------------------
# Extract oncogenic variants
#--------------------------------------------------------------------

oncogenic_variants = variants[variants['ONCOGENIC'] == 'Oncogenic']

oncogenic_genes = ( 
  oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

#--------------------------------------------------------------------
# Extract 'Likely Oncogenic' variants
#--------------------------------------------------------------------

likely_oncogenic_variants = variants[variants['ONCOGENIC'] == 'Likely Oncogenic']

likely_oncogenic_genes = ( 
  likely_oncogenic_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

#--------------------------------------------------------------------
# Extract 'Likely Neutral' variants
#--------------------------------------------------------------------

neutral_variants = variants[variants['ONCOGENIC'] == 'Likely Neutral']

neutral_genes = ( 
  neutral_variants["Hugo_Symbol"]
  .value_counts()
  .reset_index(name="Variant_Count") 
  .rename(columns={"Hugo_Symbol": "Gene"})
)

#--------------------------------------------------------------------
# Visualize top 30 genes for all oncogenicity classes 
#--------------------------------------------------------------------

sns.set_theme(style="whitegrid", context="talk") 

# Function to make consistent plots 
def plot_top_genes(df, title, color, plotname):
    """
    Create a consistent barplot for top genes based on oncogenicity class
    """
    plt.figure(figsize=(9, 5))
    
    # sort count values in descending order 
    df_sorted = df.sort_values("Variant_Count", ascending=False)

    # create barplot
    sns.barplot(
        data=df_sorted,
        x="Gene",
        y="Variant_Count",
        color=color
    )

    # Add value labels above bars
    for i, v in enumerate(df_sorted["Variant_Count"]):
        plt.text(i, v + 0.5, str(v), ha="center", va="bottom", fontsize=7, color="black")

    # Style titles and labels
    plt.title(title, fontweight="bold", fontsize=15, pad=15)
    plt.xlabel("Gene", fontsize=12)
    plt.ylabel("Number of Variants", fontsize=12)
    plt.xticks(rotation=45, ha="right", fontsize=10)
    plt.yticks(fontsize=10)

    # Clean style
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"{save_dir}/{plotname}.png", dpi=300, bbox_inches="tight")
    plt.show()
     

# ============================================================
# Plot each oncogenicity class 
# ============================================================

plot_top_genes(
    oncogenic_genes.head(30),
    "Top Genes by Number of Oncogenic Variants",
    color="#ef6f6c",
    plotname="top_oncogenic" 
)

plot_top_genes(
    likely_oncogenic_genes.head(30),
    "Top Genes by Number of Likely Oncogenic Variants",
    color="#f5ad55",
    plotname="top_likely_oncogenic" 
)

plot_top_genes(
    neutral_genes.head(30),
    "Top Genes by Number of Likely Neutral Variants",
    color="#93bdfb",
    plotname="top_likely_neutral"
)


# ============================================================
# Class distribution in the oncogenic genes 
# (Oncogenic and Likely Neutral)
# ============================================================

top_onco_genes = oncogenic_genes.head(30)["Gene"] 

top_genes_variants = variants[variants["Hugo_Symbol"].isin(top_onco_genes)]

distribution = ( 
    top_genes_variants
    .groupby(["Hugo_Symbol", "ONCOGENIC"])
    .size() 
    .reset_index(name="Count")
)

# extract only oncogenic and likely netrual variants
wanted_classes = ["Oncogenic", "Likely Neutral"]
distribution_filtered = distribution[
    distribution["ONCOGENIC"].isin(wanted_classes)]

sns.set_style(style="whitegrid") 

plt.figure(figsize=(12,6)) 

sns.barplot(
    data=distribution_filtered,
    x="Hugo_Symbol",
    y="Count", 
    hue="ONCOGENIC", 
    palette=["#5b8fdc","#ef6f6c"]
)

plt.title("Distribution of Variant Oncogenicity in Top Oncogenic Genes",
          fontsize=16,
          fontweight="bold")
plt.xlabel("Gene", fontsize=13)
plt.ylabel("Number of Variants", fontsize=13) 
plt.xticks(rotation=45, ha="right", fontsize=13)
plt.yticks(fontsize=13)
plt.legend(title="Oncogenicity Class", fontsize=13, loc="upper left") 

plt.tight_layout()

plt.savefig(f"{save_dir}/distribution_top_onco.png", dpi=300, bbox_inches="tight")

plt.show() 

# ============================================================
# Class distribution in the oncogenic genes (PIVOT)
# (Oncogenic and Likely Neutral) 
# ============================================================

pivot = distribution_filtered.pivot(
    index="Hugo_Symbol",
    columns="ONCOGENIC",
    values="Count"
).fillna(0)

pivot_pct = pivot.div(pivot.sum(axis=1), axis=0) * 100 

pivot_pct.plot(
    kind="bar",
    stacked=True, 
    figsize=(12,6),
    color=["#5b8fdc","#ef6f6c"],
    edgecolor="black",
    linewidth=0.3
)

plt.title("Percentage Distribution of Variant Oncogenicity in Top Oncogenic Genes",
          fontsize=16, fontweight="bold") 
plt.xlabel("Gene", fontsize=13)
plt.ylabel("Percentage of Variants (%)", fontsize=13) 
plt.xticks(rotation=45, ha="right", fontsize=13) 
plt.yticks(fontsize=13)
plt.legend(fontsize=13, bbox_to_anchor=(1.05,1), loc="upper left")

plt.tight_layout()

plt.savefig(f"{save_dir}/percentage_top_onco.png", dpi=300, bbox_inches="tight")

plt.show() 
