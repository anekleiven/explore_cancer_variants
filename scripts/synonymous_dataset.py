"""
Synonymous dataset
------------------------------------------------------

Script: neutral_genie.py
Author: Ane Kleiven

This script extracts possibly neutral variants from the somatic variant dataset, to enrich the number of likely neutral variants 
in analyses such as MAVEs and germline proximity. 

Filtered on the following: 
  ONCOGENIC != "Oncogenic", "Likely Oncogenic"
  gnomAD_AF >= 0.0001 
  Consequence = "synonymous_variant"

"""

# Import libraries 
import pandas as pd
import numpy as np 

# -----------------------------------------------
# Read variant files 
# -----------------------------------------------

# Load somatic variant file 
print("\nLoading somatic variants..")
somatic_var = pd.read_csv("/home/anekl/git/master/cancer_variants_annotation_pipeline/output/variants_with_maves.tsv",
    sep="\t",
    low_memory=False)
print(f"\nLoaded {len(somatic_var)} somatic variants.\n")


# -----------------------------------------------
# Filter variants
# -----------------------------------------------

# Extract variants using the following filters: 
#   Oncogenicity != "Oncogenic", "Likely Oncogenic" 
#   Consequence = synonymous_variant
#   'in_hotspot' = False 


# Consequences to include
is_synonymous = somatic_var["Consequence"].str.contains("synonymous_variant", na=False)

# Filter out the neutral reference
neutral_ref = somatic_var[
    (is_synonymous) & 
    (somatic_var["In_Hotspot"] == False) & 
    (~somatic_var["ONCOGENIC"].isin(["Oncogenic", "Likely Oncogenic"]))
].copy()

# Make log_dist column
neutral_ref["log_dist"] = np.log10(neutral_ref["Germline_Proximity"] + 1)

print(f"The filtered dataset contains {len(neutral_ref)} synonymous variants.\n") 

print("Preview of the filtered dataset:")
print(neutral_ref.head(), "\n")

# Save dataset as .tsv file 
output_path = "/home/anekl/git/master/cancer_variants_annotation_pipeline/output/synonymous_dataset.tsv"
print(f"\nSaving {len(neutral_ref)} synonymous variants to {output_path}..")
neutral_ref.to_csv(output_path, sep="\t", index=False)


print("\nSynonymous variant data complete!🎉")
