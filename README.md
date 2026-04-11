# Data Exploration Cancer Variants 👩🏽‍🔬

Scripts to explore, visualize and perform statistics on somatic cancer variants and evidence for oncogenicity classification. 
Each script contains further explanations on purpose and outputs.

The classification evidence includes: 
- Population frequencies from gnomAD
- Cancer mutation hotspots from cancerhotspots.org
- Protein domain annotations from Pfam
- Functional site annotations from UniProt
- Distance to known pathogenic germline variants from ClinVar 
- Functional score assays from MaveDB 
- Null variants in Tumor Suppressor genes 

## Requirements 💻
- Python 3.10+

## Setup Instructions 🔧

1. **Clone the repository:**
```bash
git clone https://github.com/anekleiven/explore_cancer_variants.git
``` 

2. **Create Virtual Environment:**
`python -m venv .venv`
`. .venv/bin/activate`

3. **Install Python Requirements:**
`pip install -r requirements.txt`


## Scripts included in this repository 

### Exploration of Feature Prevalence in Oncogenic vs. Neutral Variants 

`oncogenicity.py`
`top_genes.py`
`gnomAD_freq.py`
`cancerhotspots.py`
`protein_domains.py`
`functional_sites.py`
`germline_proximity.py`
`neutral_dataset_clinvar.py`
`maves.py`
`tsg_og.py`

### Statistics 

`stats_function.py`
`stats_all_variants.py`
`stats_top_genes.py`
`stats_protein_domains.py`
`stats_func_sites.py`
`stats_germline_proximity.py`
`stats_maves.py` 


## Recommended Sources 🛜

- AACR Project GENIE: https://www.aacr.org/professionals/research/aacr-project-genie/
- OncoKB: https://www.oncokb.org



