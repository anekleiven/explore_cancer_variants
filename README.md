# Data Exploration Cancer Variants 👩🏽‍🔬

Script to explore, visualize and perform statistics on cancer variants and the annotated classification evidence. 
Each script contains further explanations on content and use. 

## Requirements 💻
- Python 3.10+

## Setup Instructions 🔧

1. **Create Virtual Environment:**
`python -m venv .venv`
`. .venv/bin/activate`

2. **Install Python Requirements:**
`pip install -r requirements.txt`


## Script Descriptions

`gnomAD_freq.py` - gnomAD frequency distributions among different oncogenicity classes and statistical analyses. 

`oncogenicity.py` - Oncogenicity distributions

`top_genes.py` -  Top genes per oncogenicity class 

`variants_in_domain.py` - Variant distributions inside/outside protein domains 

`variants_in_func_sites.py` - Variant distributions inside/outside functional domains  

`variants_in_hotspots.py` - Variant distributions inside/outside cancer hotspots  

`variants_germline_proximity.py` - Distance between somatic cancer variants and known pathogenic germline variants 

`variants_with_maves.py` - Variant counts with MAVE scores and MAVE score distributions  


## Recommended Sources 🛜

- AACR Project GENIE: https://www.aacr.org/professionals/research/aacr-project-genie/
- OncoKB: https://www.oncokb.org



