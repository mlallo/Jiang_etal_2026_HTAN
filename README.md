# ZFP36L2 Orchestrates Stress-adaptive Plasticity During Intestinal Regeneration and Metastasis

This repository provides the computational workflow used for the Human Tumor Atlas Network (HTAN) single-cell RNA sequencing analyses accompanying *ZFP36L2 orchestrates stress-adaptive plasticity during intestinal regeneration and metastasis*, published in *Nature*. These analyses characterize ZFP36L2-dependent transcriptional programs associated with stress-induced dedifferentiation, epithelial plasticity, and intestinal stem cell state acquisition during intestinal regeneration and colorectal cancer metastasis.

## Repository Structure 
```bash
Jiang_etal_2026_HTAN/
├── notebooks_scripts/
│   ├── 1_AnnData_Preprocessing.ipynb
│   ├── 2a_ZFP36L2_Autocorrelation_Normal.py
│   ├── 2b_ZFP36L2_Autocorrelation_Tumor.py
│   └── 3_HotSpot_ZFP36L2_Autocorrelation_Plots.ipynb
├── input/
│   └─ pathways_for_gsea.gmt
├── scrna_py.yml
└── README.md
```
## Installation and Data Download
To run the analyses, clone and install the repository as follows:

Clone the repository
```bash
git clone https://github.com/joechanlab/Jiang_etal_2026_HTAN.git

cd Jiang_etal_2026
```
Create the Conda environment
```bash
conda env create -f scrna_py.yml

conda activate scrna_py.yml
```
Download the HTAN AnnData object as described in [*Moorman et al.* (2025, *Nature*)](https://www.nature.com/articles/s41586-024-08150-0)
```
https://dp-lab-data-public.s3.us-east-1.amazonaws.com/progressive-plasticity-crc-metastasis/h5ads/Epithelial.h5ad
```

## Directory Structure To Process HTAN Figures
<pre><code>Jiang_etal_2026_HTAN/
├── notebooks_scripts/
│   ├── 1_AnnData_Preprocessing.ipynb
│   ├── 2a_ZFP36L2_Autocorrelation_Normal.py
│   ├── 2b_ZFP36L2_Autocorrelation_Tumor.py
│   └── 3_HotSpot_ZFP36L2_Autocorrelation_Plots.ipynb
├── input/
│   ├── pathways_for_gsea.gmt
│   └── <b>Epithelial.h5ad</b> (download using the above instructions) 
└── <b>output/</b> (created automatically by sequentially running the analysis workflow in notebooks_scripts/) 
    ├── processed_adata/
    │   └──processed_adata.h5ad
    └── figures/
         ├── figures_main/ 
         │  ├── figure_1b
         │  ├── figure_1g
         │  ├── figure_1h
         │  ├── figure_1i
         │  └── figure_1j
         └── figures_extended_data/ 
             ├── figure_1a
             ├── figure_1b
             ├── figure_1c
             └── figure_1g
</code></pre>

## Reproducing the Analyses
The entire workflow can be reproduced by modifying a single variable in each script.

Open the desired notebook or script and **specify the repository location**:
```
main_directory = Path("/path/to/Jiang_etal_2026_HTAN")
```

All subsequent output directories with respective final figures and corresponding dataframes are generated automatically from this location.
