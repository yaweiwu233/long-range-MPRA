## Contents

### Python scripts and sbatch scripts
- `BC_demultiplex_supervised.py` is a python script to count eBCs, dBCs and rBCs from fastq files. It will generate a barcode count table and a barcode distribution plot for QC purpose. 
- `BC_demultiplex.sbatch` is a sbatch file to run BC_demultiplex_supervised.py on the HPC cluster

### Jupyter Notebooks
The jupyter notebooks contains scripts to furter analyze barcode count tables and generate all figures in the paper. 
- `01_GATA_library_sortseq_analysis.ipynb` contains script to analyze long-range MPRA data of GATA1/non-GATA1 library at 0,10kb of HBG promoter.
- `02_GATA_library_features.ipynb` contains code to analyze sequence and epigenetic features of GATA1 and non-GATA1 enhancers in out library. This script uses ChIP-seq bigWig files and plasmid-based MPRA data downloaded from the ENCODE portal.
- `03_enhancer_synergy_analysis.ipynb` contains code to analyze long-range MPRA data of LCR HS combination library at 0,3,10,20,50 and 100kb of HBG promoter.
- `04_compare_across_promoters.ipynb` contains code to compare long-range MPRA data across HBG, HBE and GAPDH promoter.

### Additional files
This folder also contains several reference files required for analysis.
- `sample_ID.txt` contains the names of fastq files.
- `dBCs.txt` contains distance barcode information.
- `pYWD15_eBCs.txt` contains barcode information for the LCR HS combination library.
- `total_library_final_with_BC.txt` contains barcode information and annotations for the GATA1/non-GATA1 enhancer library.
- `total_library_final_with_BC_no_head_closest_P` is a modified version where each library member is assigned to its nearest promoter-like element.
  