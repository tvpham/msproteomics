# msproteomics sitereport: reporting DIA-MS phosphoproteomics experiments at site level with ease

![sitereport](images/sitereport.png)

This repository contains a Python package to report phosphosites and phosphopeptides from a DIA-MS phosphoproteomics experiment.

## Data

* Skowronek et al. PXD034128 dataset
  * Spectronaut 16 output: data/20220530_113551_Phospho_optimal_dia-PASEF_21min_v16_Report.xls
  * DIA-NN 1.8 output: data/Phospho_EGF_diAID.tsv
* Skowronek et al. PXD034128 data re-processsed
  * Spectronaut 18 output: data/20231106_133459_optimal4_Report.tsv
  * DIA-NN 1.8.2 output: data/report.tsv
* Example fasta file
  * uniprot-reviewed_yes_AND_organism__Homo_sapiens__Human___9606___--.zip

## Installation

```
pip install msproteomics
```

On some system, you might need to install pybind11 first, for example

```
conda create -n p37 python=3.7
conda activate p37
pip install pybind11
pip install msproteomics
```

## Example usages

For Spectronaut export

```
sitereport 20231106_133459_optimal4_Report.tsv -tool sn
```

For DIA-NN export

```
read_diann -o report_msproteomics.tsv -f uniprot-reviewed_yes_AND_organism__Homo_sapiens__Human___9606___--.fasta report.tsv

sitereport report_msproteomics.tsv
```