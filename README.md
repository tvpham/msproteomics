# msproteomics sitereport: reporting DIA-MS phosphoproteomics experiments at site level with ease

![sitereport](images/sitereport.png)

This repository contains a Python package to report phosphosites and phosphopeptides from a DIA-MS phosphoproteomics experiment.

**Citation**

Pham TV, Henneman AA, Truong NX, Jimenez CR. msproteomics sitereport: reporting DIA-MS phosphoproteomics experiments at site level with ease, _Bioinformatics_ 2024 Jul 1;40(7):btae432.
[https://doi.org/10.1093/bioinformatics/btae432](https://doi.org/10.1093/bioinformatics/btae432)


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
conda create -n msproteomics
conda activate msproteomics
pip install pybind11
pip install msproteomics
```

On a Windows Cygwin system with python39, g++ (GCC) 12.4.0, and matplotlib 3.5.1 pre-installed, the following works

```
python -m pip install setuptools
python -m pip install pybind11
python -m pip install numpy==1.22.4 
python -m pip install pandas==2.2.3
python -m pip install msproteomics
```

## Example usages

For Spectronaut export

```
sitereport 20231106_133459_optimal4_Report.tsv -tool sn
```

For DIA-NN export

```
read_diann -o report_msproteomics.tsv -E Fragment.Quant.Raw -f uniprot-reviewed_yes_AND_organism__Homo_sapiens__Human___9606___--.fasta report.tsv

sitereport report_msproteomics.tsv
```