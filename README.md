# rnaseq-voila

A Jupyter Notebook / Voila App for analysis and filtering of RNA-Seq data for clinical genetics diagnostics

## Citation
Jordy Dekker, Rachel Schot, Michiel Bongaerts, Walter G. de Valk, Monique van Veghel-Plandsoen, ..., Grazia M.S. Mancini, Tjakko J. van Ham. _RNA-sequencing improves diagnosis for neurodevelopmental disorders by identifying pathogenic non-coding variants and reinterpretation of coding variants._ **medRxiv**; doi: https://doi.org/10.1101/2022.06.05.22275956

## Prerequisites - RNA-Seq Data (input data for Analysis App)

- See RNA-Seq processing Pipeline Instructions

## Prerequisites - RNA-Seq Analysis App

- A working Jupyter Notebook / Voila installation with the following Python modules
    - plotly=5.1.0
    - matplotlib-base
    - pandas
    - dominate=2.3.1
    - ipywidgets
    - jupyterlab>=3.1
    - voila>=0.2.11
    - scikit-learn

One can use the environment.yml file to install the whole environment using Anaconda / Miniconda as described below.

## Installation

First setup :

```
git clone https://github.com/KlinGenErasmusMC/rnaseq-voila.git
```

When not using conda, another path required to setup is the location of the pairix binary. peakHiC uses the 4DN-DCIC tool pairix (see **Prerequisites**) to read HiC data in the pairs format. Below we explain how to configure peakHiC to locate this tool. Please make sure this tool is installed before proceding to the next steps. After installing the pairix tool and the necessary R packages from CRAN and Bioconductor, peakHiC installation should take less than one minute. 

## RNA-Seq Data (input data for Analysis App)

For further analysis and visualization of the example data, please open en read the Tutorial.md
 
