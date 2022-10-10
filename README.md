# rnaseq-voila

A Jupyter Notebook / Voila App for analysis and filtering of RNA-Seq data for clinical genetics diagnostics

First get this repo :

```
git clone https://github.com/KlinGenErasmusMC/rnaseq-voila.git
```

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

One can use the environment.yml file to install the whole environment using Anaconda / Miniconda as described below. Alternatively you can install the required Python libraries using pip within an existing Python installation / environment 

## Installation - using pip


The `requirements.txt` file should list all Python libraries that your notebooks
depend on, and they will be installed using:

```
pip install -r requirements.txt
```


## RNA-Seq Data (input data for Analysis App)

For further analysis and visualization of the example data, please open en read the Tutorial.md
 
## Citation
Jordy Dekker, Rachel Schot, Michiel Bongaerts, Walter G. de Valk, Monique van Veghel-Plandsoen, ..., Grazia M.S. Mancini, Tjakko J. van Ham. _RNA-sequencing improves diagnosis for neurodevelopmental disorders by identifying pathogenic non-coding variants and reinterpretation of coding variants._ **medRxiv**; doi: https://doi.org/10.1101/2022.06.05.22275956
