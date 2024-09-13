## RNA-Seq Pipeline / Data (input data for Analysis App)

copied & edited from: https://github.com/KlinGenErasmusMC/rnaseq-voila

# rnaseq-voila

Below, we describe how to install our Jupyter Notebook / Voila App for analysis and filtering of RNA-Seq data for clinical genetics diagnostics locally on your laptop / PC or (institute / department) server. First you need to choose a location where to download the data. In this example we will use /home/notebook/ as the destination to install this github repository

## Prerequisites - RNA-Seq Analysis App

To run the Jupyter Notebook / Voila App, you need a working python and Jupyter / Voila installation with the following Python modules (as installed locally @Lonneke):

    - dominate=2.9.1
    - jupyter_core=5.7.2
    - jupyterlab=4.1.5
    - jupyter_server=2.25.4
    - python=3.12.1
    - plotly=5.1.0
    - nbconvert=7.16.2
    - matplotlib=3.8.3
    - pandas=2.2.1
    - ipywidgets=8.1.2
    - scikit-learn=1.4.1.post1
    - wheel=0.41.2
    - voila=0.5.5  

## Installation - using pip

The requirements.txt file should lists all required Python libraries that the notebook depend on, and they can be installed using:

```
pip install -r requirements.txt
```

## Installation - using Anaconda / Miniconda

If you prefer, you can use the environment.yml file to install the whole environment using [Anaconda / Miniconda](https://docs.conda.io/en/latest/miniconda.html) :

```
conda env create --file environment.yml
```

When the conda environment has been successfully created, it can be activated with the following command

```
conda activate rnaseq-voila
```

and used to run the Notebook server and/or Voila app as detailed below. 

## Running the Notebook / Voila App

You can start the Notebook server from the folder where you downloaded the data to (here **/home/notebook/rna-voila-main/** ).  

```
jupyter-notebook --NotebookApp.ip=localhost --no-browser --port=8888 --notebook-dir=/home/notebook/rnaseq-voila-main/
```

Start the app at: http://localhost:8888/edit/rnaseq-filtering-app.ipynb?factory=Voila+Preview

(You can interact with the Python code in the notebook at http://localhost:8888/tree)
 
## Citation
Jordy Dekker, Rachel Schot, Michiel Bongaerts, Walter G. de Valk, Monique van Veghel-Plandsoen, ..., Grazia M.S. Mancini, Tjakko J. van Ham. _RNA-sequencing improves diagnosis for neurodevelopmental disorders by identifying pathogenic non-coding variants and reinterpretation of coding variants._ **medRxiv**; doi: https://doi.org/10.1101/2022.06.05.22275956
