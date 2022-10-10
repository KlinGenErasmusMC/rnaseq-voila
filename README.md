# rnaseq-voila

Here, we describe how to install our Jupyter Notebook / Voila App for analysis and filtering of RNA-Seq data for clinical genetics diagnostics locally on your laptop / PC or (institute / department) server. First you need to choose a location where to download the data. In this example we will use **/home/notebook/rnaseq-voila** as the destination :

```
cd /home/notebook/
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

## Installation - using pip

The `requirements.txt` file should list all Python libraries that your notebooks
depend on, and they will be installed using:

```
pip install -r requirements.txt
```

## Installation - using Anaconda / Miniconda

If you prefer, you can use the [environment.yml](https://github.com/KlinGenErasmusMC/rnaseq-voila/blob/main/environment.yml) file to install the whole environment using [Anaconda / Miniconda]() :

```
conda env create --file environment.yml
```

## Running the Notebook / Voila App

To start the Notebook server, first find out the IP of the machine that will host the server (here **100.66.2.72**) and choose an available port. To more easily navigate to the notebook / app, you can start the Notebook server from the folder where you downloaded the data to (here **/home/notebook/rna-voila/** ).  

```
jupyter-notebook --NotebookApp.ip=100.66.2.72 --no-browser --port=8888 --notebook-dir=/home/notebook/rna-voila/ 
```

You can then interact with the Python code in the notebook at [this](http://100.66.2.72:8888/) URL.

Working in the Notebook gives you the freedom to do (additional) data analysis on the data yourself. A more user-friendly Voila web-app version can be started by clicking on the Voila logo :

Or alternatively, in case you do not wish to change / interact with the code, you can instead of running the Notebook server, start the Voila app directly :

```
voila --port=8866 --ExecutePreprocessor.timeout=360 --no-browser --Voila.ip=100.66.2.72 /home/notebook/rna-voila/rnaseq-filtering-app.ipynb
```
You can then connect to the App at [this](http://100.66.2.72:8866/) URL.

## RNA-Seq Pipeline / Data (input data for Analysis App)

This Notebook / Voila App comes with some demo data to highlight the features of the App. If you want to use this App for the analysis of your own (patient) RNA-Seq data, you need to first process this data. Instructions on how to process data in such a way that they can be used in the App can be found [here](https://github.com/KlinGenErasmusMC/rnaseq-voila/blob/main/rnaseq-pipeline.MD)
 
## Citation
Jordy Dekker, Rachel Schot, Michiel Bongaerts, Walter G. de Valk, Monique van Veghel-Plandsoen, ..., Grazia M.S. Mancini, Tjakko J. van Ham. _RNA-sequencing improves diagnosis for neurodevelopmental disorders by identifying pathogenic non-coding variants and reinterpretation of coding variants._ **medRxiv**; doi: https://doi.org/10.1101/2022.06.05.22275956
