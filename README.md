## RNA-Seq Pipeline / Data (input data for Analysis App)

This Notebook / Voila webbrowser App comes with some demo data to highlight the features of the App and is available [here](https://mybinder.org/v2/gh/KlinGenErasmusMC/rnaseq-clingendiag/HEAD?urlpath=voila%2Frender%2Frnaseq-filtering-app.ipynb) . If you want to use this App for the analysis of your own (patient) RNA-Seq data, you need to first process this data. Instructions on how to process data in such a way that they can be used in the App can be found [here](https://github.com/KlinGenErasmusMC/rnaseq-voila/blob/main/rnaseq-pipeline.MD). 

# rnaseq-voila

Below, we describe how to install our Jupyter Notebook / Voila App for analysis and filtering of RNA-Seq data for clinical genetics diagnostics locally on your laptop / PC or (institute / department) server. First you need to choose a location where to download the data. In this example we will use **/home/notebook/** as the destination to download and extract the zip file:

```
cd /home/notebook/
wget https://github.com/KlinGenErasmusMC/rnaseq-voila/archive/refs/heads/main.zip
unzip main.zip
```

All data from this repository should now be available at **/home/notebook/rnaseq-voila-main/**. 

## Prerequisites - RNA-Seq Analysis App

To run the Jupyter Notebook / Voila App, you need a working python and Jupyter / Voila installation with the following Python modules :

    -   plotly==5.1.0
    -   matplotlib
    -   pandas
    -   dominate==2.3.1
    -   ipywidgets
    -   jupyterlab>=3.1
    -   voila>=0.2.11
    -   scikit-learn

## Installation - using pip

The [requirements.txt](https://github.com/KlinGenErasmusMC/rnaseq-voila/blob/main/requirements.txt) file should lists all required Python libraries that the notebook depend on, and they can be installed using:

```
pip install -r requirements.txt
```

## Installation - using Anaconda / Miniconda

If you prefer, you can use the [environment.yml](https://github.com/KlinGenErasmusMC/rnaseq-voila/blob/main/environment.yml) file to install the whole environment using [Anaconda / Miniconda](https://docs.conda.io/en/latest/miniconda.html) :

```
conda env create --file environment.yml
```

When the conda environment has been successfully created, it can be activated with the following command

```
conda activate rnaseq-voila
```

and used to run the Notebook server and/or Voila app as detailed below. 

## Running the Notebook / Voila App

To start the Notebook server, first find out the IP of the machine that will host the server (here **101.66.22.62**) and choose an available port. To more easily navigate to the notebook / app, you can start the Notebook server from the folder where you downloaded the data to (here **/home/notebook/rna-voila-main/** ).  

```
jupyter-notebook --NotebookApp.ip=101.66.22.62 --no-browser --port=8888 --notebook-dir=/home/notebook/rnaseq-voila-main/
```

You can then interact with the Python code in the notebook at [this](http://101.66.22.62:8888/) URL. Please note that you will need to change this URL to suit your local server configuration. 

Working in the Notebook gives you the freedom to do (additional) data analysis on the data yourself. A more user-friendly Voila web-app version can be started by clicking on the Voila logo :

Or alternatively, in case you do not wish to change / interact with the code, you can instead of running the Notebook server, start the Voila app directly :

```
voila --port=8866 --ExecutePreprocessor.timeout=360 --no-browser --Voila.ip=101.66.22.62 /home/notebook/rna-voila-main/rnaseq-filtering-app.ipynb
```
You can then connect to the App at [this](http://101.66.22.62:8866/) URL. A user manual explaining the functionality of the App can be found [here](https://github.com/KlinGenErasmusMC/rnaseq-voila/blob/main/User_manual_RNA-seq-Voila-App.pdf). 
 
## Citation
Jordy Dekker, Rachel Schot, Michiel Bongaerts, Walter G. de Valk, Monique van Veghel-Plandsoen, ..., Grazia M.S. Mancini, Tjakko J. van Ham. _RNA-sequencing improves diagnosis for neurodevelopmental disorders by identifying pathogenic non-coding variants and reinterpretation of coding variants._ **medRxiv**; doi: https://doi.org/10.1101/2022.06.05.22275956
