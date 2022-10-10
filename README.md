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

## Running the Notebook / Voila App

To start the Notebook server, first find out the IP of the machine that will host the server (here **100.66.2.72**) and choose an available port. To more easily navigate to the notebook / app, you can start the Notebook server from the folder where you downloaded the data to (here **/home/notebook/rna-voila/** ).  

```
jupyter-notebook --NotebookApp.ip=100.66.2.72 --no-browser --port=8888 --notebook-dir=/home/notebook/rna-voila/ 
```

You can then interact with the Python code in the notebook at this URL :

Working in the Notebook gives you the freedom to do (additional) data analysis on the data yourself. A more user-friendly Voila web-app version can be started by clicking on the Voila logo :

Or alternatively, in case you do not wish to change / interact with the code, you can instead of running the Notebook server, start the Voila app directly :

```
voila --port=8866 --ExecutePreprocessor.timeout=360 --no-browser --Voila.ip=100.66.2.72 /home/notebook/rna-voila/rnaseq-filtering-app.ipynb
```
You can then connect to the App at the this URL [Voila App](http://100.66.2.72:8866/)]

## RNA-Seq Data (input data for Analysis App)

For further analysis and visualization of the example data, please open en read the Tutorial.md
 
## Citation
Jordy Dekker, Rachel Schot, Michiel Bongaerts, Walter G. de Valk, Monique van Veghel-Plandsoen, ..., Grazia M.S. Mancini, Tjakko J. van Ham. _RNA-sequencing improves diagnosis for neurodevelopmental disorders by identifying pathogenic non-coding variants and reinterpretation of coding variants._ **medRxiv**; doi: https://doi.org/10.1101/2022.06.05.22275956
