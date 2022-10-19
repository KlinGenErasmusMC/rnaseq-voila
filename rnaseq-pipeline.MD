## RNA-Seq pipeline

Here we briefly describe the steps of our RNA-Seq pipeline and the requirements of the input count data and metadata file to be used in our Juypter Notebook / Voila App for analysis and filtering of RNA-Seq data in clinical genetics diagnostics.

### Mapping

Raw RNA sequencing data (fastq.gz files) will be mapped with Hisat2 against reference genome ucsc.hg19_nohap_masked.

```
hisat2 -x {hisat_index} -1 103805-001-040_R1.fastq.gz -2 103805-001-040_R2.fastq.gz -p 4 --rg-id "RNA" --rg "LB:NEBNextUltraDirectionalRNA" --rg "PL:Illumina" --rg "PU:platform_unit" --rg "SM:103805-001-040" | gatk SortSam --INPUT /dev/stdin --OUTPUT temp("103805-001-040.bam") --TMP_DIR {"tmplocation/"+".namesort.103805-001-040"} --MAX_RECORDS_IN_RAM 500000 --SORT_ORDER queryname
```

Following by MarkDuplicates with gatk:

```
gatk MarkDuplicates --java-options \"-Djava.io.tmpdir={tmp_location} -Xmx3g\" --INPUT 103805-001-040.bam --OUTPUT 103805-001-040_dupmarked.bam --METRICS_FILE 103805-001-040_duplicate_metrics.txt; samtools flagstat 103805-001-040_dupmarked.bam > 103805-001-040_dupmarked.flagstat.txt
```

### Counts

In order to be able to increase the sensitivity of identifying (non-coding) pathogenic variants in RNA-Seq data, we analyse entire gene transcripts and individual exons and introns separately. Therefore we need 3 separate annotation (BED) files with definitions of these. Our pipeline is currently based on the GRCh37.p13 human reference (hg19) and GCF_000001405.25 transcripts. You can find these {BED} files [here](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/reference): 

- genes (GCF_000001405.25_GRCh37.p13_genomic.chr.genes.bed)
- exons (GCF_000001405.25_GRCh37.p13_genomic.chr.transcripts.bed)
- introns (GCF_000001405.25_GRCh37.p13_genomic.chr.introns.bed)

And we use bedtools to obtain counts:

```
bedtools multicov -q 1 -split -bed {BED} -bams 103805-001-040_dupmarked.bam
```

These raw counts and TPM normalised counts (see [normalisation](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/scripts/normalisation_methods.py) ) will be used as input for OUTRIDER analysis. You should end up with count files for your samples as in our [example data](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/countdata)

### OUTRIDER

We use the OUTRIDER (Outlier in RNA-Seq Finder) method for computation of statistics and assess significance of aberrant read counts in individual samples. The outrider script should be performed separately for RNA-Seq (patient) samples from different cell types (e.g. fibroblasts or blood samples) or treatments (with ('CHX') or without ('0') Cycloheximide). We also run it separately for the analysis of genes, exons and introns, with the countfiles from the previous step as input. This workflow for analyzing RNA sequencing data is mainly based on [OUTRIDER.pdf](http://bioconductor.org/packages/release/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf). The workflow runs in R and requires the OUTRIDER package and a number of other R packages to be installed (see R script below). It's recommended to run it on a server / high performance computing (HPC) infrastructure with at least 16 cores and between 64-128 GB memory available. You can run the workflow by adapting the following R script to run in your computing environment with your own data :

([OUTRIDER.R](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/scripts/outrider.R))

```
Rscript outrider.R 'experiment' 'species' 'treatment' 'fragment'
```

The result tables (3 for every cell/sample type/treatment) finally will be used as input tables for the voila application.

## Citation
Brechtmann F, Mertes C, Matusevičiūtė A, Yépez VA, Avsec Ž, Herzog M, Bader DM, Prokisch H, Gagneur J. _OUTRIDER: A Statistical Method for Detecting Aberrantly Expressed Genes in RNA Sequencing Data._ **Am J Hum Genet. 2018 Dec 6;103(6):907-917**; doi: https://doi.org/10.1016/j.ajhg.2018.10.025
