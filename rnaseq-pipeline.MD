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

The mapped reads are counted by bedtools multicov against separate bedfiles {BED} (see [ref](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/reference)): 
- genes (GCF_000001405.25_GRCh37.p13_genomic.chr.genes.bed)
- exons (GCF_000001405.25_GRCh37.p13_genomic.chr.transcripts.bed)
- introns (GCF_000001405.25_GRCh37.p13_genomic.chr.introns.bed)

```
bedtools multicov -q 1 -split -bed {BED} -bams 103805-001-040_dupmarked.bam
```

These raw counts and TPM normalised counts (see [normalisation.py](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/scripts/normalisation.py) ) will be used as input for OUTRIDER analysis.

### OUTRIDER

The outrider script will separatly be performed for every species (fibroblast or blood samples), treatment (with ('CHX') or without ('0') Cycloheximidine) and for every fragment (genes, exons, introns).
This workflow for analyzing RNA sequencing data is mainly based on [OUTRIDER.pdf](http://bioconductor.org/packages/release/bioc/vignettes/OUTRIDER/inst/doc/OUTRIDER.pdf). The workflow needs a server with a minimum memory of 128GB.

([OUTRIDER.R](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/scripts/outrider.R))

```
Rscript outrider.R 'experiment' 'species' 'treatment' 'fragment'
```

The resulttables (3 for every species/treatment) finally will be used as inputtables for the voila application.
