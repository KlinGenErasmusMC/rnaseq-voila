## RNA-Seq pipeline

Here we briefly describe the steps of our RNA-Seq pipeline and the requirements of the input count data and metadata file to be used in our Juypter Notebook / Voila App for analysis and filtering of RNA-Seq data in clinical genetics diagnostics.

### Mapping

Raw RNA sequencing data (fastq.gz files) will be mapped with Hisat2 against reference genome ucsc.hg19_nohap_masked.

`
hisat2 -x {hisat_index} -1 103805-001-040_R1.fastq.gz -2 103805-001-040_R2.fastq.gz -p 4 --rg-id "RNA" --rg "LB:NEBNextUltraDirectionalRNA" --rg "PL:Illumina" --rg "PU:platform_unit" --rg "SM:103805-001-040" | gatk SortSam --INPUT /dev/stdin --OUTPUT temp("103805-001-040.bam") --TMP_DIR {"tmplocation/"+".namesort.103805-001-040"} --MAX_RECORDS_IN_RAM 500000 --SORT_ORDER queryname
`

Following by MarkDuplicates with gatk:

`
gatk MarkDuplicates --java-options \"-Djava.io.tmpdir={tmp_location} -Xmx3g\" --INPUT 103805-001-040.bam --OUTPUT 103805-001-040_dupmarked.bam --METRICS_FILE 103805-001-040_duplicate_metrics.txt; samtools flagstat 103805-001-040_dupmarked.bam > 103805-001-040_dupmarked.flagstat.txt
`

### Counts

The mapped reads are counted by bedtools multicov against separate bedfiles for genes (GCF_000001405.25_GRCh37.p13_genomic.chr.genes.bed), exons (GCF_000001405.25_GRCh37.p13_genomic.chr.transcripts.bed) and introns (GCF_000001405.25_GRCh37.p13_genomic.chr.introns.bed) (see [ref](https://github.com/KlinGenErasmusMC/rnaseq-voila/tree/main/ref).

'bedtools multicov -q 1 -split -bed {BED} -bams 103805-001-040_dupmarked.bam'
