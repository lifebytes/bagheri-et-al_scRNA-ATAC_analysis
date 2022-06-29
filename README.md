# Single-cell transcriptome + ATAC analysis for custom barcode infected cells - Breast cancer data
Analysis tutorial for the data published under the manuscript title `Therapeutic induction of mesenchymal-epithelial transition via epigenetic reprogramming curtails metastatic progression and sensitizes breast cancers to treatment`, from the Pattabiraman lab. The paper is available on the [Manucsript](link).

## Table of Contents
* [I. scRNA+ATAC analysis](#I)
   * [1. Barcode whitelist preparation](#I.1)
   * [2. Setting up libraries file](#I.2)
   * [3. Running cellranger](#I.3)
   * [4. Standard Seurat/Signac analysis](#I.4)
* [II. Supplementary Scripts](#II)
   * [1. Barcode classifier](#II.1)
       * [1.1 Revise MULTIseqDemux assignment](#II.1.1)
       * [1.2 Denovo lentiviral barcode assignment](#II.1.2)
   * [2. Muller plot](#II.2)
   * [3. GEO](#II.3)

## I. scRNA+ATAC analysis <a name="I"></a>
As this experimental setup involves 'Custom' barcode from lentiviral barcoding, cellranger feature barcoding based methodlogy is used to analyse the lentiviral barcode assignment later with our custom scripts that uses the raw feature barcoding counts and assigns 'Custom' barcodes to each cell.

### 1. Barcode whitelist preparation <a name="I.1"></a>
We have prepared a whitelist of barcode by combining BC14 and BC30 (BC14sequence+BC30sequence) barcodes. In custom barcode only sequencing, barcodes are expected in the following format, 

`Custom barcode	=	19bp constant + BC14 sequence + TGGT + BC30 sequence` 

We have prepared the feature reference file \[as per 10X format\] to be used for feature barcode demultiplexing in downstream for custom barcodes FASTQ files. You can download the reference feature file used for this analysis from [GEO](here)

### 2. Setting up libraries file <a name="I.2"></a>
Setup the libraries file as in 10X cellranger preferred format as follows,
```
fastqs,sample,library_type
<fastq-path>,<sample>,Gene Expression
<fastq-path>,<sample>,Custom
```
```
fastqs,sample,library_type
<fastq-path>,<sample>,Chromatin Accessibility
<fastq-path>,<sample>,Gene Expression
```
You can downlaod the reference libraries file from [GEO](here)

### 3. Running cellranger <a name="I.3"></a>
For run1, where GE and ATAC libraries were generated separately, GE libraries were processed using cellranger with FB option for 'Custom' barcode counting. For ATAC libraries, cellranger-ATAC was used and CITESeq-count was used to generate 'Custom' barcode counting. Steps are given below,
```
<cellranger-3.1.0-path>/cellranger count --id <experiment_name> --transcriptome <refdata-cellranger-mm10-3.0.0> --expect-cells <#cells> --libraries <library_file_from_#2> --feature-ref <feature_reference_from_#1>

<cellranger-atac-1.2.0-path>/cellranger-atac count --id <experiment_name> --reference <refdata-cellranger-atac-mm10-1.2.0> --fastqs <fastqs_path> --sample <sample-names>
CITE-seq-Count -R1 <read1-runs-merged> -R2 <read2-runs-merged> -t <feature_reference_from_#1-customTrimmed>* -cbf 1 -cbl 16 -umif 17 -umil 26 -cells <#cells> -o <output_folder> --max-error 4 -trim 19
```
*\*As CITEseq-count was slow to handle thousands of 'Custom' barcodes, we have used only the barocdes that are identified in GE libraries to speed up the process*.

In case of multiome analysis, feature barcoding is not allowed in cellranger-arc pipeline. Hence, for standard analysis cellranger-arc pipeline was followed, but for feature barcoding GEX libraries were processed using cellranger count with feature barcoding of 'Custom' barcodes. And the 'Custom' barcodes data were brought into analysis within Seurat environment for lentiviral barcode assignment. As ATAC and GEX data are generated simultaneously, same lentiviral barcode identified in GEX can be assigned to ATAC cells as well. Steps are given below,
```
<cellranger-arc-1.0.0-path>/cellranger-arc count --id <experiment_name> --reference <refdata-cellranger-arc-mm10-2020-A> --libraries <library_file_from_#2>
<cellranger-4.0.0-path>/cellranger count --id <experiment_name> --transcriptome <refdata-cellranger-mm10-3.0.0> --expect-cells <#cells> --libraries <library_file_from_#2> --feature-ref <feature_reference_from_#1>
```

### 4. Standard Seurat/Signac analysis <a name="I.4"></a>
Standard Seurat and Signac analysis were carried out for scRNA and scATAC datasets. Multiple runs data were merged with batch effect corrections. Single cell RNA and ATAC were integrated using Seurat pipeline.

<span style="color:red"> ***Tools and reference versions are specific to the publication analysis. One can use difference versions subjective to parameters change.*** </span>

## II. Supplementary scripts <a name="II"></a>

### 1. Barcode classifier <a name="II.1"></a>
`MULTIseqDemux` function from Seurat was used to assign lentiviral barcodes for each cell. As this and other similar functions are meeant for sample level barcode [HTO] demultiplexing, where the background cutoff is estimated for each barcode across cells. Lentiviral barcoding will have variability across cells and so the above mentioned methods can overestimate the background/signal cutoffs. Hence we have developed a custom barcode classifier where background/signal is estimated for each cell by comparing the expression level of all the barcodes. Also, for other users who want to bypass other tools and run only these custom scripts for barcode assignment, we have separate scripts.


1. **Revise MULTIseqDemux assignment:**  `revise_multiseqdemux_lentiviral_barcode_non-singlets.R`: Takes MULTISeqdemux or similar tools annotaiton of 'singlet', 'doublet', 'negative' as metadata to filter only non-singlets and run the script to re-classify the cells only for non-singlets.
2. **Denovo lentiviral barcode assignment:** `lentiviral_barcode_assigner_scRNA.R`: Independent of MULTIseqdemux like tools, can just run on simple barcode matrix file

#### 1.1 Revise MULTIseqDemux assignment <a name="II.1.1"></a>
**Usage**
```
Rscript scripts/revise_multiseqdemux_lentiviral_barcode_non-singlets.R <matrix raw counts> <meta data> <output directory path>
```

**Arguments**
```
Script name: "revise_multiseqdemux_lentiviral_barcode_non-singlets.R"

Total no.of arguments: 3

Arguments description:
args[1]:raw_matrix
Description: Barcode matrix raw counts for the sample of interest: provide the file path (csv format)
View:
          cell1 cell2 cell3 cell4 cell5
Xkr4        0     0     2     0     3
Gm1992      4     0     0     4     0
Gm37381     0     0     2     0     0
Rp1         0     0     0     0     0
Sox17       0     1     0     1     0
Gm37323     0     0     0     1     0

args[2]:meta
Description: sample meta data: provide meta data file path (csv format) (make sure the rownames of meta data and column names of raw matrix are matching and "Barcode_class" column should be there in the meta data)
Additional note: "Barcode_class" column denotes the prior classification of barcodes as singlets, doublets and negative.

        Barcode_class Run_sample  Run_tech
cell1       Doublet     Run2     Run2-scRNA
cell2      Negative     Run2     Run2-scRNA
cell3       Singlet     Run2     Run2-scRNA
cell4       Doublet     Run1     Run1-ATAC
cell5       Singlet     Run1     Run1-ATAC

args[3]:out
Description: output directory: output file path to save the files and plots (without prefix name)
```
#### 1.2 Denovo lentiviral barcode assignment <a name="II.1.2"></a>
**Usage**
```
Rscript scripts/lentiviral_barcode_assigner_scRNA.R <matrix raw counts> <output directory path>
```

**Arguments**
```
Script name: "lentiviral_barcode_assigner_scRNA.R"

Total no.of arguments: 2

Arguments description:
args[1]:raw_matrix
Description: Barcode matrix raw counts for the sample of interest: provide the file path (csv format)
View:
          cell1 cell2 cell3 cell4 cell5
Xkr4        0     0     2     0     3
Gm1992      4     0     0     4     0
Gm37381     0     0     2     0     0
Rp1         0     0     0     0     0
Sox17       0     1     0     1     0
Gm37323     0     0     0     1     0

args[2]:out
Description: output directory: output file path to save the files and plots (without prefix name)
```

**sessionInfo()**
- R package (4.0.3)
- plyr (1.8.6)
- reshape2 (1.4.4)
- ggplot2 (3.3.2)
- ggpubr (0.4.0)

### 2. Muller plot <a name="II.2"></a>
This custom R script can be run to generate Muller plots as used in the publication Figure 2 I&J.

**Usage**
```
Rscript scripts/muller_plot_barcode_evolution.R <metadata> <output_file_plot_prefix>
```

**Arguments**
```
Script name: "muller_plot_barcode_evolution.R"

Total no.of arguments: 

Arguments description:
args[1]: Seurat metadata as csv
Description: MEtadata with pseudotime and other sample information used in the publication data. Refer to columsn and column names requird for the script

args[2]: output_plot_prefix
Description: Prefix for the output muller plot file. Prefix will be suffixed with 'muller_plot_' and date of plot generation.
```

**sesssionInfo()**
- ggplot2
- reshapre2
- cowplot
- MullerPlot

### 3. GEO <a name="II.3"></a>
10X ouput, UMAP projections and metadata are uploaded to GEO, and can be accessed with the ID: `GEO`
Raw data \[FASTQ\] can be downlaoded from SRA directed via the same GEO link
