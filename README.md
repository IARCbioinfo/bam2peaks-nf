# bamm2peaks -nf

## Nextflow pipeline for peaks calling with MACS

![Workflow representation](pipeline.png?raw=true "Peaks calling with MACS")

## Description

Nextflow pipeline designed for peak calling using MACS and IDR, coupled with QC generation using deeptools. The saturation option generates peaks by successively considering increasing percentages of the total reads, repeating the operation multiple times within the range of 0.05 to 0.95. 

## Dependencies

1. Nextflow : for common installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.

### MACS
2. [*MACS2*](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/06_peak_calling_macs.html) or [*MACS3*](https://github.com/macs3-project/MACS)

### IDR (Irreproducible Discovery Rate)
3. [*IDR*](https://github.com/nboley/idr)

### Deeptools
4. [*Deeptools*](https://github.com/deeptools/deepTools)

**A conda receipe, and docker and singularity containers are available with all the tools needed to run the pipeline (see "Usage")**

## Input 
 | Type      | Description   |
 |-----------|---------------|
 | --input_file    | input tabulation-separated values file with columns sample (sample name), tag (short name for figures), bam (bam file path) and group (group) |

## Parameters

* #### Mandatory

| Name | Example value | Description |
|-----------|--------------|-------------| 
|--ref | hg38 | Reference fasta file hg19, hg38 or mm10' |

* #### Optional

| Name | Default value | Description |
|-----------|--------------|-------------| 
|--mode   | atac | There is two mode : atac or chip, chip require "input" sample(s)|
|--output_folder   | bam2peaks | Output folder name |
|--cpu | 16 | number of CPUs |
|--mem | 16 | memory|
|--extsize | 150 | MACS extsize : extendsize of peaks to to fix-sized fragments. |


* #### Flags

Flags are special parameters without value.

| Name  | Description |
|-----------|-------------| 
|--help | print usage and optional parameters |
|--broad     | Compute broadpeaks instead of narrowpeaks |
|--ignoreDuplicates  | Ignore duplicates reads |
|--saturation        | Run saturation process |


## Usage
To run the pipeline for ATAC, one can type:
  
```bash
nextflow run iarcbioinfo/bam2peaks-nf -r latest -profile singularity --input_file input.tsv --ref hg38 --output_folder output --ignoreDuplicates
```

To run the pipeline without singularity just remove "-profile singularity". Alternatively, one can run the pipeline using a docker container (-profile docker) the conda receipe containing all required dependencies (-profile conda).

### Chip-seq mode
To use the pipeline for Chip-seq, add the --chip flag :

```bash
nextflow run iarcbioinfo/bam2peaks-nf -r latest -profile singularity --input_file input.tsv --ref hg38 --output_folder output --chip --broad --extsize 320
```

## Output 
  | Type      | Description     |
  |-----------|---------------|
  | bw/     | Outputs of bamCoverage in bigWig format  |
  | Counts/ | With --saturation return the number of reads for each subsets |
  | Peaks   | Peaks computed by MACS |
  | Peaks_intersect | Peaks intersections computed by idr |
  | QCs | deeptools graphics |
  | Saturation_peaks | With --saturation, all peaks files for each subsets |
  

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | Vincent Cahais | CahaisV@iarc.who.int | Developer to contact for support |
  | Claire Renard  | Renardc@iarc.who.int | Developer |
  
