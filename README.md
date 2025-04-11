# ISSAAC-seq preprocessing pipeline v2
ISSAAC-seq is a single-cell multiomic method for simultaneous profiling of chromatin accessibility and gene expression from the same cell.

## Raw Sequencing Data

The raw data (fastq files) are deposited into ArrayExpress under the accession number [E-MTAB-11264](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11264/).

## Analysing your own data

Put the path to the corresponding FASTQ files in `config.json`. To run the pipeline using all available cores, simply run:

```
snakemake --cores
```

## Softwares/Packages prerequisite:

1. __General programs__
  - Python3
  - R v4
2. __Python packages__
  - numpy (v1.19.5)
  - scipy (v1.5.3)
  - pandas (v1.1.5)
  - seaborn (v0.11.0.rc0)
  - scanpy (v1.9.0)
  - scikit-learn (v0.22.1)
3. __R packages__
  - Matrix (v1.3-3)
  - Signac (v1.6.0)
  - Seurat (v4.1.0)
  - monocle (v2.20.0)
  - GenomeInfoDb (v1.30.1)
  - EnsDb.Mmusculus.v79 (v2.99.0)
  - EnsDb.Hsapiens.v86 (v2.99.0)
  - ggplot2 (v3.35)
  - patchwork (v1.1.1)
  - dplyr (v1.0.8)
  - stringr (v1.4.0.9000)
  - mgsub (v1.7.3)
  - stringi (v1.7.6)
4. __Other Bioinformatics Utilities__
  - chromap (0.2.1-r369)
  - MACS2 (v2.2.7.1)
  - STARsolo (v2.7.11b)
  - samtools (v1.9)
  - bedtools (v2.30.0)
  - seqtk (v1.3-r106)
  - faSize, bedClip, calc, addCols from [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)

## Contact

Xi Chen  
chenx9@sustech.edu.cn
