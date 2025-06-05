# ISSAAC-seq preprocessing pipeline v2
ISSAAC-seq is a single-cell multiomic method for simultaneous profiling of chromatin accessibility and gene expression from the same cell. This repository contains an updated version of the preprocessing pipeline compared to the [original publication](https://www.nature.com/articles/s41592-022-01601-4) back in 2022.

## Raw Sequencing Data

The 10x Genomics workflow generates data in the same format as the original publication, and you can find the raw data (fastq files) from ArrayExpress under the accession number [E-MTAB-11264](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11264/). For the plate-based workflow, we now put the well barcode to the `i5` position so that it has the same library format as the 10x Genomics workflow, and the raw data generated with PBMCs from a healthy individual can be found under the accession number [E-MTAB-15131](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-15131). The Bio-Rad workflow has different library format, and the raw data can be found under the same accession number [E-MTAB-15131](https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-15131).

## Analysing your own data

1. The following programmes and packages are needed to run the full pipeline. Make sure they are executable and in your `$PATH` so that they can be invoked directly without having to specify the absolute path.

- [cutadapt](https://cutadapt.readthedocs.io) (v4.5)
- [chromap](https://github.com/haowenz/chromap) (0.2.1-r369)
- [MACS2](https://pypi.org/project/MACS2/) (v2.2.7.1)
- [STARsolo](https://github.com/alexdobin/STAR) (v2.7.11b)
- [samtools](https://www.htslib.org) (v1.9)
- [bedtools](https://bedtools.readthedocs.io/en/latest/) (v2.30.0)
- [seqtk](https://github.com/lh3/seqtk) (v1.3-r106)
- faSize, bedClip, calc, addCols from [UCSC utilities](http://hgdownload.soe.ucsc.edu/admin/exe/)

The pipeline also uses `python3` and the following python packages, all of which can be installed via [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) or the like:

```
scipy v1.11.4
numpy v1.26.0
pandas v2.1.2
matplotlib v3.9.1
```


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

## Contact

Xi Chen  
chenx9@sustech.edu.cn
