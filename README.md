# Multi-Omics Data Analysis Project

An end-to-end analysis of microbiome and metabolomics data from a mouse dietary intervention study investigating whether different stereoisomers and mixtures of **2,3-butanediol** modify the biological effects of a high-fat diet.

This project was completed as part of **BIOC0023: Specialist Research Project in Multiomics & Data Science at UCL**. The analysis combines 16S rRNA sequencing, metabolite profiling, statistical testing, and data visualisation in Python and R.

## Study design

Male C57BL/6J mice were assigned to a normal-chow control or a high-fat diet. High-fat-diet groups received one of four 2,3-butanediol interventions:

| Code | Experimental group |
| --- | --- |
| `NC` | Normal-chow control |
| `HFD` | High-fat diet |
| `HFD_M` / `HFDM` | HFD + meso-2,3-butanediol |
| `HFD_R` / `HFDR` | HFD + R-2,3-butanediol |
| `HFD_S` / `HFDS` | HFD + S-2,3-butanediol |
| `HFD_MRS` / `HFDMRS` | HFD + mixed 2,3-butanediol stereoisomers |

The two notebooks examine complementary molecular layers:

- **Microbiome:** processing and analysis of 16S rRNA amplicon-sequencing data
- **Metabolomics:** differential and multivariate analysis of liver metabolite profiles

## Repository contents

| File | Description |
| --- | --- |
| `Multi_Omics_Data_Analysis_Project_(microbiome)(1).ipynb` | 16S read processing, denoising, taxonomic classification, diversity analysis, abundance testing, and visualisation |
| `Multi_Omics_Data_Analysis_Project_(metabolites)(1).ipynb` | Metabolite-level hypothesis testing, multiple-testing correction, ordination, PERMANOVA, and selected-metabolite plots |

## Analysis workflow

### 1. Microbiome analysis

The microbiome notebook implements the following workflow:

1. Imports paired-end FASTQ files and removes primers.
2. Merges reads and examines expected-error profiles.
3. Filters reads at a maximum expected error of 1 and truncates them to 350 bp.
4. Dereplicates sequences at both sample and pooled levels.
5. Denoises reads with the UNOISE algorithm.
6. Detects and removes chimeric sequences.
7. Maps reads back to zero-radius OTU centroids at 97% identity.
8. Assigns taxonomy with USEARCH SINTAX and the RDP 16S reference database at a confidence threshold of 0.8.
9. Rarefies the feature table to 23,957 reads per sample.
10. Calculates Shannon alpha diversity and Jaccard beta diversity.
11. Tests overall and pairwise group differences using Mann–Whitney U tests, PERMANOVA, and pairwise PERMANOVA with false-discovery-rate correction.
12. Summarises phylum- and genus-level composition, examines the Firmicutes-to-Bacteroidetes ratio, and performs differential-abundance testing.

### 2. Metabolomics analysis

The metabolomics notebook:

1. Reshapes the metabolite table into tidy long format and joins sample metadata.
2. Compares `NC` and `HFD` abundance for each metabolite using independent-samples t-tests.
3. Calculates fold changes and Cohen's *d* effect sizes.
4. Controls the false discovery rate using the Benjamini–Hochberg procedure.
5. Visualises differential metabolites with a volcano plot.
6. Calculates Bray–Curtis dissimilarities across samples.
7. Uses principal coordinate analysis (PCoA) to visualise global metabolic variation.
8. Tests overall differences among the six experimental groups with PERMANOVA.
9. Produces box-and-strip plots for selected metabolites, including indolepropionate, beta-sitosterol, N-delta-acetylornithine, and 1-myristoyl-2-arachidonoyl-GPC.

## Selected results

- Microbiome composition differed significantly across the six experimental groups (Jaccard PERMANOVA: pseudo-*F* = 4.99, *p* = 0.001; 56 samples).
- Shannon diversity differed between the normal-chow and high-fat-diet groups (Mann–Whitney U: Holm-adjusted *p* = 0.033).
- Significant phylum-level differences were detected for Actinobacteria, Bacteroidetes, Firmicutes, Tenericutes, and Verrucomicrobia, among others.
- Global liver metabolite profiles also differed across groups (Bray–Curtis PERMANOVA: pseudo-*F* = 5.47, *p* = 0.001; 36 samples).
- The metabolomics workflow tested 693 metabolites and reported effect sizes, fold changes, raw *p*-values, and FDR-adjusted *p*-values.

These results describe the outputs saved in the supplied notebooks. They should be interpreted alongside the full experimental report and its biological context.

## Software and dependencies

The notebooks use Python 3, shell commands, and a short R section. Principal dependencies include:

```text
Python: pandas, numpy, scipy, matplotlib, seaborn, scikit-bio,
        statsmodels, pingouin, qiime2pandas, rpy2

Command line: cutadapt, VSEARCH, USEARCH 12, gsutil

R: vegan, tidyverse, devtools, pairwiseAdonis
```

`USEARCH 12` is separately distributed software and is not included here. The `qiime2pandas` package is installed directly from its GitHub repository in the microbiome notebook.

## Running the notebooks

The current notebooks were developed in **Google Colab** and contain Colab-specific paths such as `/content/` and Google Drive copy commands. To reproduce the analysis:

1. Open the notebooks in Google Colab.
2. Mount Google Drive when prompted.
3. Provide the required input files and update the Drive paths to match their locations.
4. Install the Python, command-line, and R dependencies listed above.
5. Run the cells sequentially, checking that feature-table sample IDs match the metadata indices before diversity or statistical analysis.

Required inputs include:

- Paired-end 16S FASTQ files
- Primer-removal/merging script
- Sample metadata (`metadata.tsv` and `metadataM.csv`)
- RDP 16S reference database
- `all_metabolites.csv`
- USEARCH 12 executable

The raw data, metadata, reference database, helper script, and USEARCH executable are not included with the notebooks. Paths and sample-column selections may need to be generalised before running the workflow in another environment.

## Outputs

The analysis generates:

- Rarefaction curves and a rarefied feature table
- Shannon-diversity boxplots
- PCoA plots of microbiome and metabolome profiles
- PERMANOVA and pairwise-comparison results
- Phylum- and genus-level abundance plots
- Firmicutes-to-Bacteroidetes ratio plots
- Microbiome and metabolomics volcano plots
- A table of metabolite-level statistical results
- Boxplots for selected metabolites

Figures are exported in PNG and PDF formats by the relevant notebook cells.

## Notes

- The notebooks preserve the original exploratory analysis, including installation cells and alternative plotting approaches.
- Some filenames and figure names are reused across cells and may overwrite earlier outputs.
- The microbiome notebook includes both Python and R, so a working `rpy2`/R configuration is required for pairwise PERMANOVA.
- Statistical results depend on the supplied data, filtering thresholds, random seeds, and software versions.

## Author

Jocelyn Cheung  
UCL BSc Biochemistry — BIOC0023 Specialist Research Project in Multiomics & Data Science

