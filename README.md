## Economic and Ecological Multifunctionality Analysis – Jena Experiment

This repository contains the full analysis pipeline for the manuscript examining economic and ecological multifunctionality in the Jena Experiment, including the dominance experiment sensitivity analysis and insurance value calculations.

All analyses are implemented in R and fully reproducible via a project-specific package environment managed with renv.

_______________________________________________________________________________________________________________________________
## Repository structure:
.
├── 1_Data/                # Input data
├── 2_Code/                # R Markdown analysis scripts
├── 3_Graphs/              # Output figures (generated)
├── functions/             # Custom plotting/theme functions
├── run_all.R              # Entry point to run full pipeline
├── renv.lock              # Reproducible package environment
└── Data Analysis.Rproj
_______________________________________________________________________________________________________________________________
## Reproducibility

### Step 1 – Clone the repository

git clone YOUR_GITHUB_REPO_URL
cd yourrepo

### Step 2 – Restore the R package environment

Open R in the project directory and run:

install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore(prompt = FALSE)

--> This installs the exact package versions used to generate the results.

### Step 3 – Run the full analysis pipeline

From the project root:

Rscript run_all.R

--> Activates the project-specific renv environment
--> Sets a fixed random seed (set.seed(1234))
--> Renders all analysis scripts in 2_Code/
--> Writes all figures to 3_Graphs/

## Bayesian Model Reproducibility

All Bayesian models were run using fixed random seeds to ensure reproducibility.
Minor numerical differences may occur across platforms due to Monte Carlo sampling and system-level differences in Stan compilation, but results should be statistically equivalent.

## Data Availability

Input data required to reproduce the analyses are stored in the 1_Data/ directory.

If any datasets are subject to usage restrictions, they are either:

 - Provided in processed form sufficient to reproduce results, or

 - Available from the original data providers upon reasonable request.

## Archived Version

The exact version of the code used for submission has been archived as:

Version: v1.0.0
DOI: (to be added after Zenodo archiving)

## Contact

Marius Munschek