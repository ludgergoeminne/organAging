### Here, we set the paths to all directories ###

# Specify your working directory
wd <- "../../"
# Specify the path where the scripts are saved
scripts.dir <- paste0(wd, "scripts/")

# Specify the path where the UK Biobank phenotypic data is located (downloaded from UK Biobank)
input.UKBB.dir <- paste0(wd, "data/input_UKBB/")
# Specify the path where the Olink data is located (downloaded from UK Biobank)
input.olink.dir <- paste0(wd, "/data/input_Olink/")
# Specify the path where the coding files are located (downloaded from UK Biobank)
# For example, download "coding143.tsv" from https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=143
input.coding.dir <- paste0(wd, "data/input_coding/")
# Specify the path where the GTEx data is located (downloaded from https://gtexportal.org/home/downloads/adult-gtex/overview)
input.GTEx.dir <- paste0(wd, "data/input_GTEx/")
# Specify the path where the Uniprot data is located (downloaded from https://www.uniprot.org/)
input.uniprot.dir <- paste0(wd, "data/input_Uniprot/")
# Specify the path where the COVID data (Filbin et al., (2021)) is located (downloaded from https://doi.org/10.17632/nf853r8xsj.2)
input.COVID.dir <- paste0(wd, "data/input_Filbin/")
# Specify the path where the Parkinson/MSA data (Dammmer et al., (2022)) is located (downloaded from https://www.synapse.org/#!Synapse:syn30549757/files/)
input.Parkinson.dir <- paste0(wd, "data/input_Dammmer/")
# Specify the path where the MESA data (Bild et al., (2002)) is located (downloaded from https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001416.v3.p1)
input.MESA.dir <- "data/input_Bild/"
# Specify the path where the sarcoidosis data (Damsky et al., (2022)) is located (downloaded from https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE169148)
input.GSE169148.dir <- paste0(wd, "data/input_Damsky/")

# Specify the path where you will save the rds files
rds.dir <- paste0(wd, "data/rds/")
# Specify the path where you will save the plots
plot.dir <- paste0(wd, "plots/")
# Specify the path where you will save the tables
table.dir <- paste0(wd, "tables/")
# Specify the directory where other aging models (e.g. PhenoAge, LocoAge), and Oh et al. statistics are stored
external.models.dir <- paste0(wd, "data/external_models/")

# Specify the name of your UK Biobank .tab file
ukbb.tab.file <- "ukb678668.tab"
# Specify the name of your Olink input file
olink.txt.file <- "olink_data.txt"
# Specify the name of the original Uniprot file
uniprot.file <- "uniprotkb_AND_model_organism_9606_AND_r_2024_02_27.tsv"
# Specify the name of the file where the phenotypic ages are stored
phenoage.file <- "ukb_phenoage.csv"
# Specify the name of the file where the locomotor ages are stored
locoage.file <- "ukb_locoage.csv"
# Specify the name of the original file that has the hazard ratio comparisons
# Downloaded from https://www.nature.com/articles/s41586-023-06802-1#MOESM5
ST8_Oh_et_al.file <- "Oh_et_al_ST8.xlsx"
# Specify the name of the original file that has the mortality summary statistics for other aging models
mortality.sumstats.file <- "other_papers_mortality_sumstats.csv"

# Specify where the GTEx_4x_FC_genes.json and GTEx_4x_FC_genes_feature_reduced.json object will be saved
dir.4x.FC.genes <- paste0(wd, "data/output_Python/")

# Specify the path where the input for training the full models (based on the data from the first visit, i.e. "instance 0") will be saved
dir.Python.input <- paste0(wd, "data/input_Python/instance_0/")
# Specify the path where the 1st-generation models (based on chronological age) will be saved
dir.gen1.models <- paste0(wd, "data/output_Python/instance_0/chronological_models/")
# Specify the path where the 2nd-generation models (based on mortality) will be saved
dir.gen2.models <- paste0(wd, "data/output_Python/instance_0/mortality_based_models/")
# Specify the path where the 2nd-generation AFT models will be saved
dir.AFT.models <- paste0(wd, "data/output_Python/instance_0/AFT_models/")
# Specify the path where the log files for the models will be saved
dir.log.files <- paste0(wd, "data/output_Python/instance_0/log_files/")

# Specify the path where the input for training the feature-reduced models ("longitudinal", as they will be used to assess longitudinality) will be saved
dir.Python.input.longitudinal <- paste0(wd, "data/input_Python/feature_reduced/")
# Specify the path where the 1st-generation feature-reduced models (based on chronological age) will be saved
dir.gen1.models.longitudinal <- paste0(wd, "data/output_Python/feature_reduced/chronological_models/")
# Specify the path where the 2nd-generation feature-reduced models (based on mortality) will be saved
dir.gen2.models.longitudinal <- paste0(wd, "data/output_Python/feature_reduced/mortality_based_models/")
# Specify the path where the log files for the feature-reduced models will be saved
dir.log.files.longitudinal <- paste0(wd, "data/output_Python/feature_reduced/log_files/")

# Specify the path where the input for training the models will be saved
dir.Python.input.CSF <- paste0(wd, "data/input_Python/CSF_first_visit/")
# Specify the path where the 1st-generation CSF models (based on chronological age) will be saved
dir.gen1.models.CSF <- paste0(wd, "data/output_Python/CSF_chronological_models/")
# Specify the path where the leave-one-out imputed data will be saved
dir.loo.imputed.data.CSF <- paste0(wd, "data/output_Python/CSF_LOO_imputed_data/")

# Here, we load the necessary R packages ###
library(tidyr)
library(caret)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(MatrixGenerics)
library(jsonlite)
library(kMeansEqual) # devtools::install_github("ludgergoeminne/kMeansEqual")
library(multcomp)
library(gratia)
library(ggpubr)
library(openxlsx)
library(geograbi)
library(reshape2)
library(broom)
library(ggbeeswarm)
library(GEOquery)
library(survival)
library(ggrepel)
library(nonnest2)
library(nonnestcox) # https://github.com/thomashielscher/nonnestcox
library(viridis)
library(fgsea)
library(annotate)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(data.table)
library(enrichR)
library(hypeR)
library(ggdendroplot)

#Source this helperscript to improve the svg files
source(paste0(scripts.dir, "Preprocess_SVG_file.R"))
