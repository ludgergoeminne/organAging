
###########################################################################################################################################
### This script helps with the installation of some of the packages that are needed if you would like to dive deep into this repository ###
###########################################################################################################################################

# Please do let us know if anything is missing in this script
# Credits to Seth Paulson for creating this script

# 1. Install CRAN packages
cran_packages <- c(
  "tidyr", "caret", "msigdbr", "dplyr", "ggplot2", "jsonlite", 
  "multcomp", "gratia", "openxlsx", "reshape2", "broom", "ggbeeswarm", 
  "survival", "ggrepel", "nonnest2", "data.table", "enrichR", "ggpubr"
)
missing_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(missing_cran)) install.packages(missing_cran)

# 2. Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
bioc_packages <- c("biomaRt", "MatrixGenerics", "GEOquery", "fgsea", "annotate", 
                   "org.Mm.eg.db", "org.Hs.eg.db", "hypeR")
BiocManager::install(bioc_packages, update = FALSE)

# 3. Install GitHub packages
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if(!("kMeansEqual" %in% installed.packages())){devtools::install_github("ludgergoeminne/kMeansEqual")}
if(!("nonnestcox" %in% installed.packages())){devtools::install_github("thomashielscher/nonnestcox")}
if(!("geograbi" %in% installed.packages())){devtools::install_github("yousefi138/geograbi")}
if(!("ggdendroplot" %in% installed.packages())){devtools::install_github("NicolasH2/ggdendroplot")}

