
#################################################################################################
### This script provides an example of how you can apply our aging models to your own dataset ###
#################################################################################################

#################
### 1. Set-up ###
#################

# Specify your working directory
wd <- ".../.../"

### 1.1. Object "organ.proteins" ###

# "organ.proteins" contains the names of all organ-specific aging models
# Download "organ.proteins" from https://github.com/ludgergoeminne/organAging/blob/main/data/rds/organ_proteins.rds:
organ.proteins <- readRDS(file = paste0(wd, "/data/rds/organ_proteins.rds"))
# The corresponding Python file is "GTEx_4x_FC_genes.json"
# This can be downloaded here: https://github.com/ludgergoeminne/organAging/blob/main/data/output_Python/GTEx_4x_FC_genes.json

### 1.2. Object "coefficients" ###

# "coefficients" contains the model coefficients
# Download "coefficients" from https://github.com/ludgergoeminne/organAging/blob/main/data/rds/coefficients.rds
coefficients <- readRDS(file = paste0(wd, "data/rds/coefficients.rds"))
# Alternatively, this object can also be reconstructed based on the csv files.
# The 1st-generation models csv files can be downloaded here: https://github.com/ludgergoeminne/organAging/tree/main/data/output_Python/instance_0/chronological_models
# The mortality-based models csv files can be downloaded here: https://github.com/ludgergoeminne/organAging/tree/main/data/output_Python//instance_0/mortality_based_models

# Alternative way to construct the object "coefficients":

### These are the full models: ###

# Specify the path where the 1st-generation models (based on chronological age) will be saved
dir.gen1.models <- paste0(wd, "data/output_Python/instance_0/chronological_models/")
# Specify the path where the 2nd-generation models (based on mortality) will be saved
dir.gen2.models <- paste0(wd, "data/output_Python/instance_0/mortality_based_models/")

### These are the feature-reduced models: ###

# # Specify the path where the 1st-generation feature-reduced models (based on chronological age) will be saved
# dir.gen1.models.longitudinal <- paste0(wd, "data/output_Python/longitudinal/chronological_models/")
# # Specify the path where the 2nd-generation feature-reduced models (based on mortality) will be saved
# dir.gen2.models.longitudinal <- paste0(wd, "data/output_Python/longitudinal/mortality_based_models/")

coefficients <- vector(mode = "list", length = 2)
names(coefficients) <- c("gen1", "gen2")
coefficients <- lapply(coefficients, function(x){
  x <- vector(mode = "list", length = length(organ.proteins))
  names(x) <- names(organ.proteins)
  return(x)
})

### We use the first fold of the full models ###
fold <- 1

for(g in 1:2){
  for(k in 1:length(organ.proteins)){
    if(!(names(organ.proteins)[k] %in% c("Bladder"))){
      if(g == 1){
        coefficients[["gen1"]][[k]] <- read.csv(paste0(dir.gen1.models, names(organ.proteins)[k], "_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
        coefficients[["gen1"]][[k]] <- unlist(coefficients[["gen1"]][[k]][fold, , drop = FALSE])
      } else if(g == 2){
        coefficients[["gen2"]][[k]] <- read.csv(paste0(dir.gen2.models, names(organ.proteins)[k], "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
        coefficients[["gen2"]][[k]] <- unlist(coefficients[["gen2"]][[k]][fold, , drop = FALSE])
      }
    }
  }
}

### 1.3. Standard deviations needed to rescale the protein NPX values ###

# "sds" contains the standard deviations for the full models needed to rescale the data.
# Download "sds" from: https://github.com/ludgergoeminne/organAging/blob/main/data/rds/standard_deviations.rds
sds <- readRDS(file = paste0(wd, "data/rds/standard_deviations.rds"))

# # Alternatively, the "sds" object can be downloaded in Excel format as Table S3 in our publication here: https://ars.els-cdn.com/content/image/1-s2.0-S1550413124004017-mmc4.xlsx
# sds <- openxlsx::read.xlsx("https://ars.els-cdn.com/content/image/1-s2.0-S1550413124004017-mmc4.xlsx", sheet = "Table S3")
# sds <- openxlsx::read.xlsx(paste0(wd, "UKBB/tables/Supplementary_Table_3.xlsx"), sheet = "Table S3")
# # Only keep the full models
# sds <- sds[sds$Type == "Full models", ]
# sds$Type <- NULL
# sds <- unlist(sds)
# sds <- sds[!is.na(sds)]

# "sds.longitudinal" contains the standard deviations for the feature-reduced models needed to rescale the data.
# Download "sds.longitudinal" from: https://github.com/ludgergoeminne/organAging/blob/main/data/rds/standard_deviations_longitudinal.rds
# sds.longitudinal <- readRDS(file = paste0(wd, "data/rds/standard_deviations_longitudinal.rds"))

# # Alternatively, the "sds" object can be downloaded in Excel format as Table S3 in our publication here: https://ars.els-cdn.com/content/image/1-s2.0-S1550413124004017-mmc4.xlsx
# sds.longitudinal <- openxlsx::read.xlsx("https://ars.els-cdn.com/content/image/1-s2.0-S1550413124004017-mmc4.xlsx", sheet = "Table S3")
# # sds.longitudinal <- openxlsx::read.xlsx(paste0(wd, "tables/Supplementary_Table_3.xlsx"), sheet = "Table S3")
# # Only keep the feature-reduced models
# sds.longitudinal <- sds.longitudinal[sds.longitudinal$Type == "Feature-reduced models", ]
# sds.longitudinal$Type <- NULL
# sds.longitudinal <- unlist(sds.longitudinal)
# sds.longitudinal <- sds.longitudinal[!is.na(sds.longitudinal)]

### 1.4. Values needed to convert the output from the mortality-based models to years ###

# From https://github.com/ludgergoeminne/organAging/, to convert rel. log(mortality hazards) to ages in years
avg.rel.log.mort.hazard <- -4.801912
intercept <- -9.94613787413831
slope <- 0.0897860500778604

### 1.5. Import your data ###

# Your data should contain the protein names as columns
df.full <- openxlsx::read.xlsx(paste0(wd, "path_to_your_data/yourdata.xlsx"))
df.full[1:10, 1:10] # Inspect your data

# Specify the name of the column that contains the chronological ages in your dataset
# This is only needed to calculate the age deviations; indicated as "residual" below
age.col <- "age"

# Check protein column names, in case there are issues
colnames(df.full)[grepl(".", colnames(df.full), fixed = TRUE)]

# Fix possible issues
# colnames(df.full)[colnames(df.full) == "NT.proBNP"] <- "NTproBNP"
# colnames(df.full)[colnames(df.full) == "HLA.DRA"] <- "HLA-DRA"
# colnames(df.full)[colnames(df.full) == "HLA.E"] <- "HLA-E"
# colnames(df.full)[colnames(df.full) == "ERVV.1"] <- "ERVV-1"
# colnames(df.full)[colnames(df.full) == "HLA.A"] <- "HLA-A"

# Check if you have missing proteins
organ.proteins$Conventional[!(organ.proteins$Conventional %in% colnames(df.full))]

# Check if you have some proteins that we don't have
colnames(df.full)[!(colnames(df.full) %in% organ.proteins$Conventional)]

# Only keep the proteins on which the models were trained, drop possible non-Olink columns
df.olink <- df.full[, colnames(df.full) %in% names(sds)]

###############################################################################################
### 2. Optional, but worth checking for possible better performance: standardize Olink data ###
###############################################################################################

# If you don't have Olink data, standardization is highly recommended, as it will almost always improve performance.

tmp.sds <- sds[names(sds) %in% colnames(df.olink)]
tmp.sds <- tmp.sds[colnames(df.olink)]
all(names(tmp.sds) == colnames(df.olink))
# TRUE

df.olink <- apply(df.olink, 2, function(x){
  x <- x/sd(x, na.rm = TRUE)
})

df.olink <- sweep(df.olink, 2, tmp.sds, FUN="*")

#####################################
### 3. Calculate the aging models ###
#####################################

df.full[, colnames(df.olink)] <- df.olink

predicted.ages <- vector(mode = "list", length = 2)
names(predicted.ages) <- c("gen1", "gen2")
predicted.ages <- lapply(predicted.ages, function(x){
  x <- vector(mode = "list", length = 2)
  names(x) <- c("predicted", "residual")
  x <- lapply(x, function(z){
    z <- vector(mode = "list", length = length(organ.proteins))
    names(z) <- names(organ.proteins)
    return(z)
  })
  return(x)
})

### 3.1. 1st generation aging models ("gen1") ###

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder"))){
    coefs <- coefficients[["gen1"]][[k]][names(coefficients[["gen1"]][[k]])[names(coefficients[["gen1"]][[k]]) %in% colnames(df.full)]]
    # Calculate the biological ages ("predicted")
    predicted.ages[["gen1"]][["predicted"]][[k]] <- coefficients[["gen1"]][[k]]["Intercept"]+rowSums(sweep(as.matrix(df.full[, names(coefs), drop = FALSE]), MARGIN=2, coefs, `*`))
    
    # Calculate the age deviations ("residual"); remove this line if you do not have chronological ages available; then it is not possible to calculate age deviations.
    predicted.ages[["gen1"]][["residual"]][[k]] <- c(resid(lm(predicted.ages[["gen1"]][["predicted"]][[k]]~df.full[, age.col])))
  }
}

### 3.2. Mortality-based (2nd generation) aging models ("gen2") ###

for(k in 1:length(organ.proteins)){
  if(names(organ.proteins)[k] != "Bladder"){
    coefs <- coefficients[["gen2"]][[k]][names(coefficients[["gen2"]][[k]])[names(coefficients[["gen2"]][[k]]) %in% colnames(df.full)]]
    # Calculate the biological ages ("predicted")
    predicted.ages[["gen2"]][["predicted"]][[k]] <- rowSums(sweep(as.matrix(df.full[, names(coefs), drop = FALSE]), MARGIN=2, coefs, `*`))
    
    # Calculate the age deviations ("residual"); remove this line if you do not have chronological ages available; then it is not possible to calculate age deviations.
    predicted.ages[["gen2"]][["residual"]][[k]] <- c(resid(lm(predicted.ages[["gen2"]][["predicted"]][[k]]~df.full[, age.col])))
  }
}
