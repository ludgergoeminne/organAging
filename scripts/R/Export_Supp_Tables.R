
### "Male", "Female", "Organismal", "Multi-organ" are not being discussed in the manuscript
organ.proteins.1B <- organ.proteins[!(names(organ.proteins) %in% c("Male", "Female", "Organismal", "Multi-organ"))]

### This comes from 1_Prepare_Data.R ###
sds <- readRDS(file = paste0(rds.dir, "standard_deviations.rds"))

### This comes 2_Prepare_Longitudinal_Data.R ###
sds.longitudinal <- readRDS(file = paste0(rds.dir, "standard_deviations_longitudinal.rds"))

### This comes from 3_Processe_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))
organ.proteins.longitudinal <- readRDS(paste0(rds.dir, "organ_proteins_longitudinal.rds"))

### This comes from Plot_Volcanos_GSEA. R ###
res.age <- readRDS(file = paste0(rds.dir, "res_age.rds"))
res.mortality <- readRDS(file = paste0(rds.dir, "res_mortality.rds"))
fgsea.age <- readRDS(file = paste0(rds.dir, "fgsea_age.rds"))
fgsea.mortality <- readRDS(file = paste0(rds.dir, "fgsea_mortality.rds"))

### 1. Table S1: 1st-generation models and mortality-based models ###

### A. Coefficients of the 1st-generation models ###
tmp.list <- vector(mode = "list", length = length(organ.proteins.1B))
names(tmp.list) <- names(organ.proteins.1B)

tmp.list <- lapply(tmp.list, function(x){
  x <- data.frame(matrix(ncol = (length(organ.proteins$Conventional)+1), nrow = 5))
  colnames(x) <- c("Intercept", organ.proteins$Conventional)
  return(x)
})

for(k in 1:length(organ.proteins.1B)){
  if(!(names(organ.proteins.1B)[k] %in% c("Bladder"))){
    tmp.coefs <- read.csv(paste0(dir.gen1.models, names(organ.proteins.1B)[k], "_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    if(!all(colnames(tmp.coefs) %in% colnames(tmp.list[[k]]))){stop("There are some missing proteins!!!")}
    tmp.list[[k]][, colnames(tmp.coefs)] <- tmp.coefs
    tmp.list[[k]] <- cbind(data.frame(
      Model = names(organ.proteins.1B)[k],
      Fold = 1:5
    ), tmp.list[[k]])
  }
}
tmp.list[["Bladder"]] <- NULL

tmp.A <- do.call("rbind", tmp.list)

# Get range of non-zero coefficients
rowSums(tmp.A[, -c(1:3)] != 0, na.rm = TRUE)
range(c(2315, 2288, 2306, 2310, 2308))
# 2288 2315

### B. Model parameters of the 1st-generation models ###

tmp.list <- vector(mode = "list", length = 5)
for(i in 1:length(tmp.list)){
  temp <- readLines(paste0(dir.log.files, "chronological_kfold_", i, ".out"))
  tmp.list[[i]] <- data.frame(Model = temp[seq(1, length(temp), by = 2)],
                              Fold = i,
                              `Number of proteins` = NA,
                              `Best L1 ratio` = as.numeric(gsub("Best l1_ratio ", "", temp[seq(2, length(temp), by = 2)])),
                              check.names = FALSE)
}

tmp.B <- do.call("rbind", tmp.list)
tmp.B <- tmp.B[tmp.B$Model %in% names(organ.proteins.1B),]
tmp.B$tmp.col <- paste0(tmp.B$Model, ".", tmp.B$Fold)

tmp.df <- data.frame(
  tmp.col = names(rowSums(tmp.A[, -c(1:3)] != 0, na.rm = TRUE)),
  `Number of proteins` = rowSums(tmp.A[, -c(1:3)] != 0, na.rm = TRUE),
  check.names = FALSE
)
tmp.B <- tmp.B[match(tmp.df$tmp.col, tmp.B$tmp.col),]

tmp.B$`Number of proteins` <- tmp.df$`Number of proteins`
tmp.B$tmp.col <- NULL

### C. Coefficients of the mortality-based models ###

tmp.list <- vector(mode = "list", length = length(organ.proteins.1B))
names(tmp.list) <- names(organ.proteins.1B)

tmp.list <- lapply(tmp.list, function(x){
  x <- data.frame(matrix(ncol = length(organ.proteins$Conventional), nrow = 5))
  colnames(x) <- organ.proteins$Conventional
  return(x)
})

for(k in 1:length(organ.proteins.1B)){
  if(!(names(organ.proteins.1B)[k] %in% c("Bladder"))){
    tmp.coefs <- read.csv(paste0(dir.gen2.models, names(organ.proteins.1B)[k], "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    if(!all(colnames(tmp.coefs) %in% colnames(tmp.list[[k]]))){stop("There are some missing proteins!!!")}
    tmp.list[[k]][, colnames(tmp.coefs)] <- tmp.coefs
    tmp.list[[k]] <- cbind(data.frame(
      Model = names(organ.proteins.1B)[k],
      Fold = 1:5
    ), tmp.list[[k]])
  }
}
tmp.list[["Bladder"]] <- NULL

tmp.C <- do.call("rbind", tmp.list)

### D. Model parameters of the mortality-based models ###

tmp.list <- vector(mode = "list", length = 5)
for(i in 1:length(tmp.list)){
  temp <- readLines(paste0(dir.log.files, "mortality_kfold_", i, ".out"))
  tmp.list[[i]] <- data.frame(Model = temp[seq(1, length(temp), by = 3)],
                              Fold = i,
                              `Number of proteins` = NA,                              
                              `Best L1 ratio` = as.numeric(gsub("Best l1_ratio ", "", temp[seq(3, length(temp), by = 3)])),
                              `Best alpha` = as.numeric(gsub("\\]", "", gsub("Best alpha \\[", "", temp[seq(2, length(temp), by = 3)]))), 
                              check.names = FALSE)
}

tmp.D <- do.call("rbind", tmp.list)
tmp.D <- tmp.D[tmp.D$Model %in% names(organ.proteins.1B),]
tmp.D$tmp.col <- paste0(tmp.D$Model, ".", tmp.D$Fold)

tmp.df <- data.frame(
  tmp.col = names(rowSums(tmp.C[, -c(1:2)] != 0, na.rm = TRUE)),
  `Number of proteins` = rowSums(tmp.C[, -c(1:2)] != 0, na.rm = TRUE),
  check.names = FALSE
)
tmp.D <- tmp.D[match(tmp.df$tmp.col, tmp.D$tmp.col),]

tmp.D$`Number of proteins` <- tmp.df$`Number of proteins`
tmp.D$tmp.col <- NULL

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "Table S1. Model coefficients and hyperparameters for the full conventional and organ-specific aging models, related to STAR Methods. 
A.) Model coefficients for the 5 folds of the 1st-generation conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. 
B.) Model hyperparameters for the 5 folds of the 1st-generation conventional and organ-specific models. 
C.) Model coefficients for the 5 folds of the mortality-based conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. 
D.) Model hyperparameters for the 5 folds of the mortality-based conventional and organ-specific models.
") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Table S1A")
openxlsx::writeData(wb, "Table S1A", tmp.A, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Table S1B")
setColWidths(wb, sheet ="Table S1B", cols = 1:ncol(tmp.B), widths = c(10.78, 7.78, 16.89, 10.78))
openxlsx::writeData(wb, "Table S1B", tmp.B, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Table S1C")
openxlsx::writeData(wb, "Table S1C", tmp.C, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Table S1D")
setColWidths(wb, sheet ="Table S1D", cols = 1:ncol(tmp.D), widths = c(10.78, 7.78, 16.89, 10.78, 10.78))
openxlsx::writeData(wb, "Table S1D", tmp.D, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Table_S1.xlsx"))

### 2. Table S2: model metrics ###

### This comes from Plot_Heatmaps_r_MAE.R ###
summary.df.list <- readRDS(file = paste0(rds.dir, "summary_df_list.rds"))

### This comes from Plot_Heatmaps_r_MAE_longitudinal.R ###
summary.df.list.longitudinal <- readRDS(file = paste0(rds.dir, "summary_df_list_longitudinal.rds"))

for(i in 1:length(summary.df.list)){
  summary.df.list[[i]] <- cbind(data.frame(model = rep(names(summary.df.list)[i], nrow(summary.df.list[[i]]))), summary.df.list[[i]])
}

tmp.A <- do.call("rbind", summary.df.list)

for(i in 1:length(summary.df.list)){
  summary.df.list.longitudinal[[i]] <- cbind(data.frame(model = rep(names(summary.df.list.longitudinal)[i], nrow(summary.df.list.longitudinal[[i]]))), summary.df.list.longitudinal[[i]])
}

tmp.B <- do.call("rbind", summary.df.list.longitudinal)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "Table S2. Correlation and accuracy measures of biological ages predicted by conventional and organ-specific models, related to Figures 1, 3, and 6.
A.) Correlation and accuracy measures of biological ages predicted by the full conventional and organ-specific 1st-generation and mortality-based models. 
B.) Correlation and accuracy measures of biological ages predicted by the feature-reduced conventional and organ-specific 1st-generation and mortality-based models.
r: correlation coefficient, r²: coefficient of determination, MAE: mean absolute error of the residuals, MSE: mean squared error, pval: p-value.
") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Table S2A")
openxlsx::writeData(wb, "Table S2A", tmp.A, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Table S2B")
openxlsx::writeData(wb, "Table S2B", tmp.B, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Table_S2.xlsx"))

### 3. Table S3: standard deviations ###

tmp <- data.frame(matrix(ncol = length(unique(c(names(sds), names(sds.longitudinal)))), nrow = 2))
colnames(tmp) <- sort(unique(c(names(sds), names(sds.longitudinal))))

tmp[1, names(sds)] <- sds
tmp[2, names(sds.longitudinal)] <- sds.longitudinal

tmp <- cbind(data.frame(`Type` = c("Full models", "Feature-reduced models")), tmp)

wb <- createWorkbook()
# Add caption worksheet
openxlsx::addWorksheet(wb, "Caption")
# Add column widths to Spreadsheet
# setColWidths(wb, sheet = "caption", cols = 1:2, widths = c(34.56, 82.44))
# create style, in this case bold header
# header_st <- createStyle(textDecoration = "Bold")
# Write descrition data frame to new worksheet
openxlsx::writeData(wb, "Caption", "Table S3. Standard deviations of the Olink NPX values in the datasets used to train the full models and the feature-reduced models, related to STAR Methods. 
Empty values indicate that the protein was not considered for model building.") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")
openxlsx::addWorksheet(wb, "Table S3")
openxlsx::writeData(wb, "Table S3", tmp, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Table_S3.xlsx"))

### 4. Table S4: Associations with chronological age and associations with mortality ###

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "Table S4. Effects on chronological age and mortality for 2,923 proteins and their associated pathways in 53,016 UK Biobank participants, related to Figure 2. 
A.) Effects on chronological age after correction for sex.
B.) Gene set enrichment analysis using MSigDB pathways on the Z-values in Table S4A.
C.) Effects on mortality after correction for sex, chronological age, and the interaction between sex and chronological age.
D.) Gene set enrichment analysis using MSigDB pathways on the Z-values in Table S4C.
ln: natural logarithm, pval: enrichment p-value, padj: Benjamini-Hochberg-adjusted p-value, log2err: the expected error for the standard deviation of the log2(p-value), ES: enrichment score, NES: normalized enrichment score, size: size of the pathway after removing genes not present in the input, leadingEdge: vector with indexes of leading edge genes that drive the enrichment, zval: Z-value.
") # , headerStyle = header_st

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
header_st <- createStyle(textDecoration = "Bold")

### A. Associations with chronological age ###

openxlsx::addWorksheet(wb, "Table S4A")

tmp.A <- res.age
tmp.A$significant <- tmp.A$qval < 0.05
tmp.A$zval <- limma::zscoreT(tmp.A$`t value`, df = tmp.A$df, approx=FALSE)
tmp.A <- tmp.A %>% arrange(-abs(zval))
tmp.A$zval <- NULL
colnames(tmp.A)[3] <- "protein name"
colnames(tmp.A)[4] <- "estimate"
colnames(tmp.A)[5] <- "standard error"
colnames(tmp.A)[6] <- "degrees of freedom"
colnames(tmp.A)[8] <- "p-value"
colnames(tmp.A)[9] <- "lower bound 95% CI"
colnames(tmp.A)[10] <- "upper bound 95% CI"
colnames(tmp.A)[11] <- "correlation coefficient (r)"
colnames(tmp.A)[12] <- "q-value"
colnames(tmp.A)

# Add numeric formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "NUMBER"), rows = 2:(nrow(tmp.A)+1), cols = c(4:5, 7, 9:11), sheet = "Table S4A", gridExpand = TRUE)

# Add scientific formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "SCIENTIFIC"), rows = 2:(nrow(tmp.A)+1), cols = c(8, 12), sheet = "Table S4A", gridExpand = TRUE)

# Add column widths to Spreadsheet
setColWidths(wb, sheet ="Table S4A", cols = 1:ncol(tmp.A), widths = c(rep(10.78, 2), 40.11, rep(10.78, 10)))

# Add conditional formatting
conditionalFormatting(wb, sheet = "Table S4A", cols = 4, rows = 2:(nrow(tmp.A)+1), rule = "<0", style = negStyle)
conditionalFormatting(wb, sheet = "Table S4A", cols = 4, rows = 2:(nrow(tmp.A)+1), rule = ">0", style = posStyle)

openxlsx::writeData(wb, "Table S4A", tmp.A, headerStyle = header_st)

### B. Gene set enrichment analysis for the associations with chronological age ###

openxlsx::addWorksheet(wb, "Table S4B")

tmp.B <- fgsea.age
tmp.B$significant <- tmp.B$padj < 0.05

# Add numeric formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "NUMBER"), rows = 2:(nrow(tmp.B)+1), cols = c(4:7, 9), sheet = "Table S4B", gridExpand = TRUE)

# Add scientific formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "SCIENTIFIC"), rows = 2:(nrow(tmp.B)+1), cols = c(2, 3), sheet = "Table S4B", gridExpand = TRUE)

# Add column widths to Spreadsheet
setColWidths(wb, sheet ="Table S4B", cols = 1:ncol(tmp.B), widths = c(54.11, rep(7.78, 6), 30.89, 7.78, 10.78))

# Add conditional formatting
conditionalFormatting(wb, sheet = "Table S4B", cols = 6, rows = 2:(nrow(tmp.B)+1), rule = "<1", style = negStyle)
conditionalFormatting(wb, sheet = "Table S4B", cols = 6, rows = 2:(nrow(tmp.B)+1), rule = ">1", style = posStyle)

openxlsx::writeData(wb, "Table S4B", tmp.B, headerStyle = header_st)

### C. Associations with with mortality ###

openxlsx::addWorksheet(wb, "Table S4C")

tmp.C <- res.mortality
tmp.C <- tmp.C %>% arrange(-abs(z))
tmp.C$significant <- tmp.C$qval < 0.05
colnames(tmp.C)[3] <- "protein name"
colnames(tmp.C)[4] <- "hazard ratio"
colnames(tmp.C)[5] <- "ln(hazard ratio)"
colnames(tmp.C)[6] <- "z-value"
colnames(tmp.C)[7] <- "p-value"
colnames(tmp.C)[8] <- "lower bound 95% CI"
colnames(tmp.C)[9] <- "upper bound 95% CI"
colnames(tmp.C)[10] <- "q-value"
colnames(tmp.C)

# Add numeric formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "NUMBER"), rows = 2:(nrow(tmp.C)+1), cols = c(4:6, 8:9), sheet = "Table S4C", gridExpand = TRUE)

# Add scientific formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "SCIENTIFIC"), rows = 2:(nrow(tmp.C)+1), cols = c(7, 10), sheet = "Table S4C", gridExpand = TRUE)

# Add column widths to Spreadsheet
setColWidths(wb, sheet ="Table S4C", cols = 1:ncol(tmp.C), widths = c(rep(10.78, 2), 40.11, 10.78, 13.33, rep(10.78, 6)))

# Add conditional formatting
conditionalFormatting(wb, sheet = "Table S4C", cols = 4, rows = 2:(nrow(tmp.C)+1), rule = "<1", style = negStyle)
conditionalFormatting(wb, sheet = "Table S4C", cols = 4, rows = 2:(nrow(tmp.C)+1), rule = ">1", style = posStyle)

openxlsx::writeData(wb, "Table S4C", tmp.C, headerStyle = header_st)

### D. Gene set enrichment analysis for the associations with chronological age ###

openxlsx::addWorksheet(wb, "Table S4D")

tmp.D <- fgsea.mortality
tmp.D$significant <- tmp.D$padj < 0.05

# Add numeric formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "NUMBER"), rows = 2:(nrow(tmp.D)+1), cols = c(4:7, 9), sheet = "Table S4D", gridExpand = TRUE)

# Add scientific formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "SCIENTIFIC"), rows = 2:(nrow(tmp.D)+1), cols = c(2, 3), sheet = "Table S4D", gridExpand = TRUE)

# Add column widths to Spreadsheet
setColWidths(wb, sheet ="Table S4D", cols = 1:ncol(tmp.D), widths = c(54.11, rep(7.78, 6), 30.89, 7.78, 10.78))

# Add conditional formatting
conditionalFormatting(wb, sheet = "Table S4D", cols = 6, rows = 2:(nrow(tmp.D)+1), rule = "<1", style = negStyle)
conditionalFormatting(wb, sheet = "Table S4D", cols = 6, rows = 2:(nrow(tmp.D)+1), rule = ">1", style = posStyle)

openxlsx::writeData(wb, "Table S4D", tmp.D, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Table_S4.xlsx"))

### 5. Table S5: Feature-reduced 1st-generation and mortality-based models ###

### A. Coefficients of the feature-reduced 1st-generation models ###

tmp.list <- vector(mode = "list", length = length(organ.proteins.1B))
names(tmp.list) <- names(organ.proteins.1B)

tmp.list <- lapply(tmp.list, function(x){
  x <- data.frame(matrix(ncol = (length(organ.proteins.longitudinal$Conventional)+1), nrow = 5))
  colnames(x) <- c("Intercept", organ.proteins.longitudinal$Conventional)
  return(x)
})

for(k in 1:length(organ.proteins.1B)){
  if(!(names(organ.proteins.1B)[k] %in% c("Bladder", "Thyroid"))){
    tmp.coefs <- read.csv(paste0(dir.gen1.models.longitudinal, names(organ.proteins.1B)[k], "_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
    if(!all(colnames(tmp.coefs) %in% colnames(tmp.list[[k]]))){stop("There are some missing proteins!!!")}
    tmp.list[[k]][, colnames(tmp.coefs)] <- tmp.coefs
    tmp.list[[k]] <- cbind(data.frame(
      Model = names(organ.proteins.1B)[k],
      Fold = 1:5
    ), tmp.list[[k]])
  }
}
tmp.list[["Bladder"]] <- NULL
tmp.list[["Thyroid"]] <- NULL

tmp.A <- do.call("rbind", tmp.list)

### B. Model hyperparameters of the feature-reduced 1st-generation models ###

tmp.list <- vector(mode = "list", length = 5)
for(i in 1:length(tmp.list)){
  temp <- readLines(paste0(dir.log.files.longitudinal, "reduced_chrono_kfold_", i, ".out"))
  tmp.list[[i]] <- data.frame(Model = temp[seq(1, length(temp), by = 2)],
                              Fold = i,
                              `Number of proteins` = NA,
                              `Best L1 ratio` = as.numeric(gsub("Best l1_ratio ", "", temp[seq(2, length(temp), by = 2)])),
                              check.names = FALSE)
}

tmp.B <- do.call("rbind", tmp.list)
tmp.B <- tmp.B[tmp.B$Model %in% names(organ.proteins.1B),]
tmp.B$tmp.col <- paste0(tmp.B$Model, ".", tmp.B$Fold)

tmp.df <- data.frame(
  tmp.col = names(rowSums(tmp.A[, -c(1:3)] != 0, na.rm = TRUE)),
  `Number of proteins` = rowSums(tmp.A[, -c(1:3)] != 0, na.rm = TRUE),
  check.names = FALSE
)
tmp.B <- tmp.B[match(tmp.df$tmp.col, tmp.B$tmp.col),]

tmp.B$`Number of proteins` <- tmp.df$`Number of proteins`
tmp.B$tmp.col <- NULL

### C. Coefficients of the feature-reduced mortality-based models ###

tmp.list <- vector(mode = "list", length = length(organ.proteins.1B))
names(tmp.list) <- names(organ.proteins.1B)

tmp.list <- lapply(tmp.list, function(x){
  x <- data.frame(matrix(ncol = length(organ.proteins.longitudinal$Conventional), nrow = 5))
  colnames(x) <- organ.proteins.longitudinal$Conventional
  return(x)
})

for(k in 1:length(organ.proteins.1B)){
  if(!(names(organ.proteins.1B)[k] %in% c("Bladder", "Thyroid"))){
    tmp.coefs <- read.csv(paste0(dir.gen2.models.longitudinal, names(organ.proteins.1B)[k], "_mortality_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
    if(!all(colnames(tmp.coefs) %in% colnames(tmp.list[[k]]))){stop("There are some missing proteins!!!")}
    tmp.list[[k]][, colnames(tmp.coefs)] <- tmp.coefs
    tmp.list[[k]] <- cbind(data.frame(
      Model = names(organ.proteins.1B)[k],
      Fold = 1:5
    ), tmp.list[[k]])
  }
}
tmp.list[["Bladder"]] <- NULL
tmp.list[["Thyroid"]] <- NULL

tmp.C <- do.call("rbind", tmp.list)

### D. Model hyperparameters of the feature-reduced mortality-based models ###

tmp.list <- vector(mode = "list", length = 5)
for(i in 1:length(tmp.list)){
  temp <- readLines(paste0(dir.log.files.longitudinal, "reduced_mortality_kfold_", i, ".out"))
  tmp.list[[i]] <- data.frame(Model = temp[seq(1, length(temp), by = 3)],
                              Fold = i,
                              `Number of proteins` = NA,                              
                              `Best L1 ratio` = as.numeric(gsub("Best l1_ratio ", "", temp[seq(3, length(temp), by = 3)])),
                              `Best alpha` = as.numeric(gsub("\\]", "", gsub("Best alpha \\[", "", temp[seq(2, length(temp), by = 3)]))), 
                              check.names = FALSE)
}

tmp.D <- do.call("rbind", tmp.list)
tmp.D <- tmp.D[tmp.D$Model %in% names(organ.proteins.1B),]
tmp.D$tmp.col <- paste0(tmp.D$Model, ".", tmp.D$Fold)

tmp.df <- data.frame(
  tmp.col = names(rowSums(tmp.C[, -c(1:2)] != 0, na.rm = TRUE)),
  `Number of proteins` = rowSums(tmp.C[, -c(1:2)] != 0, na.rm = TRUE),
  check.names = FALSE
)
tmp.D <- tmp.D[match(tmp.df$tmp.col, tmp.D$tmp.col),]

tmp.D$`Number of proteins` <- tmp.df$`Number of proteins`
tmp.D$tmp.col <- NULL

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "Table S5. Model coefficients and hyperparameters for the feature-reduced conventional and organ-specific aging models, related to STAR Methods. 
A.) Model coefficients for the 5 folds of the feature-reduced 1st-generation conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. 
B.) Model hyperparameters for the 5 folds of the feature-reduced 1st-generation conventional and organ-specific models. 
C.) Model coefficients for the 5 folds of the feature-reduced mortality-based conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. 
D.) Model hyperparameters for the 5 folds of the feature-reduced mortality-based conventional and organ-specific models.
") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Table S5A")
openxlsx::writeData(wb, "Table S5A", tmp.A, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Table S5B")
# Add column widths to Spreadsheet for S5B
setColWidths(wb, sheet ="Table S5B", cols = 1:ncol(tmp.B), widths = c(10.78, 7.78, 16.89, 10.78))
openxlsx::writeData(wb, "Table S5B", tmp.B, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Table S5C")
openxlsx::writeData(wb, "Table S5C", tmp.C, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Table S5D")
# Add column widths to Spreadsheet for S5D
setColWidths(wb, sheet ="Table S5D", cols = 1:ncol(tmp.D), widths = c(10.78, 7.78, 16.89, 10.78, 10.78))
openxlsx::writeData(wb, "Table S5D", tmp.D, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Table_S5.xlsx"))

### 6. Table S6: 1st-generation CSF model ###
tmp.list <- vector(mode = "list", length = length(organ.proteins.1B))
names(tmp.list) <- names(organ.proteins.1B)

tmp.list <- lapply(tmp.list, function(x){
  x <- data.frame(matrix(ncol = (length(organ.proteins.1B$Conventional)+1), nrow = 184))
  colnames(x) <- c("Intercept", organ.proteins.1B$Conventional)
  return(x)
})

for(k in 1:length(organ.proteins.1B)){
  if(!(names(organ.proteins.1B)[k] %in% c("Bladder", "Thyroid"))){
    tmp.list[[k]][, c("Intercept", organ.proteins.1B[[k]])] <- read.csv(paste0(dir.gen1.models.CSF, names(organ.proteins.1B)[k], "_coefs.csv"), header = TRUE, check.names = FALSE)
    tmp.list[[k]] <- cbind(data.frame(
      Model = names(organ.proteins.1B)[k],
      `Sample left out` = 1:184, check.names = FALSE
    ), tmp.list[[k]])
  }
}
tmp.list[["Bladder"]] <- NULL
tmp.list[["Thyroid"]] <- NULL

tmp <- do.call("rbind", tmp.list)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "Table S6. Model coefficients for the 1st-generation conventional and organ-specific models trained on 1,031 proteins measured in the CSF in the Dammer et al. (2022) dataset, related to STAR Methods.
Model coefficients for the 1st-generation conventional and organ-specific models trained on 1,031 proteins measured in the CSF in the Dammer et al. (2022) dataset. Proteins that were not reliably measured in the CSF or are not specific for the corresponding organ have empty values.
") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Table S6")
openxlsx::writeData(wb, "Table S6", tmp, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Table_S6.xlsx"))

