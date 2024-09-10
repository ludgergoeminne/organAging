
### "Male", "Female", "Organismal", "Multi-organ" are not being discussed in the manuscript
organ.proteins.1B <- organ.proteins[!(names(organ.proteins) %in% c("Male", "Female", "Organismal", "Multi-organ"))]

### This comes from 3_Processe_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))
organ.proteins.longitudinal <- readRDS(paste0(rds.dir, "organ_proteins_longitudinal.rds"))

### 1. Supp. Table 1: 1st-generation models and mortality-based models ###

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
openxlsx::writeData(wb, "Caption", "Supplementary Table 1. A. Model coefficients for the 5 folds of the 1st-generation conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. B. Model parameters for the 5 folds of the 1st-generation conventional and organ-specific models. C. Model coefficients for the 5 folds of the mortality-based conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. D. Model parameters for the 5 folds of the mortality-based conventional and organ-specific models.") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Supp. Table 1A")
openxlsx::writeData(wb, "Supp. Table 1A", tmp.A, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Supp. Table 1B")
openxlsx::writeData(wb, "Supp. Table 1B", tmp.B, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Supp. Table 1C")
openxlsx::writeData(wb, "Supp. Table 1C", tmp.C, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Supp. Table 1D")
openxlsx::writeData(wb, "Supp. Table 1D", tmp.D, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Supplementary_Table_1.xlsx"))

### 2. Supp. Table 2: model metrics ###

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
openxlsx::writeData(wb, "Caption", "Supplementary Table 2. A. Correlation and accuracy measures of biological ages and relative log(mortality hazards) predicted by conventional and organ-specific models. B. Correlation and accuracy measures of biological ages and relative log(mortality hazards) predicted by conventional and organ-specific feature-reduced models.") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Supp. Table 2A")
openxlsx::writeData(wb, "Supp. Table 2A", tmp.A, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Supp. Table 2B")
openxlsx::writeData(wb, "Supp. Table 2B", tmp.B, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Supplementary_Table_2.xlsx"))

### 3. Supp. Table 3: standard deviations ###

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
openxlsx::writeData(wb, "Caption", "Supplementary Table 3. Standard deviations of the Olink NPX values in the datasets used to train the full models and the feature-reduced models. Empty values indicate that the protein was not considered for model building.") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")
openxlsx::addWorksheet(wb, "Supp. Table 3")
openxlsx::writeData(wb, "Supp. Table 3", tmp, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Supplementary_Table_3.xlsx"))

### 4. Supp. Table 4: associations with chronological age ###

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "Supplementary Table 4. Effects on chronological age after correction for sex for 2,923 proteins in 53,016 UK Biobank participants.") # , headerStyle = header_st

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
header_st <- createStyle(textDecoration = "Bold")
openxlsx::addWorksheet(wb, "Supp. Table 4")

tmp <- res.age
tmp$significant <- tmp$qval < 0.05
tmp$zval <- limma::zscoreT(tmp$`t value`, df = tmp$df, approx=FALSE)
tmp <- tmp %>% arrange(-abs(zval))
tmp$zval <- NULL
colnames(tmp)[3] <- "protein name"
colnames(tmp)[5] <- "standard error"
colnames(tmp)[6] <- "degrees of freedom"
colnames(tmp)[8] <- "p-value"
colnames(tmp)[9] <- "lower bound 95% CI"
colnames(tmp)[10] <- "upper bound 95% CI"
colnames(tmp)[11] <- "correlation coefficient (r)"
colnames(tmp)[12] <- "q-value"
colnames(tmp)

# Add numeric formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "NUMBER"), rows = 2:(nrow(tmp)+1), cols = c(4:5, 7, 9:11), sheet = "Supp. Table 4", gridExpand = TRUE)

# Add scientific formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "SCIENTIFIC"), rows = 2:(nrow(tmp)+1), cols = c(8, 12), sheet = "Supp. Table 4", gridExpand = TRUE)

# Add column widths to Spreadsheet
setColWidths(wb, sheet ="Supp. Table 4", cols = 1:ncol(tmp), widths = c(rep(10.78, 2), 40.11, rep(10.78, 10)))

# Add conditional formatting
conditionalFormatting(wb, sheet = "Supp. Table 4", cols = 4, rows = 2:(nrow(tmp)+1), rule = "<0", style = negStyle)
conditionalFormatting(wb, sheet = "Supp. Table 4", cols = 4, rows = 2:(nrow(tmp)+1), rule = ">0", style = posStyle)

openxlsx::writeData(wb, "Supp. Table 4", tmp, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Supplementary_Table_4.xlsx"))

### 5. Supp. Table 5: associations with mortality ###

wb <- createWorkbook()
# Add caption worksheet
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "Supplementary Table 5. Effects on mortality after correction for sex, chronological age, and the interaction between sex and chronological age for 2,923 proteins in 53,016 UK Biobank participants. \"ln\" denotes the natural logarithm. ") # , headerStyle = header_st

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
header_st <- createStyle(textDecoration = "Bold")
openxlsx::addWorksheet(wb, "Supp. Table 5")

tmp <- res.mortality
tmp <- tmp %>% arrange(-abs(z))
tmp$significant <- tmp$qval < 0.05
colnames(tmp)[3] <- "protein name"
colnames(tmp)[4] <- "Hazard ratio"
colnames(tmp)[5] <- "ln(Hazard ratio)"
colnames(tmp)[6] <- "z-value"
colnames(tmp)[7] <- "p-value"
colnames(tmp)[8] <- "lower bound 95% CI"
colnames(tmp)[9] <- "upper bound 95% CI"
colnames(tmp)[10] <- "q-value"
colnames(tmp)

# Add numeric formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "NUMBER"), rows = 2:(nrow(tmp)+1), cols = c(4:6, 8:9), sheet = "Supp. Table 5", gridExpand = TRUE)

# Add scientific formatting to Spreadsheet
addStyle(wb = wb, style = createStyle(numFmt = "SCIENTIFIC"), rows = 2:(nrow(tmp)+1), cols = c(7, 10), sheet = "Supp. Table 5", gridExpand = TRUE)

# Add column widths to Spreadsheet
setColWidths(wb, sheet ="Supp. Table 5", cols = 1:ncol(tmp), widths = c(rep(10.78, 2), 40.11, 10.78, 13.33, rep(10.78, 6)))

# Add conditional formatting
conditionalFormatting(wb, sheet = "Supp. Table 5", cols = 4, rows = 2:(nrow(tmp)+1), rule = "<1", style = negStyle)
conditionalFormatting(wb, sheet = "Supp. Table 5", cols = 4, rows = 2:(nrow(tmp)+1), rule = ">1", style = posStyle)

openxlsx::writeData(wb, "Supp. Table 5", tmp, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Supplementary_Table_5.xlsx"))


### 6. Supp. Table 6: Feature-reduced 1st-generation and mortality-based models ###

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

### B. Model parameters of the feature-reduced 1st-generation models ###

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

### D. Model parameters of the feature-reduced mortality-based models ###

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
openxlsx::writeData(wb, "Caption", "Supplementary Table 6. A. Model coefficients for the 5 folds of the feature-reduced 1st-generation conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. B. Model parameters for the 5 folds of the feature-reduced 1st-generation conventional and organ-specific models. C. Model coefficients for the 5 folds of the feature-reduced mortality-based conventional and organ-specific models. Proteins that are not specific for the corresponding organ have empty values. D. Model parameters for the 5 folds of the feature-reduced mortality-based conventional and organ-specific models.") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Supp. Table 6A")
openxlsx::writeData(wb, "Supp. Table 6A", tmp.A, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Supp. Table 6B")
openxlsx::writeData(wb, "Supp. Table 6B", tmp.B, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Supp. Table 6C")
openxlsx::writeData(wb, "Supp. Table 6C", tmp.C, headerStyle = header_st)

openxlsx::addWorksheet(wb, "Supp. Table 6D")
openxlsx::writeData(wb, "Supp. Table 6D", tmp.D, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Supplementary_Table_6.xlsx"))

### 7. Supp. Table 7: 1st-generation CSF model ###
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
openxlsx::writeData(wb, "Caption", "Supplementary Table 7. Model coefficients for the 1st-generation conventional and organ-specific models trained on 1,031 proteins measured in the CSF in the Dammer et al. (2022) dataset. Proteins that were not reliably measured in the CSF or are not specific for the corresponding organ have empty values.") # , headerStyle = header_st

header_st <- createStyle(textDecoration = "Bold")

openxlsx::addWorksheet(wb, "Supp. Table 7")
openxlsx::writeData(wb, "Supp. Table 7", tmp, headerStyle = header_st)

openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Supplementary_Table_7.xlsx"))

