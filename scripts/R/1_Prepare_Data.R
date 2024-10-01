source(paste0(scripts.dir, "UKBB_Processing.R")) # function get_ukb_participant_annotation_data
source(paste0(scripts.dir, "Olink_Processing.R"))
source(paste0(scripts.dir, "Outcomes_Processing.R"))

bd <- get_ukb_participant_annotation_data(input.UKBB.dir, ukbb.tab.file)
# NAs introduced by coercion
rownames(bd) <- bd$eid

### Process data at first visit ###

olink_data_unfiltered_0 <- load_olink_data(paste0(input.olink.dir, olink.txt.file), paste0(input.coding.dir, "coding143.tsv"), instance = 0)
olink_data_unfiltered_0 <- olink_data_unfiltered_0[rownames(olink_data_unfiltered_0) %in% bd$eid,]
olink_data_0 <- filter_olink_data(olink_data_unfiltered_0, participant_cutoff = 0.49948682860075, protein_cutoff = 0.1)

included_rownames <- rownames(olink_data_0)
olink_bd_annotation_0 <- bd[included_rownames,]
rm(included_rownames)

### Annotate the data at first visit ###
general_hazard_outcomes <- get_general_outcomes()
olink_bd_annotation_0 <- annotate_outcomes(olink_bd_annotation_0, general_hazard_outcomes, visit = "first")

# Create, impute and save train-test splits for olink data FIRST visit
# Required for K-fold validation (python knn is relatively slow compared to R)
kfold_split_save(data = olink_data_0, save_dir = dir.Python.input, olink_bd_annotation_0, n_splits = 5, visit = "first", seed = 5049857)

# Calculate the standard deviations (useful for rescaling on external data)
sds <- apply(olink_data_0, 2, sd, na.rm = TRUE)

# saveRDS(olink_bd_annotation_0, file = paste0(rds.dir, "olink_bd_annotation_0.rds"))
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

# saveRDS(olink_data_unfiltered_0, file = paste0(rds.dir, "olink_data_unfiltered_0.rds"))
olink_data_unfiltered_0 <- readRDS(file = paste0(rds.dir, "olink_data_unfiltered_0.rds"))

# saveRDS(olink_data_0, file = paste0(rds.dir, "olink_data_0.rds"))
olink_data_0 <- readRDS(file = paste0(rds.dir, "olink_data_0.rds"))

# saveRDS(sds, file = paste0(rds.dir, "standard_deviations.rds"))
sds <- readRDS(file = paste0(rds.dir, "standard_deviations.rds"))

# if you intend to perform GTEx tissue specificity analysis, store coding separately
coding143 <- get_olink_coding(paste0(input.coding.dir, "coding143.tsv"))
# saveRDS(coding143, file = paste0(rds.dir, "coding143.rds"))
coding143 <- readRDS(paste0(rds.dir, "coding143.rds"))

# olink_data_unfiltered_list is only useful for "find_filtering_thresholds.R"
# It can be removed if you do not want to run that script
# rm(olink_data_unfiltered_list)


