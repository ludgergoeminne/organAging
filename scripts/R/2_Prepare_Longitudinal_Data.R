source(paste0(scripts.dir, "UKBB_Processing.R")) # function get_ukb_participant_annotation_data
source(paste0(scripts.dir, "Olink_Processing.R"))
source(paste0(scripts.dir, "Outcomes_Processing.R"))

# bd comes from 1_Prepare_Data.R
# bd <- get_ukb_participant_annotation_data(input.UKBB.dir, ukbb.tab.file)
# # NAs introduced by coercion
# rownames(bd) <- bd$eid

### Process longitudinal data ###

all.longitudinal.proteins <- colnames(load_olink_data(paste0(input.olink.dir, "olink_data.txt"), paste0(input.coding.dir, "coding143.tsv"), instance = 3))

instances <- c(0,2,3)

olink_data_unfiltered_list <- vector(mode = "list", length = length(instances))
names(olink_data_unfiltered_list) <- paste0("instance_", instances)
olink_data_list <- vector(mode = "list", length = length(instances))
names(olink_data_list) <- paste0("instance_", instances)
olink_bd_annotation_list <- vector(mode = "list", length = length(instances))
names(olink_bd_annotation_list) <- paste0("instance_", instances)

for (instance in instances){
  olink_data_unfiltered <- load_olink_data(paste0(input.olink.dir, olink.txt.file), paste0(input.coding.dir, "coding143.tsv"), instance = instance)
  olink_data_unfiltered <- olink_data_unfiltered[rownames(olink_data_unfiltered) %in% bd$eid,]
  olink_data_unfiltered <- olink_data_unfiltered[, colnames(olink_data_unfiltered) %in% all.longitudinal.proteins]
  olink_data <- filter_olink_data(olink_data_unfiltered, participant_cutoff = 0.49948682860075, protein_cutoff = 0.1)
  
  included_rownames <- rownames(olink_data)
  olink_bd_annotation_list[[paste0("instance_", instance)]] <- bd[included_rownames,]
  
  olink_data_unfiltered_list[[paste0("instance_", instance)]] <- olink_data_unfiltered
  olink_data_list[[paste0("instance_", instance)]] <- olink_data
}
rm(included_rownames)

longitudinal.proteins <- names(which(table(c(colnames(olink_data_list[["instance_0"]]), colnames(olink_data_list[["instance_2"]]), colnames(olink_data_list[["instance_3"]]))) == 3))
longitudinal.eids <- names(which(table(c(rownames(olink_data_list[["instance_0"]]), rownames(olink_data_list[["instance_2"]]), rownames(olink_data_list[["instance_3"]]))) == 3))

length(longitudinal.proteins)
# 1451
length(longitudinal.eids)
# 985

# saveRDS(longitudinal.proteins, file = paste0(rds.dir, "longitudinal_proteins.rds"))
longitudinal.proteins <- readRDS(file = paste0(rds.dir, "longitudinal_proteins.rds"))

# saveRDS(longitudinal.eids, file = paste0(rds.dir, "longitudinal_eids.rds"))
longitudinal.eids <- readRDS(file = paste0(rds.dir, "longitudinal_eids.rds"))

# Filter based on longitudinal proteins: only <10% missing values *in each instance* should be retained, 
# Do not filter yet on longitudinal eids
for (instance in instances){
  olink_data_list[[paste0("instance_", instance)]] <- olink_data_list[[paste0("instance_", instance)]][, longitudinal.proteins]
}

### Annotate the longitudinal data ###
general_hazard_outcomes <- get_general_outcomes()
for (i in 1:length(instances)){
  olink_bd_annotation_list[[paste0("instance_", instances[i])]] <- annotate_outcomes(olink_bd_annotation_list[[paste0("instance_", instances[i])]], 
                                                                                 general_hazard_outcomes, c("first", "third", "fourth")[i])
}

# Calculate the standard deviations (useful for rescaling on external data)
sds.longitudinal <- apply(olink_data_list[[paste0("instance_0")]], 2, sd, na.rm = TRUE)

# saveRDS(olink_bd_annotation_list, file = paste0(rds.dir, "olink_bd_annotation_list.rds"))
olink_bd_annotation_list <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_list.rds"))

# saveRDS(olink_data_list, file = paste0(rds.dir, "olink_data_list.rds"))
olink_data_list <- readRDS(file = paste0(rds.dir, "olink_data_list.rds"))

# saveRDS(sds.longitudinal, file = paste0(rds.dir, "standard_deviations_longitudinal.rds"))
sds.longitudinal <- readRDS(file = paste0(rds.dir, "standard_deviations_longitudinal.rds"))

# Create, impute and save train-test splits for longitudinal Olink data FIRST visit
# Required for K-fold validation (python knn is relatively slow compared to R)
kfold_split_save(data = olink_data_list[["instance_0"]], save_dir = dir.Python.input.longitudinal, olink_bd_annotation_list[["instance_0"]], n_splits = 5, visit = "first", seed = 5049857)


