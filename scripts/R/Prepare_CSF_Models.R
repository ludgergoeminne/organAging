
### This comes from 6_Prepare_Parkinson_MSA_Data.R
PD_olink_data_list_raw <- readRDS(file = paste0(rds.dir, "PD_olink_data_list_raw.rds"))
PD_annotation_list <- readRDS(file = paste0(rds.dir, "PD_annotation_list.rds"))

### This comes from 3_Process_GTEx.R
subsets.Olink <- readRDS(file = paste0(rds.dir, "subsets_Olink.rds"))

### Sanity checks ###
all(rownames(PD_olink_data_list_raw[[1]]) == PD_annotation_list[[1]]$GUID)
# TRUE
all(rownames(PD_olink_data_list_raw[[2]]) == PD_annotation_list[[2]]$GUID)
# TRUE

# Train CSF model only for people at the first visit #

CSF.annotation <- PD_annotation_list[["CSF"]][which(PD_annotation_list[["CSF"]]$visit_month == 0),]
tmp <- PD_olink_data_list_raw[["CSF"]][which(PD_annotation_list[["CSF"]]$visit_month == 0),]

# Get dimensions
dim(tmp)
# 184 1031

# Add the "Age" column
CSF.Olink <- cbind(tmp, Age = CSF.annotation[["age_at_visit"]])

### Sanity check ###
any(is.na(CSF.Olink[, "Age"]))
# FALSE

# Export the results so that they can be loaded into Python
write.csv(CSF.Olink, paste0(dir.Python.input.CSF, "full_CSF.csv"))

### 2. Get the set of organ-enriched proteins per organ for the chronological models trained on the CSF data from Dammer et al. (2022) ###

proteins <- coding143[match(colnames(PD_olink_data_list_raw[["CSF"]]), coding143$protein), "protein"]
organ.proteins.CSF <- lapply(subsets.Olink, function(x){
  x <- sort(x[x %in% coding143$protein]) # Sort alphabetically; otherwise, we have the random order from GTEx
  x <- x[x %in% proteins]
  return(x)
})
organ.proteins.CSF$Organismal <- proteins[!(proteins %in% unlist(subsets.Olink)[unlist(subsets.Olink) %in% coding143$protein])]
organ.proteins.CSF$`Multi-organ` <- proteins[proteins %in% unlist(subsets.Olink)[unlist(subsets.Olink) %in% coding143$protein]]
organ.proteins.CSF$Conventional <- proteins

# saveRDS(organ.proteins.CSF, paste0(rds.dir, "organ_proteins_CSF.rds"))
organ.proteins.CSF <- readRDS(paste0(rds.dir, "organ_proteins_CSF.rds"))

ListJSON.CSF <- toJSON(organ.proteins.CSF, pretty = TRUE, auto_unbox = TRUE)
write(ListJSON.CSF, paste0(dir.4x.FC.genes, "GTEx_4x_FC_genes_CSF.json"))

# After this, run the Python scripts
# Where you will find the output:
# dir.gen1.models.CSF

