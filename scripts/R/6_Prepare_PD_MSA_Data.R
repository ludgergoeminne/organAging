
### The Dammer et al. (2022) Parkinson's Disease / MSA dataset ###

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
coefficients.longitudinal <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal.rds"))

### This comes from 2_Prepare_Longitudinal_data.R ###
longitudinal.proteins <- readRDS(file = paste0(rds.dir, "longitudinal_proteins.rds"))
# sds.longitudinal <- readRDS(file = "standard_deviations_longitudinal.rds")

load(paste0(input.Parkinson.dir, "13g2_Olink-Final_cleanDats_aligned_traits.RData"))

PD_olink_data_list_raw <- vector(mode = "list", length = 2)
names(PD_olink_data_list_raw) <- c("plasma", "CSF")

PD_olink_data_list_imputed <- vector(mode = "list", length = 2)
names(PD_olink_data_list_imputed) <- c("plasma", "CSF")

for(m in 1:length(PD_olink_data_list_imputed)){
  
  if(names(PD_olink_data_list_imputed)[m] == "plasma"){
    ### cleanDat.Olink: no need for log2:
    cleanDat.Olink <- cleanDat.Olink.plasmaNeat # abun.Olink.plasmaNeat
    numericMeta <- numericMeta.plasma
  } else if(names(PD_olink_data_list_imputed)[m] == "CSF"){
    ### cleanDat.Olink: no need for log2:
    cleanDat.Olink <- cleanDat.Olink.CSFneat # abun.final.Olink.CSFneat
    numericMeta <- numericMeta.CSF
  }
  
  cleanDat.Olink <- t(cleanDat.Olink)
  cleanDat.Olink <- cleanDat.Olink[rownames(cleanDat.Olink) %in% numericMeta$GUID, ]
  colnames(cleanDat.Olink) <- gsub("(.+?)\\|.+", "\\1", colnames(cleanDat.Olink))
  
  # No need to rescale with the standard deviations, this is Olink data!
  # cleanDat.Olink <- apply(cleanDat.Olink, 2, function(x){return((x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))})
  # cleanDat.Olink <- sweep(as.matrix(cleanDat.Olink[, names(sds.longitudinal)[names(sds.longitudinal) %in% colnames(cleanDat.Olink)]]), MARGIN=2, sds.longitudinal[names(sds.longitudinal) %in% colnames(cleanDat.Olink)], `*`)
  
  # Only filter out bad participants, not proteins, as we don't want to lose features (we will only use PD_olink_data_list_imputed to test existing models)
  cleanDat.Olink <- filter_olink_data(cleanDat.Olink, participant_cutoff = 0.49948682860075, protein_cutoff = 1)
  cleanDat.Olink <- cleanDat.Olink[, sort(colnames(cleanDat.Olink))]
  
  ### Needed for building the chronological CSF models, here we do filter on protein cut off
  PD_olink_data_list_raw[[m]] <- filter_olink_data(cleanDat.Olink, participant_cutoff = 0.49948682860075, protein_cutoff = 0.1)
  
  cleanDat.Olink <- cleanDat.Olink[, colnames(cleanDat.Olink) %in% longitudinal.proteins]
  PD_olink_data_list_imputed[[m]] <- impute::impute.knn(as.matrix(cleanDat.Olink), k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed = 1714933057)$data
}

# saveRDS(PD_olink_data_list_raw, file = paste0(rds.dir, "PD_olink_data_list_raw.rds"))
PD_olink_data_list_raw <- readRDS(file = paste0(rds.dir, "PD_olink_data_list_raw.rds"))

# saveRDS(PD_olink_data_list_imputed, file = paste0(rds.dir, "cleanDat_Olink_imputed_list.rds"))
PD_olink_data_list_imputed <- readRDS(file = paste0(rds.dir, "cleanDat_Olink_imputed_list.rds"))

numericMeta.list <- vector(mode = "list", length = 2)
names(numericMeta.list) <- c("plasma", "CSF")

for(m in 1:length(numericMeta.list)){
  if(names(numericMeta.list)[m] == "plasma"){
    numericMeta.list[[m]] <- numericMeta.plasma
  } else if(names(numericMeta.list)[m] == "CSF"){
    numericMeta.list[[m]] <- numericMeta.CSF
  }
}

PD_annotation_list <- vector(mode = "list", length = 2)
names(PD_annotation_list) <- c("plasma", "CSF")

for(m in 1:length(PD_annotation_list)){

  PD_annotation_list[[m]] <- data.frame(
    GUID = rownames(PD_olink_data_list_imputed[[m]])
  )
  
  PD_annotation_list[[m]] <- dplyr::inner_join(numericMeta.list[[m]], PD_annotation_list[[m]], by = c("GUID" = "GUID"))
  
  # Add the age at that specific visit
  PD_annotation_list[[m]]$age_at_visit <- PD_annotation_list[[m]]$age_at_baseline + PD_annotation_list[[m]]$visit_month/12
  
  # Add the diagnosis: "Control", "Parkinson's Disease", "Essential Tremor", "Multiple System Atrophy"
  PD_annotation_list[[m]]$diagnosis <- PD_annotation_list[[m]]$diagnosis_latest
  PD_annotation_list[[m]]$diagnosis[PD_annotation_list[[m]]$diagnosis_latest %in% c("Idiopathic PD", "Parkinson's Disease")] <- "Parkinson's Disease"
  PD_annotation_list[[m]]$diagnosis[PD_annotation_list[[m]]$diagnosis_latest == "No PD Nor Other Neurological Disorder"] <- "Control"
}

rm(numericMeta.list)

# saveRDS(PD_annotation_list, file = paste0(rds.dir, "PD_annotation_list.rds"))
PD_annotation_list <- readRDS(file = paste0(rds.dir, "PD_annotation_list.rds"))

### Sanity checks ###
all(rownames(PD_olink_data_list_imputed[[1]]) == PD_annotation_list[[1]]$GUID)
# TRUE
all(rownames(PD_olink_data_list_imputed[[2]]) == PD_annotation_list[[2]]$GUID)
# TRUE
all(rownames(PD_olink_data_list_raw[[1]]) == PD_annotation_list[[1]]$GUID)
# TRUE
all(rownames(PD_olink_data_list_raw[[2]]) == PD_annotation_list[[2]]$GUID)
# TRUE

predicted.ages.PD <- vector(mode = "list", length = 2)
names(predicted.ages.PD) <- c("plasma", "CSF")

predicted.ages.PD <- lapply(predicted.ages.PD, function(x){
  x <- vector(mode = "list", length = 2)
  names(x) <- c("gen1", "gen2")
  x <- lapply(x, function(y){
    y <- vector(mode = "list", length = 2)
    names(y) <- c("predicted", "residual")
    y <- lapply(y, function(z){
      z <- vector(mode = "list", length = length(organ.proteins))
      names(z) <- names(organ.proteins)
      return(z)
    })
    return(y)
  })
  return(x)
})

coefficients.longitudinal.Parkinson <- vector(mode = "list", length = 2)
names(coefficients.longitudinal.Parkinson) <- c("plasma", "CSF")

coefficients.longitudinal.Parkinson <- lapply(coefficients.longitudinal.Parkinson, function(x){
  x <- vector(mode = "list", length = 2)
  names(x) <- c("gen1", "gen2")
  x <- lapply(x, function(z){
      z <- vector(mode = "list", length = length(organ.proteins))
      names(z) <- names(organ.proteins)
      return(z)
    })
  return(x)
})


for(m in 1:length(PD_annotation_list)){
  for(g in 1:length(coefficients.longitudinal)){
    for(k in 1:length(coefficients.longitudinal[[g]])){
      
      if(!(names(coefficients.longitudinal[[g]])[k] %in% c("Bladder", "Thyroid"))){
        
        if(names(coefficients.longitudinal)[g] == "gen1"){
 
          coefficients <- coefficients.longitudinal[[g]][[k]]
          coefficients.longitudinal.Parkinson[[m]][[g]][[k]] <- c(coefficients["Intercept"], coefficients[names(coefficients) %in% colnames(PD_olink_data_list_imputed[[m]]), drop = FALSE])
          coefficients <- coefficients.longitudinal.Parkinson[[m]][[g]][[k]][coefficients.longitudinal.Parkinson[[m]][[g]][[k]] != 0]
          tmp2 <- PD_olink_data_list_imputed[[m]][, names(coefficients)[-1], drop = FALSE]
          
          # Add the residual ages
          predicted.ages.PD[[m]][[g]][["predicted"]][[k]] <- coefficients["Intercept"] + rowSums(sweep(tmp2, 2, coefficients[-1], FUN="*"))
          predicted.ages.PD[[m]][[g]][["residual"]][[k]] <- resid(lm(predicted.ages.PD[[m]][[g]][["predicted"]][[k]] ~ PD_annotation_list[[m]][, "age_at_visit"]))
          
        } else if(names(coefficients.longitudinal)[g] == "gen2"){
          
          coefficients <- coefficients.longitudinal[[g]][[k]]
          coefficients.longitudinal.Parkinson[[m]][[g]][[k]] <- coefficients[names(coefficients) %in% colnames(PD_olink_data_list_imputed[[m]]), drop = FALSE]
          coefficients <- coefficients.longitudinal.Parkinson[[m]][[g]][[k]][coefficients.longitudinal.Parkinson[[m]][[g]][[k]] != 0]
          tmp2 <- PD_olink_data_list_imputed[[m]][, names(coefficients), drop = FALSE]
          
          # Add the residual ages
          predicted.ages.PD[[m]][[g]][["predicted"]][[k]] <- rowSums(sweep(tmp2, 2, coefficients, FUN="*"))
          predicted.ages.PD[[m]][[g]][["residual"]][[k]] <- resid(lm(predicted.ages.PD[[m]][[g]][["predicted"]][[k]] ~ PD_annotation_list[[m]][, "age_at_visit"]))
          
        }
        
      }
    }
  }
}



# saveRDS(predicted.ages.PD, file = paste0(rds.dir, "predicted_ages_PD.rds"))
predicted.ages.PD <- readRDS(file = paste0(rds.dir, "predicted_ages_PD.rds"))

# saveRDS(coefficients.longitudinal.Parkinson, file = paste0(rds.dir, "coefficients_longitudinal_Parkinson.rds"))
coefficients.longitudinal.Parkinson <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal_Parkinson.rds"))



