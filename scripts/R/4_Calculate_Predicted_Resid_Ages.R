
### Basic directories ###
# rds.dir
# dir.gen1.models
# dir.gen2.models

### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from 2_Prepare_Longitudinal_Data.R
olink_data_list <- readRDS(file = paste0(rds.dir, "olink_data_list.rds"))
olink_bd_annotation_list <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_list.rds"))

### 1. Calculate for the non-longitudinal models ###

### Create a list with 1st- and 2nd-generation age and residual predictions for all organs in training and test datasets ###
predicted.ages.1 <- vector(mode = "list", length = 2)
names(predicted.ages.1) <- c("gen1", "gen2")
predicted.ages.1 <- lapply(predicted.ages.1, function(x){
  x <- vector(mode = "list", length = 2)
  names(x) <- c("predicted", "residual")
  x <- lapply(x, function(z){
    z <- vector(mode = "list", length = length(organ.proteins))
    names(z) <- names(organ.proteins)
    return(z)
  })
  return(x)
})

### Create a list with 1st- and 2nd-generation coefficients for all organs ###
coefficients <- vector(mode = "list", length = 2)
names(coefficients) <- c("gen1", "gen2")
coefficients <- lapply(coefficients, function(x){
  x <- vector(mode = "list", length = length(organ.proteins))
  names(x) <- names(organ.proteins)
  return(x)
})

read_csv_fast <- function(...){
  out <- as.data.frame(data.table::fread(...))
  rownames(out) <- out[, 1]
  out <- out[, -1, drop = FALSE]
}

df.train.1 <- read_csv_fast(paste0(dir.Python.input, "train_imputed_1.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.1 <- read_csv_fast(paste0(dir.Python.input, "test_imputed_1.csv"), header = TRUE, sep = ",", quote = "\"")
df.full.1 <- rbind(df.train.1, df.test.1)
df.full.1 <- df.full.1[as.character(olink_bd_annotation_0$eid),]

# Inspect the data
df.train.1[1:10, 1:10]
df.test.1[1:10, 1:10]

training.eids <- rownames(df.train.1)
test.eids <- rownames(df.test.1)

# saveRDS(training.eids, file = paste0(rds.dir, "training_eids.rds"))
training.eids <- readRDS(file = paste0(rds.dir, "training_eids.rds"))
# saveRDS(test.eids, file = paste0(rds.dir, "test_eids.rds"))
test.eids <- readRDS(file = paste0(rds.dir, "test_eids.rds"))

### Generation 1 ###

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder"))){
    coefficients[["gen1"]][[k]] <- read.csv(paste0(dir.gen1.models, names(organ.proteins)[k], "_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    coefficients[["gen1"]][[k]] <- unlist(coefficients[["gen1"]][[k]][1, , drop = FALSE])
    predicted.ages.1[["gen1"]][["predicted"]][[k]] <- coefficients[["gen1"]][[k]]["Intercept"]+rowSums(sweep(as.matrix(df.full.1[, names(coefficients[["gen1"]][[k]])[-1], drop = FALSE]), MARGIN=2, coefficients[["gen1"]][[k]][-1], `*`))
    
    predicted.ages.1[["gen1"]][["residual"]][[k]] <- c(resid(lm(predicted.ages.1[["gen1"]][["predicted"]][[k]][training.eids]~olink_bd_annotation_0[training.eids, "age_first_visit"])),
                                                       resid(lm(predicted.ages.1[["gen1"]][["predicted"]][[k]][test.eids]~olink_bd_annotation_0[test.eids, "age_first_visit"])))[as.character(olink_bd_annotation_0$eid)]
  }
}

### Generation 2 ###

for(k in 1:length(organ.proteins)){
  if(names(organ.proteins)[k] != "Bladder"){
    coefficients[["gen2"]][[k]] <- read.csv(paste0(dir.gen2.models, names(organ.proteins)[k], "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    coefficients[["gen2"]][[k]] <- unlist(coefficients[["gen2"]][[k]][1, , drop = FALSE])
    predicted.ages.1[["gen2"]][["predicted"]][[k]] <- rowSums(sweep(as.matrix(df.full.1[, names(coefficients[["gen2"]][[k]]), drop = FALSE]), MARGIN=2, coefficients[["gen2"]][[k]], `*`))
    
    predicted.ages.1[["gen2"]][["residual"]][[k]] <- c(resid(lm(predicted.ages.1[["gen2"]][["predicted"]][[k]][training.eids]~olink_bd_annotation_0[training.eids, "age_first_visit"])), 
                                                       resid(lm(predicted.ages.1[["gen2"]][["predicted"]][[k]][test.eids]~olink_bd_annotation_0[test.eids, "age_first_visit"])))[as.character(olink_bd_annotation_0$eid)]
  }
}

# saveRDS(coefficients, file = paste0(rds.dir, "coefficients.rds"))
coefficients <- readRDS(file = paste0(rds.dir, "coefficients.rds"))

# saveRDS(predicted.ages.1, file = paste0(rds.dir, "predicted_ages_1.rds"))
predicted.ages.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_1.rds"))

### 2. Calculate out-of-fold predictions for the non-longitudinal models ###

predicted.ages.oof <- vector(mode = "list", length = 2)
names(predicted.ages.oof) <- c("gen1", "gen2")
predicted.ages.oof <- lapply(predicted.ages.oof, function(x){
  x <- vector(mode = "list", length = 2)
  names(x) <- c("predicted", "residual")
  x <- lapply(x, function(z){
    z <- vector(mode = "list", length = length(organ.proteins))
    names(z) <- names(organ.proteins)
    return(z)
  })
  return(x)
})

df.test.1 <- read_csv_fast(paste0(dir.Python.input, "test_imputed_1.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.2 <- read_csv_fast(paste0(dir.Python.input, "test_imputed_2.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.3 <- read_csv_fast(paste0(dir.Python.input, "test_imputed_3.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.4 <- read_csv_fast(paste0(dir.Python.input, "test_imputed_4.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.5 <- read_csv_fast(paste0(dir.Python.input, "test_imputed_5.csv"), header = TRUE, sep = ",", quote = "\"")

df.test.list <- list(df.test.1,
                     df.test.2,
                     df.test.3,
                     df.test.4,
                     df.test.5)

# saveRDS(df.test.list, file = paste0(rds.dir, "df_test_list.rds"))
df.test.list <- readRDS(file = paste0(rds.dir, "df_test_list.rds"))

### Generation 1 ###

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder"))){
    coef.tmp <- read.csv(paste0(dir.gen1.models, names(organ.proteins)[k], "_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    
    tmp <- vector(mode = "list", length = 5)
    for(j in 1:5){
      tmp[[j]] <- coef.tmp[j, "Intercept"]+rowSums(sweep(as.matrix(df.test.list[[j]][, colnames(coef.tmp)[-1], drop = FALSE]), MARGIN=2, unlist(coef.tmp[j, -1, drop = FALSE]), `*`))
    }
    predicted.ages.oof[["gen1"]][["predicted"]][[k]] <- do.call("c", tmp)
    predicted.ages.oof[["gen1"]][["predicted"]][[k]] <- predicted.ages.oof[["gen1"]][["predicted"]][[k]][as.character(olink_bd_annotation_0$eid)]
    
    predicted.ages.oof[["gen1"]][["residual"]][[k]] <- resid(lm(predicted.ages.oof[["gen1"]][["predicted"]][[k]]~olink_bd_annotation_0[, "age_first_visit"]))
  }
}

### Generation 2 ###

for(k in 1:length(organ.proteins)){
  if(names(organ.proteins)[k] != "Bladder"){
    coef.tmp <- read.csv(paste0(dir.gen2.models, names(organ.proteins)[k], "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    
    tmp <- vector(mode = "list", length = 5)
    for(j in 1:5){
      tmp[[j]] <- rowSums(sweep(as.matrix(df.test.list[[j]][, colnames(coef.tmp), drop = FALSE]), MARGIN=2, unlist(coef.tmp[j, , drop = FALSE]), `*`))
    }
    
    predicted.ages.oof[["gen2"]][["predicted"]][[k]] <- do.call("c", tmp)
    predicted.ages.oof[["gen2"]][["predicted"]][[k]] <- predicted.ages.oof[["gen2"]][["predicted"]][[k]][as.character(olink_bd_annotation_0$eid)]
    predicted.ages.oof[["gen2"]][["residual"]][[k]] <- resid(lm(predicted.ages.oof[["gen2"]][["predicted"]][[k]]~olink_bd_annotation_0[, "age_first_visit"]))
  }
}

# saveRDS(predicted.ages.oof, file = paste0(rds.dir, "predicted_ages_oof.rds"))
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

### 3. Calculate predictions in the first fold for the feature-reduced ("longitudinal") models ###
# I.e. the models re-trained on the proteins that are available longitudinally

df.test.1.longitudinal <- read_csv_fast(paste0(dir.Python.input.longitudinal, "test_imputed_1.csv"), header = TRUE, sep = ",", quote = "\"")
df.train.1.longitudinal <- read_csv_fast(paste0(dir.Python.input.longitudinal, "train_imputed_1.csv"), header = TRUE, sep = ",", quote = "\"")
df.full.1.longitudinal.i0 <- rbind(df.train.1.longitudinal, df.test.1.longitudinal)
df.full.1.longitudinal.i0 <- df.full.1.longitudinal.i0[as.character(olink_bd_annotation_list[["instance_0"]]$eid),]

longitudinal.training.eids <- rownames(df.train.1.longitudinal)
longitudinal.test.eids <- rownames(df.test.1.longitudinal)

# saveRDS(longitudinal.training.eids, file = paste0(rds.dir, "longitudinal_training_eids.rds"))
longitudinal.training.eids <- readRDS(file = paste0(rds.dir, "longitudinal_training_eids.rds"))
# saveRDS(longitudinal.test.eids, file = paste0(rds.dir, "longitudinal_test_eids.rds"))
longitudinal.test.eids <- readRDS(file = paste0(rds.dir, "longitudinal_test_eids.rds"))

set.seed(1921157651)
df.full.longitudinal.i2 <- impute_olink(olink_data_list[["instance_2"]])
df.full.longitudinal.i2 <- add_outcome_columns_to_olink(df.full.longitudinal.i2, olink_bd_annotation_list[["instance_2"]], "third")

set.seed(1921157683)
df.full.longitudinal.i3 <- impute_olink(olink_data_list[["instance_3"]])
df.full.longitudinal.i3 <- add_outcome_columns_to_olink(df.full.longitudinal.i3, olink_bd_annotation_list[["instance_3"]], "fourth")

df.full.longitudinal.list <- list(df.full.1.longitudinal.i0,
                                  df.full.longitudinal.i2,
                                  df.full.longitudinal.i3)

predicted.ages.longitudinal.1 <- vector(mode = "list", length = length(olink_data_list))
names(predicted.ages.longitudinal.1) <- names(olink_data_list)

predicted.ages.longitudinal.1 <- lapply(predicted.ages.longitudinal.1, function(x){
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

### Create a list with 1st- and 2nd-generation coefficients for all organs ###
coefficients.longitudinal <- vector(mode = "list", length = 2)
names(coefficients.longitudinal) <- c("gen1", "gen2")
coefficients.longitudinal <- lapply(coefficients.longitudinal, function(x){
  x <- vector(mode = "list", length = length(organ.proteins))
  names(x) <- names(organ.proteins)
  return(x)
})

### Age visit, training, test ids!!!

visit.names <- c("first", "third", "fourth")

### Generation 1 ###

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder", "Thyroid"))){
    coefficients.longitudinal[["gen1"]][[k]] <- read.csv(paste0(dir.gen1.models.longitudinal, names(organ.proteins)[k], "_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
    coefficients.longitudinal[["gen1"]][[k]] <- unlist(coefficients.longitudinal[["gen1"]][[k]][1, , drop = FALSE])
    
    for(i in 1:length(predicted.ages.longitudinal.1)){
      predicted.ages.longitudinal.1[[i]][["gen1"]][["predicted"]][[k]] <- coefficients.longitudinal[["gen1"]][[k]]["Intercept"]+rowSums(sweep(as.matrix(df.full.longitudinal.list[[i]][, names(coefficients.longitudinal[["gen1"]][[k]])[-1], drop = FALSE]), MARGIN=2, coefficients.longitudinal[["gen1"]][[k]][-1], `*`))
      
      if(i == 1){
        predicted.ages.longitudinal.1[[i]][["gen1"]][["residual"]][[k]] <- c(resid(lm(predicted.ages.longitudinal.1[[i]][["gen1"]][["predicted"]][[k]][longitudinal.training.eids]~olink_bd_annotation_list[[i]][longitudinal.training.eids, paste0("age_", visit.names[i], "_visit")])),
                                                                             resid(lm(predicted.ages.longitudinal.1[[i]][["gen1"]][["predicted"]][[k]][longitudinal.test.eids]~olink_bd_annotation_list[[i]][longitudinal.test.eids, paste0("age_", visit.names[i], "_visit")])))[as.character(olink_bd_annotation_list[[i]]$eid)]
      } else if(i %in% c(2, 3)){
        predicted.ages.longitudinal.1[[i]][["gen1"]][["residual"]][[k]] <- resid(lm(predicted.ages.longitudinal.1[[i]][["gen1"]][["predicted"]][[k]]~olink_bd_annotation_list[[i]][, paste0("age_", visit.names[i], "_visit")]))
      }
    }
  }
}

### Generation 2 ###

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder", "Thyroid"))){
    coefficients.longitudinal[["gen2"]][[k]] <- read.csv(paste0(dir.gen2.models.longitudinal, names(organ.proteins)[k], "_mortality_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
    coefficients.longitudinal[["gen2"]][[k]] <- unlist(coefficients.longitudinal[["gen2"]][[k]][1, , drop = FALSE])
    
    for(i in 1:length(predicted.ages.longitudinal.1)){
      predicted.ages.longitudinal.1[[i]][["gen2"]][["predicted"]][[k]] <- rowSums(sweep(as.matrix(df.full.longitudinal.list[[i]][, names(coefficients.longitudinal[["gen2"]][[k]]), drop = FALSE]), MARGIN=2, coefficients.longitudinal[["gen2"]][[k]], `*`))
      
      if(i == 1){
        predicted.ages.longitudinal.1[[i]][["gen2"]][["residual"]][[k]] <- c(resid(lm(predicted.ages.longitudinal.1[[i]][["gen2"]][["predicted"]][[k]][longitudinal.training.eids]~olink_bd_annotation_list[[i]][longitudinal.training.eids, paste0("age_", visit.names[i], "_visit")])), 
                                                                             resid(lm(predicted.ages.longitudinal.1[[i]][["gen2"]][["predicted"]][[k]][longitudinal.test.eids]~olink_bd_annotation_list[[i]][longitudinal.test.eids, paste0("age_", visit.names[i], "_visit")])))[as.character(olink_bd_annotation_list[[i]]$eid)]
      } else if(i %in% c(2, 3)){
        predicted.ages.longitudinal.1[[i]][["gen2"]][["residual"]][[k]] <- resid(lm(predicted.ages.longitudinal.1[[i]][["gen2"]][["predicted"]][[k]]~olink_bd_annotation_list[[i]][, paste0("age_", visit.names[i], "_visit")]))
      }
    }
  }
}

# saveRDS(coefficients.longitudinal, file = paste0(rds.dir, "coefficients_longitudinal.rds"))
coefficients.longitudinal <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal.rds"))

# saveRDS(predicted.ages.longitudinal.1, file = paste0(rds.dir, "predicted_ages_longitudinal_1.rds"))
predicted.ages.longitudinal.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_longitudinal_1.rds"))


### 4. Calculate out-of-fold predictions for the longitudinal models at the first instance ###
# I.e. the models re-trained on the proteins that are available longitudinally

predicted.ages.longitudinal.oof <- vector(mode = "list", length = 2)
names(predicted.ages.longitudinal.oof) <- c("gen1", "gen2")

predicted.ages.longitudinal.oof <- lapply(predicted.ages.longitudinal.oof, function(y){
  y <- vector(mode = "list", length = 2)
  names(y) <- c("predicted", "residual")
  y <- lapply(y, function(z){
    z <- vector(mode = "list", length = length(organ.proteins))
    names(z) <- names(organ.proteins)
    return(z)
  })
  return(y)
})

df.test.1.longitudinal <- read_csv_fast(paste0(dir.Python.input.longitudinal, "test_imputed_1.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.2.longitudinal <- read_csv_fast(paste0(dir.Python.input.longitudinal, "test_imputed_2.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.3.longitudinal <- read_csv_fast(paste0(dir.Python.input.longitudinal, "test_imputed_3.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.4.longitudinal <- read_csv_fast(paste0(dir.Python.input.longitudinal, "test_imputed_4.csv"), header = TRUE, sep = ",", quote = "\"")
df.test.5.longitudinal <- read_csv_fast(paste0(dir.Python.input.longitudinal, "test_imputed_5.csv"), header = TRUE, sep = ",", quote = "\"")

df.test.list.longitudinal <- list(df.test.1.longitudinal,
                                  df.test.2.longitudinal,
                                  df.test.3.longitudinal,
                                  df.test.4.longitudinal,
                                  df.test.5.longitudinal)

# saveRDS(df.test.list.longitudinal, file = paste0(rds.dir, "df_test_list_longitudinal.rds"))
df.test.list.longitudinal <- readRDS(file = paste0(rds.dir, "df_test_list_longitudinal.rds"))

### Generation 1 ###

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder", "Thyroid"))){
    coef.tmp <- read.csv(paste0(dir.gen1.models.longitudinal, names(organ.proteins)[k], "_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
    
    tmp <- vector(mode = "list", length = 5)
    for(j in 1:5){
      tmp[[j]] <- coef.tmp[j, "Intercept"]+rowSums(sweep(as.matrix(df.test.list.longitudinal[[j]][, colnames(coef.tmp)[-1], drop = FALSE]), MARGIN=2, unlist(coef.tmp[j, -1, drop = FALSE]), `*`))
    }
    predicted.ages.longitudinal.oof[["gen1"]][["predicted"]][[k]] <- do.call("c", tmp)
    predicted.ages.longitudinal.oof[["gen1"]][["predicted"]][[k]] <- predicted.ages.longitudinal.oof[["gen1"]][["predicted"]][[k]][as.character(olink_bd_annotation_list[["instance_0"]]$eid)]
    predicted.ages.longitudinal.oof[["gen1"]][["residual"]][[k]] <- resid(lm(predicted.ages.longitudinal.oof[["gen1"]][["predicted"]][[k]]~olink_bd_annotation_list[["instance_0"]][, "age_first_visit"]))
  }
}

### Generation 2 ###

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder", "Thyroid"))){
    coef.tmp <- read.csv(paste0(dir.gen2.models.longitudinal, names(organ.proteins)[k], "_mortality_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
    
    tmp <- vector(mode = "list", length = 5)
    for(j in 1:5){
      tmp[[j]] <- rowSums(sweep(as.matrix(df.test.list.longitudinal[[j]][, colnames(coef.tmp), drop = FALSE]), MARGIN=2, unlist(coef.tmp[j, , drop = FALSE]), `*`))
    }
    
    predicted.ages.longitudinal.oof[["gen2"]][["predicted"]][[k]] <- do.call("c", tmp)
    predicted.ages.longitudinal.oof[["gen2"]][["predicted"]][[k]] <- predicted.ages.longitudinal.oof[["gen2"]][["predicted"]][[k]][as.character(olink_bd_annotation_list[["instance_0"]]$eid)]
    predicted.ages.longitudinal.oof[["gen2"]][["residual"]][[k]] <- resid(lm(predicted.ages.longitudinal.oof[["gen2"]][["predicted"]][[k]]~olink_bd_annotation_list[["instance_0"]][, "age_first_visit"]))
  }
}

# saveRDS(predicted.ages.longitudinal.oof, file = paste0(rds.dir, "predicted_ages_longitudinal_oof.rds"))
predicted.ages.longitudinal.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_longitudinal_oof.rds"))

