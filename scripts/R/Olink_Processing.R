get_olink_coding <- function(protein_file){
  # Import the coding of the Olink proteins; we will use this for sorting
  coding143 <- read.table(protein_file, quote = "", comment.char = "", sep = "\t", header = TRUE)
  coding143$coding <- as.character(coding143$coding)
  coding143$protein <- unlist(lapply(strsplit(coding143$meaning, ";"), function(x){return(x[1])}))
  coding143$protein.name <- unlist(lapply(strsplit(coding143$meaning, ";"), function(x){return(x[2])}))
  coding143 <- coding143[, c("coding", "protein", "protein.name")]
  return(coding143)
}

load_olink_data <- function(olink_file, protein_file, instance=0){
  olink <- data.table::fread(olink_file, header=TRUE, sep="\t")
  olink <- as.data.frame(olink)
  
  coding143 <- get_olink_coding(protein_file)
  olink.wide <- tidyr::pivot_wider(olink[olink$ins_index == instance, -2], names_from = protein_id, values_from = result)
  
  olink.wide <- olink.wide[, c("eid", coding143$coding[coding143$coding %in% colnames(olink.wide)])]
  olink.wide <- as.data.frame(olink.wide)
  rownames(olink.wide) <- olink.wide$eid
  olink.wide <- olink.wide[,-1]
  colnames(olink.wide) <- coding143[match(colnames(olink.wide), coding143$coding), "protein"]
  return(olink.wide)
}

filter_olink_data <- function(olink.wide, participant_cutoff = 0.49948682860075, protein_cutoff = 0.1){
  olink.wide <- olink.wide[which(rowMeans(is.na(olink.wide)) < participant_cutoff), ] 
  olink.wide <- olink.wide[, which(colMeans(is.na(olink.wide)) < protein_cutoff)]
  return(olink.wide)
}

impute_olink <- function(olink.wide){
  rownames_olink <- rownames(olink.wide)
  olink.values <- as.matrix(olink.wide)
  values.imputed <- impute::impute.knn(olink.values, k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed = 1714933057)$data
  
  df_imputed <- as.data.frame(values.imputed)
  rownames(df_imputed) <- rownames_olink
  return(df_imputed)
}

add_outcome_columns_to_olink <- function(olink_data, bd.olink, visit = "first"){
  
  olink_data[["Age"]] <- bd.olink[[paste0("age_", visit, "_visit")]]
  olink_data[["Status"]] <- bd.olink[["status_mortality"]]
  olink_data[["Time"]] <- bd.olink[["time_mortality"]]
  
  return(olink_data)
}

kfold_split_save <- function(data, save_dir, bd, n_splits = 5, visit = "first", seed = 5049857){
  set.seed(seed)
  index <- 1:nrow(data)
  subset_size <- length(index) %/% n_splits
  shuffled_index <- sample(index)
  
  subsets <- list()
  
  # Create the subsets
  for (i in 1:(n_splits-1)) {
    subsets[[i]] = shuffled_index[((i-1)*subset_size+1): (i*subset_size)]
  }
  subsets[[n_splits]] = shuffled_index[((n_splits-1)*subset_size+1):length(shuffled_index)]
  
  for (i in 1:n_splits){
    subset_index <- sort(subsets[[i]])
    train_data <- data[-subset_index,]
    test_data <- data[subset_index,]
    
    train_data <- impute_olink(train_data)
    test_data <- impute_olink(test_data)
    
    train_data <- add_outcome_columns_to_olink(train_data, bd[-subset_index,], visit)
    test_data <- add_outcome_columns_to_olink(test_data, bd[subset_index,], visit)
    
    write.csv(train_data, paste0(save_dir, "train_imputed_", i, ".csv"))
    write.csv(test_data, paste0(save_dir, "test_imputed_", i, ".csv"))
  }
  
}