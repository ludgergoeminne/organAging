get_general_outcomes <- function(){
  cox.outcomes <- list(
    c("p130836",
      "p130838",
      "p130840",
      "p130842"),
    
    c("p130708"),
    
    c("p131368"), 
    
    c("p131298", 
      "p131300"),
    
    c("p131492"), 
    
    c("p131354"),
    
    c("p131286", 
      "p131294"),
    
    c("p131666"),
    
    c("p40000_i0", 
      "p40000_i1"), 
    
    c("p132030",
      "p132032",
      "p132034"), 
    
    paste0("p40005_i", 0:21)
  )
  names(cox.outcomes) <- c("dementia", "type_2_diabetes", "stroke", 
                           "myocardial_infarction", "COPD", "heart_failure", 
                           "hypertension", "cirrhosis_fibrosis", "mortality", 
                           "kidney_failure", "cancer")
  return(cox.outcomes)
}

annotate_outcomes <- function(bd, cox.outcomes, visit = "first"){
  
  if (visit == "first")
    instance <- 0
  else if (visit == "third")
    instance <- 2
  else if (visit == "fourth")
    instance <- 3
  else
    stop("Check if the specified visit is correct. Extend implementation if needed.")
  
  # https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=819:
  # These have the disease, but date is incorrect, so we remove them for Cox PH model
  special.values <- as.numeric(as.Date(c("1900-01-01", "1901-01-01", "1902-02-02", "1903-03-03", "1909-09-09", "2037-07-07"), format = "%Y-%m-%d"))
  latest.date <- rep(NA, length(cox.outcomes))
  names(latest.date) <- names(cox.outcomes)
  
  dead <- bd[, c("p40000_i0", "p40000_i1"), drop = FALSE]
  dimnames <- dimnames(dead)
  dead <- apply(dead, 2, function(x){return(as.numeric(as.Date(x, format = "%Y-%m-%d")))})
  dimnames(dead) <- dimnames # Useful to keep for debugging purposes
  # We can keep this here because there are no special values for mortality, but do check if you add new fields!!!
  special.indices <- rowMeans(apply(dead, 2, function(x){x %in% special.values})) > 0
  dead[dead %in% special.values] <- NA
  
  # Take the average of the age of death if both are close to each other
  # range(dead[, 1] - dead[, 2], na.rm = TRUE)
  # 0 0: so no problem, but just in case, let's remove it if there would be over a year's difference:
  dead[which((dead[, 1] - dead[, 2]) > 365), ] <- c(NA, NA)
  
  event.time.dead <- rowMeans(dead, na.rm = TRUE)
  event.time.dead[is.infinite(event.time.dead)] <- NA
  
  status.dead <- ifelse(is.na(event.time.dead), 0, 1)
  status.dead[which((dead[, 1] - dead[, 2]) > 365)] <- 1 # These would still be dead, but without event time
  
  for(i in 1:length(cox.outcomes)){
    
    tmp <- bd[, cox.outcomes[[i]], drop = FALSE]
    dimnames <- dimnames(tmp)
    tmp <- apply(tmp, 2, function(x){return(as.numeric(as.Date(x, format = "%Y-%m-%d")))})
    dimnames(tmp) <- dimnames # Useful to keep for debugging purposes
    # We can keep this here because there are no special values for mortality, but do check if you add new fields!!!
    special.indices <- rowMeans(apply(tmp, 2, function(x){x %in% special.values})) > 0
    tmp[tmp %in% special.values] <- NA
    
    latest.date[i] <- max(tmp, na.rm = TRUE)
    
    if(names(cox.outcomes)[i] %in% c("mortality")){
      
      event.time <- event.time.dead
      status <- status.dead
      
      # For mortality, exclude "Chapter XX External causes of morbidity and mortality" as *primary* cause of death
      # We will put these at 0, but with the drop-out point being their time of death
      status[which(grepl("^(V|W|X|Y)", bd$p40001_i0) & (grepl("^(V|W|X|Y)", bd$p40001_i1) | is.na(bd$p40001_i1)))] <- 0
      
    } else{
      
      # Take the first occurrence of the disease
      event.time <- rowMins(tmp, na.rm = TRUE)
      event.time[is.infinite(event.time)] <- NA
      
      status <- ifelse(is.na(event.time), 0, 1)
    }
    
    # NA for event.time with at least one "special date" should get status 1, and as a consequence have no date imputation for status 0
    # A non-NA value with at least one "special date" will have the other non-special date with status 1
    status[which(is.na(event.time) & special.indices)] <- 1
    
    # For people with the disease, "time" is the time of the event minus the time of first visit
    time <- event.time-as.numeric(as.Date(bd[,paste0("p53_i", instance)]))
    
    # For people still alive without the disease for whom the event time has not been overruled because of a previous exception, 
    # take the latest date the disease was reported minus the date of first visit
    time[((status == 0) & is.na(event.time) & (status.dead == 0))] <- (latest.date[i]-as.numeric(as.Date(bd[,paste0("p53_i", instance)])))[((status == 0) & is.na(event.time) & (status.dead == 0))]
    
    # For people who died without the disease for whom the event time has not been overruled because of a previous exception, 
    # take the latest date of death minus the date of first visit
    time[((status == 0) & is.na(event.time) & (status.dead == 1))] <- (event.time.dead-as.numeric(as.Date(bd[,paste0("p53_i", instance)])))[((status == 0) & is.na(event.time) & (status.dead == 1))]
    
    bd[, paste0("status_", names(cox.outcomes)[i])] <- status
    bd[, paste0("time_", names(cox.outcomes)[i])] <- time
    
  }
  return(bd)
}