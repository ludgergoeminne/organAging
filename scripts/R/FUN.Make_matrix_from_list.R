
##########################################################################
### This script creates matrix from list with data on individual genes ###
##########################################################################

make_matrix_from_list <- function(data_list,column_name){
  temp_rownames <- c()
  for (i in names(data_list)){
    temp_rownames <- union(temp_rownames,rownames(data_list[[i]]))
  }
  temp_colnames <- names(data_list)
  temp_matrix <- matrix(NA,nrow=length(temp_rownames),ncol=length(temp_colnames))
  rownames(temp_matrix) <- temp_rownames
  colnames(temp_matrix) <- temp_colnames
  for (j in colnames(temp_matrix)){
    temp_dataset <- data_list[[j]]
    temp_matrix[rownames(temp_dataset),j] <- temp_dataset[,column_name]
  }
  temp_matrix <- as.data.frame(temp_matrix)
  return(temp_matrix)
}


