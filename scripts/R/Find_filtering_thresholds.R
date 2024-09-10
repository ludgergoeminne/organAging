
### This script contains the justification for our filtering thresholds ###
### You can run this optionally after 1_Prepare_Data.R and 2_Prepare_Longitudinal_Data.R ###

print(paste0("Dimension for instance 0 (rows: participants, cols: proteins) before filtering."))
print(dim(olink_data_unfiltered_0))
print(paste0("Dimension for instance 0 (rows: participants, cols: proteins) after filtering."))
print(dim(olink_data_0))

# "Dimension for instance 0 (rows: participants, cols: proteins) before filtering."
# 53014  2923
# "Dimension for instance 0 (rows: participants, cols: proteins) after filtering."
# 44952  2916

instances <- c(0,2,3)
for (instance in instances){
print(paste0("Dimension for longitudinal instance ", instance, " (rows: participants, cols: proteins) before filtering."))
print(dim(olink_data_unfiltered_list[[paste0("instance_", instance)]]))

print(paste0("Dimension for longitudinal instance ", instance, " (rows: participants, cols: proteins) after filtering."))
print(dim(olink_data_list[[paste0("instance_", instance)]]))
}

# "Dimension for instance 0 (rows: participants, cols: proteins) before filtering."
# 53014  1463
# "Dimension for instance 0 (rows: participants, cols: proteins) after filtering."
# 985 1451
# "Dimension for instance 2 (rows: participants, cols: proteins) before filtering."
# 1172 1463
# "Dimension for instance 2 (rows: participants, cols: proteins) after filtering."
# 985 1451
# "Dimension for instance 3 (rows: participants, cols: proteins) before filtering."
# 1123 1463
# "Dimension for instance 3 (rows: participants, cols: proteins
# 985 1451

values.0 <- as.matrix(olink_data_unfiltered_0)

hist(rowMeans(is.na(values.0)), breaks = 10000, main = "First visit", xlab = "Fraction of missing values")
abline(v = 0.499486828600753, col = "red")

table(rowMeans(is.na(values.0)))[920:940]
# Starting from "0.499486828600753", there is a marked increase in missing values, i.e. around 50%, likely due to enrichment in the smaller panel that was measured longitudinally

values.0.filtered <- values.0[rowMeans(is.na(values.0)) < 0.499486828600753,]

# Make a histogram of the fraction of missing values in the training dataset
hist(colMeans(is.na(values.0.filtered)), breaks = 100, main = "Fraction of missing values", xlab = "")
abline(v = 0.1, col = "red")

coding143[which(colMeans(is.na(values.0.filtered)) > 0.1), ]
# coding protein                                 protein.name
# 115     115   AMY2B                             Alpha-amylase 2B
# 709     709    CST1                                  Cystatin-SN
# 732     732    CTSS                                  Cathepsin S
# 1173   1173  GLIPR1        Glioma pathogenesis-related protein 1
# 1889   1889    NPM1                                Nucleophosmin
# 1991   1991  PCOLCE       Procollagen C-endopeptidase enhancer 1
# 2614   2614 TACSTD2 Tumor-associated calcium signal transducer 2


