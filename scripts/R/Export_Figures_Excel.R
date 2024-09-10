### This script contains the code to export the figure data to Excel ###
# See Arrange_Figures.R and Arrange_Supp_Figures.R for more information on how the figures were composed

##################################
### 1. Export the main figures ###
##################################

### 1.1. Export Figure 1 ###

tmp.list <- list(
  attr(prediction.plots.gen1[["training dataset"]]$Conventional, "data"), 
  attr(prediction.plots.gen1[["test dataset"]]$Conventional, "data"), 
  attr(violins.sex.gen1$Conventional, "data"), 
  attr(heatmap.traits, "data"), 
  attr(HR.plot.gen1, "data"), 
  attr(prediction.plots.gen2[["training dataset"]]$Conventional, "data"), 
  attr(prediction.plots.gen2[["test dataset"]]$Conventional, "data"), 
  attr(violins.sex.gen2$Conventional, "data"), 
  attr(HR.plot.gen2, "data"), 
  attr(density.plot.gen2, "data")
)
names(tmp.list) <- paste0("Fig. 1", LETTERS[1:10])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Figure 1. Note that UK Biobank does not allow us to provide participant-level data. Therefore, we provide nearly-equal k-means summarized data (k = 10) for Fig. 1A, B, F, G, and quantiles for Fig. 1C, H.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_1.xlsx"))

### 1.2. Export Figure 2 ###

tmp.list <- list(
  attr(volcano.age, "data"),
  attr(barplot.fgsea.age, "data"),
  attr(volcano.mortality, "data"),
  attr(barplot.fgsea.mortality, "data"),
  attr(heatmap.genes, "data"),
  attr(heatmap.functions, "data"),
  attr(violin.secreted.aging, "data"),
  attr(violin.secreted.mortality, "data"),
  attr(corrplot.nonsecreted, "data"),
  attr(corrplot.secreted, "data")
)
names(tmp.list) <- paste0("Fig. 2", LETTERS[1:10])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Figure 2. 'df': degrees of freedom, 'll': lower limit of the 95% confidence interval, 'ul': upper limit of the 95% confidence interval, 'r': correlation coefficient, 'qval' and 'padj': Benjamini-Hochberg FDR-adjusted p-value. 'coef' in Fig. 2C is the estimated log(mortality hazard ratio).") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_2.xlsx"))

### 1.3. Export Figure 3 ###

tmp.list <- list(
  attr(barplots.age.r[["gen1"]], "data"),
  attr(barplots.age.r[["gen2"]], "data"),
  attr(prediction.plots.gen1[["test dataset"]]$Brain, "data"),
  attr(prediction.plots.gen1[["test dataset"]]$Artery, "data"),
  attr(prediction.plots.gen2[["test dataset"]]$Brain, "data"),
  attr(prediction.plots.gen2[["test dataset"]]$Artery, "data"),
  attr(heatmap.corr.organs$gen1$residual, "data"),
  attr(heatmap.corr.organs$gen2$residual, "data"),
  attr(heatmap.sex.deviations.gen1, "data"),
  attr(heatmap.sex.deviations.gen2, "data"),
  attr(HR.plot.mortality.other.papers, "data")
)
names(tmp.list) <- paste0("Fig. 3", LETTERS[1:11])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Figure 3. Note that UK Biobank does not allow us to provide participant-level data. Therefore, we provide nearly-equal k-means summarized data (k = 10) for Fig. 3C, D, E, F. 'df': degrees of freedom, 'll': lower limit of the 95% confidence interval, 'ul': upper limit of the 95% confidence interval, 'r': correlation coefficient, 'qval' and 'padj': Benjamini-Hochberg FDR-adjusted p-value. 'coef' in Fig. 3K is the estimated log(mortality hazard ratio).") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_3.xlsx"))

### 1.4. Export Figure 4 ###

tmp.list <- list(
  attr(HR.plots.organs.combined$Liver, "data"), 
  attr(HR.plots.organs.combined$Kidney, "data"), 
  attr(HR.plots.organs.combined$Lung, "data"), 
  attr(heatmaps.HR.disease$gen1, "data"), 
  attr(heatmaps.HR.disease$gen2, "data"), 
  attr(scatterplot.leaveout[[1]], "data"), 
  attr(scatterplot.leaveout[[2]], "data"), 
  attr(scatterplot.leaveout[[3]], "data"), 
  attr(scatterplot.leaveout[[4]], "data")
)
names(tmp.list) <- paste0("Fig. 4", LETTERS[1:9])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Figure 4. 'df': degrees of freedom, 'll': lower limit of the 95% confidence interval, 'ul': upper limit of the 95% confidence interval, 'r': correlation coefficient, 'qval': Benjamini-Hochberg FDR-adjusted p-value. 'coef' in Fig. 4A - C is the estimated log(hazard ratio) per unit of standard deviation; 'coef' in Fig. 4D - E is the estimated log(hazard ratio) per unit of predicted biological age (D) and predicted log(mortality hazard) (E).") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_4.xlsx"))

### 1.5. Export Figure 5 ###

tmp.list <- list(
  attr(heatmap.top.jobs.gen2, "data"), 
  attr(heatmap.top.foods.gen2, "data"), 
  attr(heatmap.top.medications.gen2, "data"), 
  attr(heatmap.metformin.simvastatin, "data"), 
  attr(barplot.smoking.current.never.gen2, "data"), 
  attr(barplot.drinking.gen2, "data"), 
  attr(scatterplot.yoghurt.drinking.leaveout[[1]], "data"), 
  attr(scatterplot.yoghurt.drinking.leaveout[[3]], "data"), 
  attr(scatterplot.yoghurt.drinking.leaveout[[2]], "data"), 
  attr(scatterplot.yoghurt.drinking.leaveout[[4]], "data")
)
names(tmp.list) <- paste0("Fig. 5", LETTERS[1:10])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Figure 5. 'df': degrees of freedom, 'll': lower limit of the 95% confidence interval, 'ul': upper limit of the 95% confidence interval, 'r': correlation coefficient, 'qval': Benjamini-Hochberg FDR-adjusted p-value.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_5.xlsx"))

### 1.6. Export Figure 6 ###

tmp.list <- list(
  attr(hist.slopes.gen2.longitudinal.gen2[["Conventional"]], "data"), 
  attr(heatmap.sex.deviations.gen2.longitudinal, "data"), 
  attr(violin.plots.slopes.gen2.category[["Conventional"]], "data"), 
  attr(correlation.slopes.gen2.n.diseases, "data")
)
names(tmp.list) <- paste0("Fig. 6", LETTERS[1:4])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Figure 6. 'qval': Benjamini-Hochberg FDR-adjusted p-value. Note that UK Biobank does not allow us to provide participant-level data. Therefore, we provide quantiles for Fig. 6C and nearly-equal k-means summarized data (k = 10) for Fig. 6D.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_6.xlsx"))

### 1.7. Export Figure 7 ###

tmp.list <- list(
  attr(violin.plots.COVID$COVID$Conventional, "data"), 
  attr(violin.plots.COVID$COVID$Lung, "data"), 
  attr(violin.plots.COVID$Acuity.max$Conventional, "data"), 
  attr(violin.plots.COVID$Acuity.max$Lung, "data"), 
  attr(scatterplot.COVID.leaveout[[1]], "data"), 
  attr(scatterplot.COVID.leaveout[[2]], "data"), 
  attr(scatterplot.COVID.leaveout[[3]], "data"), 
  attr(scatterplot.COVID.leaveout[[4]], "data"), 
  attr(prediction.plots.Parkinson[["plasma"]]$Conventional, "data"), 
  attr(prediction.plots.Parkinson[["CSF"]]$Conventional, "data"), 
  attr(violin.Parkinson.model[["plasma"]]$Conventional, "data"), 
  attr(violin.Parkinson.model[["CSF"]]$Conventional, "data"), 
  attr(prediction.plots.Parkinson[["plasma"]]$Brain, "data"), 
  attr(prediction.plots.Parkinson[["CSF"]]$Brain, "data"), 
  attr(violin.Parkinson.model[["plasma"]]$Brain, "data"), 
  attr(violin.Parkinson.model[["CSF"]]$Brain, "data"), 
  attr(scatterplot.PD.MSA.leaveout[[1]], "data"), 
  attr(scatterplot.PD.MSA.leaveout[[2]], "data"), 
  attr(scatterplot.PD.MSA.leaveout[[3]], "data"), 
  attr(scatterplot.PD.MSA.leaveout[[4]], "data")
)
names(tmp.list) <- paste0("Fig. 7", LETTERS[c(1:8, seq(9, 15, by = 2), seq(10, 16, by = 2), 17:20)])
tmp.list <- tmp.list[sort(names(tmp.list))]
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Figure 7. 'qval': Benjamini-Hochberg FDR-adjusted p-value. Note that the individual-level data in the Dammer et al. (2022) dataset can only be accessed by registered Synapse users. Therefore, we provide quantiles for Fig. 7 M - P and nearly-equal k-means summarized data (k = 5) for Fig. 7I - L.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_7.xlsx"))

###########################################
### 2. Export the supplementary figures ###
###########################################

### 2.1. Export Supplementary Figure S1 ###

tmp.list <- list(
  attr(heatmaps.r[["1st-generation models"]], "data"), 
  attr(heatmaps.r[["1st-generation models"]], "data"), # !Exceptional!: should be heatmaps.MAE, but heatmaps.r contains the same info and more (q-value and stars) 
  attr(heatmaps.r[["Mortality-based models"]], "data"), 
  attr(heatmaps.r[["Mortality-based models"]], "data"), # !Exceptional!: should be heatmaps.MAE, but heatmaps.r contains the same info and more (q-value and stars)
  attr(barplot.other.models, "data"), 
  attr(barplots.comparison.Oh, "data"), 
  attr(heatmap.corr.organs$gen1$predicted, "data"), 
  attr(heatmap.corr.organs$gen2$predicted, "data")
)
names(tmp.list) <- paste0("Fig. S1", LETTERS[1:8])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Supplementary Figure S1. 'r': correlation coefficient, 'r²': coefficient of determination, 'MAE': median absolute error of the residuals, 'MSE': mean squared error of the residuals. 'qval': Benjamini-Hochberg FDR-adjusted p-value.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_S1.xlsx"))

### 2.2. Export Supplementary Figure S2 ###

tmp.list <- list(
  attr(heatmaps.r.longitudinal[["1st-generation models"]], "data"), 
  attr(heatmaps.r.longitudinal[["1st-generation models"]], "data"), # !Exceptional!: should be heatmaps.MAE.longitudinal, but heatmaps.r.longitudinal contains the same info and more (q-value and stars) 
  attr(heatmaps.r.longitudinal[["Mortality-based models"]], "data"), 
  attr(heatmaps.r.longitudinal[["Mortality-based models"]], "data") # !Exceptional!: should be heatmaps.MAE.longitudinal, but heatmaps.r.longitudinal contains the same info and more (q-value and stars) 
)
names(tmp.list) <- paste0("Fig. S2", LETTERS[1:4])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Supplementary Figure S2. 'r': correlation coefficient, 'r²': coefficient of determination, 'MAE': median absolute error of the residuals, 'MSE': mean squared error of the residuals. 'qval': Benjamini-Hochberg FDR-adjusted p-value.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_S2.xlsx"))

### 2.3. Export Supplementary Figure S3 ###

tmp.list <- list(
  attr(prediction.plots.longitudinal.gen2[["training dataset"]], "data"), 
  attr(prediction.plots.longitudinal.gen2[["test dataset"]], "data"), 
  attr(prediction.plots.longitudinal.gen2[["third visit"]], "data"), 
  attr(prediction.plots.longitudinal.gen2[["fourth visit"]], "data")
)
names(tmp.list) <- paste0("Fig. S3", LETTERS[1:4])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Supplementary Figure S3.  Note that UK Biobank does not allow us to provide participant-level data. Therefore, we provide nearly-equal k-means summarized data (k = 10) for Supplementary Figure S3A - S3D.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_S3.xlsx"))

### 2.4. Export Supplementary Figure S4 ###

tmp.list <- list(
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[1]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[2]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[3]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[4]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[5]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[6]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[7]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[8]]], "data"), 
  attr(hist.slopes.gen2.longitudinal.gen2[[names(organ.proteins.selected)[9]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[1]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[2]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[3]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[4]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[5]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[6]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[7]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[8]]], "data"), 
  attr(plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[9]]], "data")
)
names(tmp.list) <- paste0("Fig. S4", LETTERS[1:18])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Supplementary Figure S4.  Note that UK Biobank does not allow us to provide participant-level data. Therefore, we provide nearly-equal k-means summarized data (k = 10) for Supplementary Fig. 4J - R.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_S4.xlsx"))

### 2.5. Export Supplementary Figure S5 ###

tmp.list <- list(
  attr(violin.plots.COVID$Age.cat[["Conventional"]], "data"), 
  attr(violin.plots.COVID$Age.cat[["Lung"]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[1]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[2]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[3]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[4]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[5]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[6]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[7]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[8]]], "data"), 
  attr(violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[9]]], "data")
)
names(tmp.list) <- paste0("Fig. S5", LETTERS[1:11])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Supplementary Figure S5.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_S5.xlsx"))

### 2.6. Export Supplementary Figure S6 ###

tmp.list <- list(
  attr(prediction.plots.CSF[["first visit"]]$Conventional, "data"), 
  attr(prediction.plots.CSF[["last visit"]]$Conventional, "data"), 
  attr(violin.model.CSF[["Conventional"]], "data"), 
  attr(prediction.plots.CSF[["first visit"]]$Brain, "data"), 
  attr(prediction.plots.CSF[["last visit"]]$Brain, "data"), 
  attr(violin.model.CSF[["Brain"]], "data")
)
names(tmp.list) <- paste0("Fig. S6", LETTERS[1:6])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Supplementary Figure S6. Note that the individual-level data in the Dammer et al. (2022) dataset can only be accessed by registered Synapse users. Therefore, we provide quantiles for Supplementary Fig. S6C and S6F and nearly-equal k-means summarized data (k = 5) for Supplementary Fig. S6A, S6B, S6D, and S6E.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_S6.xlsx"))

### 2.7. Export Supplementary Figure S7 ###

tmp.list <- list(
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[1]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[2]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[3]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[4]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[5]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[6]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[7]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[8]]], "data"), 
  attr(violin.plots.GSE169148[[names(organ.proteins.selected)[9]]], "data")
)
names(tmp.list) <- paste0("Fig. S7", LETTERS[1:9])
tmp.list <- lapply(tmp.list, as.data.frame)
# To quickly check if the list is complete
lapply(tmp.list, head)

wb <- createWorkbook()
openxlsx::addWorksheet(wb, "Caption")
openxlsx::writeData(wb, "Caption", "This Excel file contains the underlying data of Supplementary Figure S7.") # , headerStyle = header_st
header_st <- createStyle(textDecoration = "Bold")

for(i in 1:length(tmp.list)){
  openxlsx::addWorksheet(wb, names(tmp.list)[i])
  openxlsx::writeData(wb, names(tmp.list)[i], tmp.list[[i]], headerStyle = header_st)
}

rm(tmp.list)
openxlsx::saveWorkbook(wb, file = paste0(table.dir, "Fig_S7.xlsx"))

