
empty.plot <- ggplot_gtable(ggplot_build(ggplot() + theme_void()))

### Arrange, plot, and save Fig. 1 ###
# Plot_Age_Correlations.R
# Plot_Violins_Sex.R
# Plot_Heatmap_Traits.R
# Plot_Hazard_Ratios.R
# Plot_Densities_Exreme_Agers.R

Fig1 <- ggarrange(
  ggarrange(prediction.plots.gen1[["training dataset"]]$Conventional,
                            prediction.plots.gen1[["test dataset"]]$Conventional,
                            violins.sex.gen1$Conventional,
                            ncol = 3, nrow = 1, labels = c(LETTERS[1:3]), widths = c(1/3, 1/3, 1/3), font.label = list(size = 17)),
  ggarrange(heatmap.traits,
            HR.plot.gen1,
            ncol = 2, nrow = 1, labels = c(LETTERS[4:5]), widths = c(1/2, 1/2), font.label = list(size = 17)),
  ggarrange(prediction.plots.gen2[["training dataset"]]$Conventional,
            prediction.plots.gen2[["test dataset"]]$Conventional,
            violins.sex.gen2$Conventional,
            ncol = 3, nrow = 1, labels = c(LETTERS[6:8]), widths = c(1/3, 1/3, 1/3), font.label = list(size = 17)),
  ggarrange(HR.plot.gen2,
            density.plot.gen2,
            ncol = 2, nrow = 1, labels = c(LETTERS[9:10]), widths = c(1/2, 1/2), font.label = list(size = 17)),
  nrow = 4,  # Total number of rows!
  heights = c(1, 1, 1, 1)
)

### Plot Figure 1 ###
plot(Fig1)

### Save Figure 1 ###
ggsave(Fig1, file = paste0(plot.dir, "Fig1.svg"), width = 11.20, height = 14.90, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig1.svg"))

### Arrange, plot, and save Fig. 2 ###
# Plot_Volcanos_GSEA.R
# Plot_Heatmap_Genes_Functions.R
# Plot_Violins_Scatter_Secretion.R

Fig2 <- ggarrange(
  ggarrange(empty.plot,
            volcano.age,
            barplot.fgsea.age,
            ncol = 3, nrow = 1, labels = c(LETTERS[1], "", LETTERS[2]), widths = c(2/90, 88*1/3/90, 88*2/3/90), font.label = list(size = 17)),
  ggarrange(empty.plot,
            volcano.mortality,
            barplot.fgsea.mortality,
            ncol = 3, nrow = 1, labels = c(LETTERS[3], "", LETTERS[4]), widths = c(2/90, 88*1/3/90, 88*2/3/90), font.label = list(size = 17)),
  ggarrange(heatmap.genes,
            heatmap.functions,
            ncol = 2, nrow = 1, labels = c(LETTERS[5:6]), widths = c(1/2, 1/2), font.label = list(size = 17)),
  ggarrange(empty.plot,
            violin.secreted.aging,
            empty.plot,
            empty.plot,
            violin.secreted.mortality,
            empty.plot,
            ncol = 6, nrow = 1, labels = c(LETTERS[7], "", "", LETTERS[8], "", ""), widths = c(2/90, 38/90, 5/90, 2/90, 38/90, 5/90), font.label = list(size = 17)),
  ggarrange(empty.plot,
            corrplot.nonsecreted,
            empty.plot,
            empty.plot,
            corrplot.secreted,
            empty.plot,
            ncol = 6, nrow = 1, labels = c(LETTERS[9], "", "", LETTERS[10], "", ""), widths = c(2/90, 38/90, 5/90, 2/90, 38/90, 5/90), font.label = list(size = 17)),
  nrow = 5,  # Total number of rows!
  heights = c(1, 1, 1.575, 1, 1)
)

### Plot Figure 2 ###
plot(Fig2)

### Save Figure 2 ###
ggsave(Fig2, file = paste0(plot.dir, "Fig2.svg"), width = 14.10, height = 18.82, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig2.svg"))

### Arrange, plot, and save Fig. 3 ###
# Plot_Barplots_Correlations.R
# Plot_Age_Correlations.R
# Plot_Heatmaps_Organ_Ages.R
# Plot_Heatmaps_Sex_Deviations.R
# Plot_External_Hazard_Ratios_and_Barplot.R

Fig3 <- ggarrange(ggarrange(barplots.age.r[["gen1"]],
                            barplots.age.r[["gen2"]],
                            ncol = 3, nrow = 1, labels = c(LETTERS[1:2], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                  ggarrange(prediction.plots.gen1[["test dataset"]]$Brain,
                            prediction.plots.gen1[["test dataset"]]$Artery,
                            empty.plot,
                            ncol = 3, nrow = 1, labels = c(LETTERS[3:4], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                  ggarrange(prediction.plots.gen2[["test dataset"]]$Brain,
                            prediction.plots.gen2[["test dataset"]]$Artery,
                            empty.plot,
                            ncol = 3, nrow = 1, labels = c(LETTERS[5:6], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                  ggarrange(heatmap.corr.organs$gen1$residual,
                            heatmap.corr.organs$gen2$residual,
                            ncol = 2, nrow = 1, labels = c(LETTERS[7:8]), widths = c(0.5, 0.5), font.label = list(size = 17)),
                  ggarrange(heatmap.sex.deviations.gen1,
                            heatmap.sex.deviations.gen2,
                            HR.plot.mortality.other.papers,
                            ncol = 3, nrow = 1, labels = c(LETTERS[9:11]), widths = c(0.25, 0.3, 0.45), font.label = list(size = 17)),
                  nrow = 5,  # Total number of rows!
                  heights = c(1, 1, 1, 1, 1))

### Plot Figure 3 ###
plot(Fig3)

### Save Figure 3 ###
ggsave(Fig3, file = paste0(plot.dir, "Fig3.svg"), width = 13.5, height = 18.75, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig3.svg"))

### Arrange, plot, and save Fig. 4 ###
# Plot_Hazard_Ratios.R
# Plot_Heatmaps_HR_Diseases.R
# Plot_Scatter_Leave_Proteins_Out.R

Fig4 <- ggarrange(ggarrange(HR.plots.organs.combined$Liver,
                            HR.plots.organs.combined$Kidney,
                            HR.plots.organs.combined$Lung,
                            ncol = 3, nrow = 1, labels = LETTERS[1:3], widths = c(0.3, 0.3, 0.3), font.label = list(size = 17)),
                  ggarrange(heatmaps.HR.disease$gen1,
                            heatmaps.HR.disease$gen2,
                            scatterplot.leaveout[[1]],
                            ncol = 3, nrow = 1, labels = LETTERS[4:6], widths = c(0.3, 0.3, 0.3), font.label = list(size = 17)), 
                  ggarrange(scatterplot.leaveout[[2]],
                            scatterplot.leaveout[[3]],
                            scatterplot.leaveout[[4]],
                            ncol = 3, nrow = 1, labels = LETTERS[7:9], widths = c(0.3, 0.3, 0.3), font.label = list(size = 17)),
                  nrow = 3,  # Total number of rows!
                  heights = c(1.2, 1, 1))
plot(Fig4)

### Save Figure 4 ###
ggsave(Fig4, file = paste0(plot.dir, "Fig4.svg"), width = 20, height = 11.60, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig4.svg"))

### Arrange, plot, and save Fig. 5 ###
# Plot_Heatmaps_Professions.R
# Plot_Heatmaps_Foods.R
# Plot_Heatmaps_Medications.R
# Plot_Heatmap_Metformin_Simvastatin.R
# Plot_Barplots_Smoking_Drinking.R
# Plot_Scatter_Yoghurt_Drinking_Leave_Protein_Out.R

Fig5 <- ggarrange(ggarrange(heatmap.top.jobs.gen2, ncol = 1, nrow = 1, labels = c(LETTERS[1]), widths = c(1), font.label = list(size = 17)),
                  ggarrange(heatmap.top.foods.gen2, empty.plot, ncol = 2, nrow = 1, labels = c(LETTERS[2], ""), widths = c(0.9, 0.1), font.label = list(size = 17)),
                  ggarrange(heatmap.top.medications.gen2, heatmap.metformin.simvastatin, ncol = 2, nrow = 1, labels = c(LETTERS[3:4], ""), widths = c(0.8, 0.2), font.label = list(size = 17)),
                  ggarrange(barplot.smoking.current.never.gen2, 
                            barplot.drinking.gen2, 
                            scatterplot.yoghurt.drinking.leaveout[[1]],
                            ncol = 3, nrow = 1, labels = c(LETTERS[5:7]), widths = c(0.3, 0.3, 0.3), font.label = list(size = 17)),
                  ggarrange(scatterplot.yoghurt.drinking.leaveout[[3]], 
                            scatterplot.yoghurt.drinking.leaveout[[2]], 
                            scatterplot.yoghurt.drinking.leaveout[[4]],
                            ncol = 3, nrow = 1, labels = c(LETTERS[8:10]), widths = c(0.3, 0.3, 0.3), font.label = list(size = 17)),
                  nrow = 5,  # Total number of rows!
                  heights = c(1, 1, 1, 0.7, 0.7))
# plot(Fig5)

### Save Figure 5 ###
ggsave(Fig5, file = paste0(plot.dir, "Fig5.svg"), width = 16, height = 25, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig5.svg"))

### Arrange, plot, and save Fig. 6 ###
# 5_Analyze_Longitudinal_Data.R

Fig6 <- ggarrange(ggarrange(hist.slopes.gen2.longitudinal.gen2[["Conventional"]], heatmap.sex.deviations.gen2.longitudinal, empty.plot, ncol = 3, nrow = 1, labels = c(LETTERS[1:2], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                  ggarrange(empty.plot, violin.plots.slopes.gen2.category[["Conventional"]], empty.plot, correlation.slopes.gen2.n.diseases, empty.plot, ncol = 5, nrow = 1, labels = c(LETTERS[3], "", LETTERS[4], "", ""), widths = c(2.5/100, 0.43, 2.5/100, 0.43, 0.09), font.label = list(size = 17)),
                  nrow = 2,  # Total number of rows!
                  heights = c(1, 1))
# plot(Fig6)

### Save Figure 6 ###
ggsave(Fig6, file = paste0(plot.dir, "Fig6.svg"), width = 9.5, height = 7, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig6.svg"))

### Arrange, plot, and save Fig. 7 ###

# Plot_Violins_COVID.R
# Plot_Violin_Scatter_PD_MSA.R
# Plot_Scatter_COVID_Leave_Protein_Out.R
# Plot_Scatter_PD_MSA_Leave_Protein_Out.R

Fig7 <- ggpubr::ggarrange(
  ggpubr::ggarrange(plotlist = list(violin.plots.COVID$COVID$Conventional, violin.plots.COVID$COVID$Lung, 
                                    violin.plots.COVID$Acuity.max$Conventional, violin.plots.COVID$Acuity.max$Lung), 
                    widths = c(0.5,0.5,1.2,1.2),
                    nrow = 1,
                    ncol = 4,
                    labels = LETTERS[1:4],
                    font.label = list(size = 17)), 
  ggpubr::ggarrange(plotlist = list(scatterplot.COVID.leaveout[[1]], scatterplot.COVID.leaveout[[2]], 
                                    scatterplot.COVID.leaveout[[3]], scatterplot.COVID.leaveout[[4]]), 
                    widths = c(1, 1, 1, 1),
                    nrow = 1,
                    ncol = 4,
                    labels = LETTERS[5:8],
                    font.label = list(size = 17)),
  ggpubr::ggarrange(prediction.plots.Parkinson[["plasma"]]$Conventional, 
                    prediction.plots.Parkinson[["CSF"]]$Conventional, 
                    violin.Parkinson.model[["plasma"]]$Conventional, 
                    violin.Parkinson.model[["CSF"]]$Conventional, 
                    prediction.plots.Parkinson[["plasma"]]$Brain, 
                    prediction.plots.Parkinson[["CSF"]]$Brain, 
                    violin.Parkinson.model[["plasma"]]$Brain, 
                    violin.Parkinson.model[["CSF"]]$Brain,
                    labels = LETTERS[c(seq(9, 15, by = 2), seq(10, 16, by = 2))],
                    nrow = 2, 
                    ncol = 4,
                    widths = c(1, 1, 0.65, 0.65), font.label = list(size = 17)), 
  ggpubr::ggarrange(plotlist = list(scatterplot.PD.MSA.leaveout[[1]], scatterplot.PD.MSA.leaveout[[2]], 
                                    scatterplot.PD.MSA.leaveout[[3]], scatterplot.PD.MSA.leaveout[[4]]), 
                    widths = c(1, 1, 1, 1),
                    nrow = 1,
                    ncol = 4,
                    labels = LETTERS[17:20],
                    font.label = list(size = 17)),
  nrow = 4,
  heights = c(1, 1, 2, 1)
)

plot(Fig7)

### Save Figure 7 ###
ggsave(Fig7, file = paste0(plot.dir, "Fig7.svg"), width = 4*4.12*0.9, height = 17, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig7.svg"))



