
### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

empty.plot <- ggplot_gtable(ggplot_build(ggplot() + theme_void()))

### Arrange, plot, and save Supp. Fig. 1 ###
# Plot_Heatmaps_r_MAE.R
# Plot_Barplot_Other_Models.R
# Plot_External_Hazard_Ratios_and_Barplot.R
# Plot_Heatmaps_Organ_Ages.R

Fig_S1 <- ggarrange(ggarrange(heatmaps.r[["1st-generation models"]],
                              heatmaps.MAE[["1st-generation models"]],
                            ncol = 3, nrow = 1, labels = c(LETTERS[1:2], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                  ggarrange(heatmaps.r[["Mortality-based models"]],
                            heatmaps.MAE[["Mortality-based models"]],
                            empty.plot,
                            ncol = 3, nrow = 1, labels = c(LETTERS[3:4], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                  ggarrange(barplot.other.models,
                            barplots.comparison.Oh,
                            empty.plot,
                            ncol = 3, nrow = 1, labels = c(LETTERS[5:6], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                  ggarrange(heatmap.corr.organs$gen1$predicted,
                            heatmap.corr.organs$gen2$predicted,
                            ncol = 2, nrow = 1, labels = c(LETTERS[7:8]), widths = c(0.5, 0.5), font.label = list(size = 17)),
                  nrow = 4,  # Total number of rows!
                  heights = c(1, 1, 1, 1.5))

### Plot Figure S1 ###
plot(Fig_S1)

### Save Figure S1 ###
ggsave(Fig_S1, file = paste0(plot.dir, "Fig_S1.svg"), width = 13.5, height = 13, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig_S1.svg"))

### Arrange, plot, and save Supp. Fig. 2 ###
# Plot_Heatmaps_r_MAE_longitudinal.R

Fig_S2 <- ggarrange(ggarrange(heatmaps.r.longitudinal[["1st-generation models"]],
                              heatmaps.MAE.longitudinal[["1st-generation models"]],
                              ncol = 3, nrow = 1, labels = c(LETTERS[1:2], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                    ggarrange(heatmaps.r.longitudinal[["Mortality-based models"]],
                              heatmaps.MAE.longitudinal[["Mortality-based models"]],
                              empty.plot,
                              ncol = 3, nrow = 1, labels = c(LETTERS[3:4], ""), widths = c(0.45, 0.45, 0.1), font.label = list(size = 17)),
                    nrow = 2,  # Total number of rows!
                    heights = c(1, 1))

### Plot Figure S2 ###
plot(Fig_S2)

### Save Figure S2 ###
ggsave(Fig_S2, file = paste0(plot.dir, "Fig_S2.svg"), width = 13.5, height = 8.1, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig_S2.svg"))

### Arrange, plot, and save Supp. Fig. 3 ###
# 5_Analyze_Longitudinal_Data.R

Fig_S3 <- ggarrange(
  ggarrange(empty.plot, prediction.plots.longitudinal.gen2[["training dataset"]], empty.plot,
            empty.plot, prediction.plots.longitudinal.gen2[["test dataset"]], empty.plot,
            labels = c(LETTERS[1], "", "", LETTERS[2], "", ""), font.label = list(size = 17), ncol = 6, widths = c(2.5/100, 0.435, 0.04, 2.5/100, 0.435, 0.04)),
  ggarrange(empty.plot, prediction.plots.longitudinal.gen2[["third visit"]], empty.plot,
            empty.plot, prediction.plots.longitudinal.gen2[["fourth visit"]], empty.plot,
            labels = c(LETTERS[3], "", "", LETTERS[4], "", ""), font.label = list(size = 17), ncol = 6, widths = c(2.5/100, 0.435, 0.04, 2.5/100, 0.435, 0.04)),
  nrow = 2,  # Total number of rows!
  heights = c(1, 1))

### Plot Figure S3 ###
plot(Fig_S3)

### Save Figure S3 ###
ggsave(Fig_S3, file = paste0(plot.dir, "Fig_S3.svg"), width = 13, height = 6.4, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig_S3.svg"))

### Arrange, plot, and save Supp. Fig. 4 ###
# 5_Analyze_Longitudinal_Data.R
# Plot_Violins_COVID.R

Fig_S4 <- ggarrange(
  ggarrange(plotlist = hist.slopes.gen2.longitudinal.gen2[names(organ.proteins.selected)], 
            labels = LETTERS[(1:length(names(organ.proteins.selected)))], font.label = list(size = 17)),
  ggarrange(empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[1]]], empty.plot, 
            empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[2]]], empty.plot, 
            empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[3]]], empty.plot, 
            labels = c(LETTERS[10], "", "", LETTERS[11], "", "", LETTERS[12], "", ""), widths = c(1.2/60, 29/100, 2/75, 1.2/60, 29/100, 2/75, 1.2/60, 29/100, 2/75),
            ncol = 9, nrow = 1, font.label = list(size = 17)), 
  ggarrange(empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[4]]], empty.plot, 
            empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[5]]], empty.plot, 
            empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[6]]], empty.plot, 
            labels = c(LETTERS[13], "", "", LETTERS[14], "", "", LETTERS[15], "", ""), widths = c(1.2/60, 29/100, 2/75, 1.5/60, 29/100, 2/75, 1.5/60, 29/100, 2/75),
            ncol = 9, nrow = 1, font.label = list(size = 17)), 
  ggarrange(empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[7]]], empty.plot, 
            empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[8]]], empty.plot, 
            empty.plot, plots.longitudinal.change.vs.first.visit[[names(organ.proteins.selected)[9]]], empty.plot, 
            labels = c(LETTERS[16], "", "", LETTERS[17], "", "", LETTERS[18], "", ""), widths = c(1.2/60, 29/100, 2/75, 1.2/60, 29/100, 2/75, 1.2/60, 29/100, 2/75),
            ncol = 9, nrow = 1, font.label = list(size = 17)), 
      nrow = 4,  # Total number of rows!
      heights = c(3*0.9, 1, 1, 1)
      )

### Plot Figure S4 ###
plot(Fig_S4)

### Save Figure S4 ###
ggsave(Fig_S4, file = paste0(plot.dir, "Fig_S4.svg"), width = 13.5, height = 18.5, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig_S4.svg"))

### Arrange, plot, and save Supp. Fig. 5 ###
# Plot_Violins_COVID.R

Fig_S5 <- ggarrange(
  ggarrange(plotlist = list(empty.plot, violin.plots.COVID$Age.cat[["Conventional"]],
                            empty.plot, violin.plots.COVID$Age.cat[["Lung"]]), 
            labels = c(LETTERS[1], "", LETTERS[2], ""), widths = c(2/90, 43/90, 2/90, 43/90),
            ncol = 4, nrow = 1, font.label = list(size = 17)),
  ggarrange(empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[1]]], 
                            empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[2]]], 
                            labels = c(LETTERS[3], "", LETTERS[4], ""), widths = c(2/90, 43/90, 2/90, 43/90), 
                            font.label = list(size = 17), nrow = 1, ncol = 4),
  ggarrange(empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[3]]], 
            empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[4]]], 
            labels = c(LETTERS[5], "", LETTERS[6], ""), widths = c(2/90, 43/90, 2/90, 43/90), 
            font.label = list(size = 17), nrow = 1, ncol = 4),
  ggarrange(empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[5]]], 
            empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[6]]], 
            labels = c(LETTERS[7], "", LETTERS[8], ""), widths = c(2/90, 43/90, 2/90, 43/90), 
            font.label = list(size = 17), nrow = 1, ncol = 4),
  ggarrange(empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[7]]], 
            empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[8]]], 
            labels = c(LETTERS[9], "", LETTERS[10], ""), widths = c(2/90, 43/90, 2/90, 43/90), 
            font.label = list(size = 17), nrow = 1, ncol = 4),
  ggarrange(empty.plot, violin.plots.COVID$Acuity.max[[names(organ.proteins.selected)[9]]], 
            empty.plot, empty.plot, 
            labels = c(LETTERS[11], "", "", ""), widths = c(2/90, 43/90, 2/90, 43/90), 
            font.label = list(size = 17), nrow = 1, ncol = 4),
            nrow = 6, 
            ncol = 1)

plot(Fig_S5)

### Save Figure S5 ###
ggsave(Fig_S5, file = paste0(plot.dir, "Fig_S5.svg"), width = 4.12*2*1.2, height = 15.6, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig_S5.svg"))

### Arrange, plot, and save Supp. Fig. 6 ###
# 7_Analyze_CSF_Models.R

Fig_S6 <- ggarrange(ggarrange(empty.plot,
                              prediction.plots.CSF[["first visit"]]$Conventional,
                              empty.plot,
                              prediction.plots.CSF[["last visit"]]$Conventional,
                              labels = c(LETTERS[1], "", LETTERS[2], ""), 
                              font.label = list(size = 17),
                              nrow = 1, 
                              ncol = 4, widths = c(2/90, 43/90, 2/90, 43/90)), 
                    ggarrange(empty.plot,
                              violin.model.CSF[["Conventional"]],
                              empty.plot,
                              empty.plot,
                              labels = c(LETTERS[3], "", "", ""), 
                              font.label = list(size = 17),
                              nrow = 1, 
                              ncol = 4, widths = c(2/90, 43/90, 2/90, 43/90)), 
                    ggarrange(empty.plot,
                              prediction.plots.CSF[["first visit"]]$Brain,
                              empty.plot,
                              prediction.plots.CSF[["last visit"]]$Brain,
                              labels = c(LETTERS[4], "", LETTERS[5], ""),
                              font.label = list(size = 17),
                              nrow = 1, 
                              ncol = 4, widths = c(2/90, 43/90, 2/90, 43/90)),
                    ggarrange(empty.plot,
                              violin.model.CSF[["Brain"]],
                              empty.plot,
                              empty.plot,
                              labels = c(LETTERS[6], "", "", ""),
                              font.label = list(size = 17),
                              nrow = 1, 
                              ncol = 4, widths = c(2/90, 43/90, 2/90, 43/90)),
                    nrow = 4,  # Total number of rows!
                    heights = c(1, 1, 1, 1)
)

plot(Fig_S6)

### Save Figure S6 ###
ggsave(Fig_S6, file = paste0(plot.dir, "Fig_S6.svg"), width = 10, height = 11.52, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "Fig_S6.svg"))

### Arrange, plot, and save Supp. Fig. 7 ###
# Plot_Boxplots_Sarcoidosis.R

FigS7 <- ggpubr::ggarrange(
  ggarrange(empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[1]]],
            empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[2]]],
            empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[3]]],
            ncol = 6, nrow = 1, labels = c(LETTERS[1], "", LETTERS[2], "", LETTERS[3], ""), widths = c(0.05*1/3, 0.95*1/3, 0.05*1/3, 0.95*1/3, 0.05*1/3, 0.95*1/3), font.label = list(size = 17)),
  ggarrange(empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[4]]],
            empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[5]]],
            empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[6]]],
            ncol = 6, nrow = 1, labels = c(LETTERS[4], "", LETTERS[5], "", LETTERS[6], ""), widths = c(0.05*1/3, 0.95*1/3, 0.05*1/3, 0.95*1/3, 0.05*1/3, 0.95*1/3), font.label = list(size = 17)),
  ggarrange(empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[7]]],
            empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[8]]],
            empty.plot,
            violin.plots.GSE169148[[names(organ.proteins.selected)[9]]],
            ncol = 6, nrow = 1, labels = c(LETTERS[7], "", LETTERS[8], "", LETTERS[9], ""), widths = c(0.05*1/3, 0.95*1/3, 0.05*1/3, 0.95*1/3, 0.05*1/3, 0.95*1/3), font.label = list(size = 17)),
  nrow = 3,  # Total number of rows!
  heights = c(1, 1, 1)
)

plot(FigS7)

### Save Figure S7 ###
ggsave(FigS7, file = paste0(plot.dir, "FigS7.svg"), width = 4.12*3*1.2, height = 2.88*3, limitsize = FALSE)
preprocess_SVG_file(paste0(plot.dir, "FigS7.svg"))

