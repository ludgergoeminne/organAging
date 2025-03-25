
### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### This comes from Plot_Hazard_Ratios.R ###
HR.list <- readRDS(file = paste0(rds.dir, "HR_list.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

### 1. Barplot r our study vs Oh et al. ###

### 1.1. From Supplementary Table 8 and Supplementary Table 12 in Oh et al. (https://www.nature.com/articles/s41586-023-06802-1) ###
Oh_et_al_r <- openxlsx::read.xlsx(paste0(external.models.dir, ST8_Oh_et_al.file))
Oh_et_al_r <- Oh_et_al_r[, c("organ", "All_r")]
Oh_et_al_r$`Oh et al. (2023)` <- as.numeric(gsub("\\ \\(.+", "", Oh_et_al_r$All_r))
Oh_et_al_r$All_r <- NULL

### 1.2. Calculate out of fold correlations with age in our 1st-generation models ###

correlations.gen1.oof <- data.frame(organ = names(organ.proteins), 
                                    r = rep(NA, length(organ.proteins)))

for(k in 1:nrow(correlations.gen1.oof)){
  if(!(names(organ.proteins)[k] %in% c("Bladder"))){ # , "Thyroid"
    
    df <- data.frame(
      chronological_age = olink_bd_annotation_0$age_first_visit,
      biological_age = predicted.ages.oof$gen1$predicted[[k]],
      sex = olink_bd_annotation_0$p31
    )
    
    summary <- summary(lm(biological_age ~ chronological_age, data = df))
    # resid <- resid(lm(biological_age ~ chronological_age, data = df))
    # metrics <- paste0("r = ", sprintf("%.2f", round(sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared), 2)), ", r² = ", sprintf("%.2f", round(summary$r.squared, 2)), ", MAE = ", sprintf("%.2f", round(mean(abs(resid)), 2)))
    correlations.gen1.oof[k, "r"] <- sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared)
  }
}

correlation.comparison <- dplyr::inner_join(Oh_et_al_r, correlations.gen1.oof, by = "organ")
colnames(correlation.comparison)[3] <- "Out of fold"

plot.df <- tidyr::pivot_longer(correlation.comparison, cols = 2:ncol(correlation.comparison), names_to = "Study", values_to = "r")
plot.df <- plot.df[plot.df$organ %in% names(organ.proteins.selected),]
plot.df$organ <- factor(plot.df$organ, levels = rev(names(organ.proteins.selected)))

# plot.df$y.pos <- as.numeric(plot.df$organ)
# plot.df[plot.df$Study == "Oh et al. (2023)", "y.pos"] <- plot.df[plot.df$Study == "Oh et al. (2023)", "y.pos"]+0.125
# plot.df[plot.df$Study == "Out of fold", "y.pos"] <- plot.df[plot.df$Study == "Out of fold", "y.pos"]-0.125

barplots.comparison.Oh <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = r, y = organ, fill = Study)) + 
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#3C6682", "#45A778")) +
    xlab("Correlation with chronological age") +
    ylab("") +
    theme_bw() +
    theme(legend.position = c(0.85, 0.5), # The coordinates for legend.position are x- and y- offsets from the bottom-left of the plot, ranging from 0 - 1
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          # plot.title = element_text(face = "bold"), # element_blank(),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"), # , angle = 45, hjust = 1, vjust = 1
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold",size = 12),
          strip.background = element_blank())
))

attr(barplots.comparison.Oh, "data") <- plot.df

plot(barplots.comparison.Oh)


### 2. HR plot our study vs Oh et al. ###

# 1st generation & 2nd generation

# ST12 (effect sizes, etc.) and ST8 (standard deviation in the LonGenity cohort) from https://www.nature.com/articles/s41586-023-06802-1:
# "Age gaps were z-scored per aging model to account for the differences in model variability (Supplementary Fig. 3f)."
other.papers.mortality.sumstats <- read.csv(file = paste0(external.models.dir, mortality.sumstats.file), header = TRUE, quote = "")

# From: Supplementary Table 5 in https://www.aging-us.com/article/101684/text
other.papers.mortality.sumstats <- other.papers.mortality.sumstats[other.papers.mortality.sumstats$Study != "GrimAge",]

zval <- 3.93 # From: Supplementary Table 5 in https://www.aging-us.com/article/101684/text
se <- log(1.12)/zval # From: Supplementary Table 5 in https://www.aging-us.com/article/101684/text
abs(qnorm(0.05/2))*se

exp(log(1.12)-abs(qnorm(0.05/2))*se)
exp(log(1.12)+abs(qnorm(0.05/2))*se)

# Has to be GrimAge DNAm because of the 3.81 that is reported for this!
df.GrimAge <- data.frame(
  organ = c("Conventional"),
  Study = c("GrimAge"),
  beta = log(1.12),
  expbeta = 1.12,
  expci.lb = exp(log(1.12)-abs(qnorm(0.05/2))*se),
  expci.ub = exp(log(1.12)+abs(qnorm(0.05/2))*se),
  pheno_SD = 3.81  # From: Supplementary Table 9 in https://www.aging-us.com/article/101684/text
)
### The models from Oh et al. are already expressed per unit of s.d.!
# "A standard deviation increase (approximately four years of extra organ aging, Supplementary Table 8) in heart, adipose, liver, pancreas, brain, lung, immune or muscle age gap each conferred between 15–50% increased all-cause mortality risk."
# "Age gaps were z-scored per aging model to account for the differences in model variability (Supplementary Fig. 3f)."

### For GrimAge, what is reported is the hazard ratio associated with a 1 unit increase in the variable!
df.GrimAge$beta <- df.GrimAge$beta*df.GrimAge$pheno_SD
df.GrimAge$expbeta <- exp(log(df.GrimAge$expbeta)*df.GrimAge$pheno_SD)
df.GrimAge$expci.lb <- exp(log(df.GrimAge$expci.lb)*df.GrimAge$pheno_SD)
df.GrimAge$expci.ub <- exp(log(df.GrimAge$expci.ub)*df.GrimAge$pheno_SD)

other.papers.mortality.sumstats <- rbind(other.papers.mortality.sumstats, 
                                         df.GrimAge)

df.other.papers <- data.frame(
  outcome = "mortality",
  set = other.papers.mortality.sumstats$Study,
  `exp(coef)` = other.papers.mortality.sumstats$expbeta,
  coef = other.papers.mortality.sumstats$beta,
  z = NA,
  `Pr(>|z|)` = NA,
  ll = other.papers.mortality.sumstats$expci.lb,
  ul = other.papers.mortality.sumstats$expci.ub,
  qval = NA,
  organ = other.papers.mortality.sumstats$organ,
  model.type = other.papers.mortality.sumstats$Study,
  check.names = FALSE
)

df.gen1 <- HR.list[["gen1"]][names(organ.proteins.selected)]
for(i in 1:length(df.gen1)){
  df.gen1[[i]] <- df.gen1[[i]][df.gen1[[i]]$outcome == "mortality",]
  df.gen1[[i]]$organ <- names(df.gen1)[i]
}
df.gen1 <- do.call("rbind", df.gen1)
df.gen1$model.type <- "1st-generation model"

df.gen2 <- HR.list[["gen2"]][names(organ.proteins.selected)]
for(i in 1:length(df.gen2)){
  df.gen2[[i]] <- df.gen2[[i]][df.gen2[[i]]$outcome == "mortality",]
  df.gen2[[i]]$organ <- names(df.gen2)[i]
}
df.gen2 <- do.call("rbind", df.gen2)
df.gen2$model.type <- "Mortality-based model"

plot.df <- rbind(df.gen1, df.gen2, df.other.papers)
plot.df <- plot.df[plot.df$set != "Out of fold", ] # Remove the Out of fold for the plot
plot.df <- plot.df[plot.df$organ != "Skin", ] # Remove skin, as it does not exist in Oh et al.
plot.df$organ <- factor(plot.df$organ, levels = rev(names(organ.proteins.selected)[-which(names(organ.proteins.selected) == "Skin")]))

plot.df$y.pos <- as.numeric(plot.df$organ)
plot.df[plot.df$model.type == "Mortality-based model", "y.pos"] <- plot.df[plot.df$model.type == "Mortality-based model", "y.pos"]+0.125
plot.df[plot.df$model.type == "1st-generation model", "y.pos"] <- plot.df[plot.df$model.type == "1st-generation model", "y.pos"]-0.125

plot.df[plot.df$set == "training", "y.pos"] <- plot.df[plot.df$set == "training", "y.pos"]+0.125/2
plot.df[plot.df$set == "test", "y.pos"] <- plot.df[plot.df$set == "test", "y.pos"]-0.125/2

# Since GrimAge and Oh et al. conventional have more or less the same HR/unit sd, we also shift these a bit higher and lower
plot.df[plot.df$model.type == "GrimAge", "y.pos"] <- plot.df[plot.df$model.type == "GrimAge", "y.pos"]+0.125/2
plot.df[(plot.df$model.type == "Oh et al. (2023)") & (plot.df$organ == "Conventional"), "y.pos"] <- plot.df[(plot.df$model.type == "Oh et al. (2023)") & (plot.df$organ == "Conventional"), "y.pos"]-0.125/2

# plot.df$model.type[plot.df$model.type == "1st-generation model"] <- "1st-generation model\n(UK Biobank)"
# plot.df$model.type[plot.df$model.type == "Mortality-based model"] <- "Mortality-based model\n(UK Biobank)"
# plot.df$model.type[plot.df$model.type == "Oh et al. (2023)"] <- "Oh et al. (2023)\n(LonGenity cohort)"
# plot.df$model.type[plot.df$model.type == "GrimAge"] <- "GrimAge\n(Framingham Heart Study)"
# 
# plot.df$model.type <- factor(plot.df$model.type, levels = c(
#   "1st-generation model\n(UK Biobank)",
#   "Mortality-based model\n(UK Biobank)",
#   "Oh et al. (2023)\n(LonGenity cohort)",
#   "GrimAge\n(Framingham Heart Study)"
# ))

plot.df$model.type <- factor(plot.df$model.type, levels = c(
  "1st-generation model",
  "Mortality-based model",
  "Oh et al. (2023)",
  "GrimAge"
))

plot.df$shape <- "Test dataset" # Both Oh et al. and GrimAge were validated in a test dataset: the separate LonGenity cohort for Oh et al., and a separate test set in the Framingham Heart Study for GrimAge
plot.df$shape[plot.df$set == "Training set"] <- "Training dataset"
plot.df$shape[plot.df$set == "Test set"] <- "Test dataset"

plot.df$shape <- factor(plot.df$shape, levels = c("Training dataset", "Test dataset"))

lims <- c(0.9, 2.5)

HR.plot.mortality.other.papers <- ggplot_gtable(ggplot_build(
  ggplot(plot.df) + 
    geom_pointrange(aes(x = `exp(coef)`, y = y.pos, xmin = ll, xmax = ul, color = model.type, shape = shape), size = 0.7)  +
    # facet_grid(~set, scales = "free_y", space="free") +
    scale_color_manual(values = c("#00BFC4", "#F8766D", "darkgreen", "#19196F")) +
    scale_y_continuous(breaks = as.numeric(unique(plot.df$organ)), labels = unique(plot.df$organ)) +
    scale_shape_manual(values = c(15, 17)) +
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    # scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    # guides(color = "none") + # no legend for color
    # labs(fill="NES") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    # geom_text(size = 4.5, nudge_y = 0.13, na.rm = TRUE, colour = "black") +
    geom_vline(xintercept = 1, colour = "darkgrey") +
    xlab("Hazard ratio of mortality per standard deviation increase") +
    ylab("") +
    xlim(lims) +
    ggtitle(paste0("")) +
    # geom_vline(xintercept = 0) +
    # facet_grid(antibody~diet, scales = "free_y") +
    guides(color = guide_legend(order=1),
           shape = guide_legend(order=2)) + # First color, then shape in the legend
    theme_bw() +
    theme(legend.title = element_blank(), # No legend title
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          # plot.title = element_text(face = "bold"), # element_blank(),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"), # , angle = 45, hjust = 1, vjust = 1
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold",size = 12),
          strip.background = element_blank())
))

attr(HR.plot.mortality.other.papers, "data") <- plot.df

plot(HR.plot.mortality.other.papers)
