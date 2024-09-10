

### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from 4_calculate_predicted_resid_ages.R ###
predicted.ages.1 <- readRDS(paste0(rds.dir, "predicted_ages_1.rds"))
training.eids <- readRDS(file = paste0(rds.dir, "training_eids.rds"))
test.eids <- readRDS(file = paste0(rds.dir, "test_eids.rds"))

# Only needed for export to Excel
xlsx.quantiles <- c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99, 0.999, 1)

### 1. 1st-generation model: split per sex ###

plot.df <- data.frame(
  eid = olink_bd_annotation_0$eid,
  chronological_age = olink_bd_annotation_0$age_first_visit, 
  residual_age = predicted.ages.1$gen1$residual[["Conventional"]],
  sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men")),
  set = factor(ifelse((olink_bd_annotation_0$eid %in% training.eids), "Training set", "Test set"), levels = c("Training set", "Test set"))
)

### Statistics ###
t.test(residual_age~sex, data = plot.df[plot.df$set == "Training set", ], var.equal = TRUE)
# mean in group Women   mean in group Men 
# -0.02181390          0.02572406 
# p-value = 0.07572
t.test(residual_age~sex, data = plot.df[plot.df$set == "Test set", ], var.equal = TRUE)
# mean in group Women   mean in group Men 
# -0.02524861          0.02973807 
# p-value = 0.3278

table(plot.df[plot.df$set == "Training set", "sex"])
# Women   Men 
# 19460 16502 

table(plot.df[plot.df$set == "Test set", "sex"])
# Women   Men 
# 4862  4128 

violins.sex.gen1 <- vector(mode = "list", length = length(organ.proteins))
names(violins.sex.gen1) <- names(organ.proteins)

for(k in 1:length(violins.sex.gen1)){
  if(!(names(organ.proteins)[k] %in% c("Bladder"))){
    
    plot.df <- data.frame(
      eid = olink_bd_annotation_0$eid,
      chronological_age = olink_bd_annotation_0$age_first_visit, 
      residual_age = predicted.ages.1$gen1$residual[[k]],
      sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men")),
      set = factor(ifelse((olink_bd_annotation_0$eid %in% training.eids), "Training set", "Test set"), levels = c("Training set", "Test set"))
    )

    significance.df <- data.frame(x = rep(1.5, 2), 
                    y = rep(8, 2), 
                    set = factor(c("Training set", "Test set"), levels = c("Training set", "Test set")),
                    pval = c(NA, NA),
                    stars = c("n.s.", "n.s."),
                    size = c(4.5, 4.5),
                    nudge_y = c(1, 1)
                    )
    
    significance.df[significance.df$set == "Training set", ]$pval <- t.test(residual_age~sex, data = plot.df[plot.df$set == "Training set", ], var.equal = TRUE)$p.value
    significance.df[significance.df$set == "Test set", ]$pval <- t.test(residual_age~sex, data = plot.df[plot.df$set == "Test set", ], var.equal = TRUE)$p.value
    
    significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
    significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
    significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
    significance.df$stars[which(significance.df$pval < 0.001)] <- "***"
    
    significance.df[(significance.df$stars != "n.s."), "size"] <- 10
    significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- 0.13
    
    violins.sex.gen1[[k]] <- ggplot_gtable(ggplot_build(
      ggplot(plot.df, aes(x = sex, y = residual_age)) +
        facet_wrap(~set) +
        geom_violin(aes(color = sex, fill = sex)) + # aes(color = Diet)
        scale_color_manual(values=c("magenta", "#4169E1")) +
        scale_fill_manual(values=c("magenta", "#4169E1")) +
        # facet_wrap(~protein, nrow = 1, strip.position = "bottom") +
        # geom_point(cex = 1, aes(x = jitter(as.numeric(Diet), factor = 1.5), color = Diet)) + # cex = 1.75
        
        # Too many points!
        # geom_quasirandom(cex = 1, aes(x = as.numeric(sex), color = sex, alpha = 1)) +
        
        geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
        
        # geom_text(aes(x, y, label = stars), data = significance.df, size = 10, nudge_y = 0.13, na.rm = TRUE, colour = "black") +
        geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
        
        geom_segment(aes(x = 1.1, y = 8, xend = 1.9, yend = 8)) +
        
        guides(alpha = "none") + # no legend for "alpha"
        coord_cartesian(ylim = c(-10, 10)) +
        ggtitle("1st-generation model") +
        xlab("") + # No label on the horizontal axis
        ylab(paste0("Age deviation")) +
        theme_bw() +
        theme(# plot.title = element_text(face = "bold", hjust = 0.5), # centered and bold title
          legend.title = element_text(size = 12), # legend title size
          legend.text = element_text(size = 12), # legend text size
          legend.position = "none", # No legend
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=12),
          # strip.text.x = element_text(size=12), # the size of the facet labels
          axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 12),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    ))
    
    tmp <- data.frame(
      Set = rep(c("Training set", "Test set"), each = 2*length(xlsx.quantiles)),
      Sex = rep(c(rep("Women", each = length(xlsx.quantiles)), rep("Men", each = length(xlsx.quantiles))), 2),
      Quantile = rep(paste0(100*xlsx.quantiles, "%"), times = 4)
    )
    
    tmp[, "Age deviation"] <- c(
      round(quantile(plot.df[(plot.df$set == "Training set") & (plot.df$sex == "Women"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$set == "Training set") & (plot.df$sex == "Men"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$set == "Test set") & (plot.df$sex == "Women"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$set == "Test set") & (plot.df$sex == "Men"), "residual_age"], xlsx.quantiles), 2)
    )
    
    attr(violins.sex.gen1[[k]], "data") <- tmp
    
  }
}


plot(violins.sex.gen1$Conventional)

### 2. Mortality-based model: split per sex ###

plot.df <- data.frame(
  eid = olink_bd_annotation_0$eid,
  chronological_age = olink_bd_annotation_0$age_first_visit, 
  residual_age = predicted.ages.1$gen2$residual[["Conventional"]],
  sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men")),
  set = factor(ifelse((olink_bd_annotation_0$eid %in% training.eids), "Training set", "Test set"), levels = c("Training set", "Test set"))
)

### Statistics ###
t.test(residual_age~sex, data = plot.df[plot.df$set == "Training set", ], var.equal = TRUE)
# mean in group Women   mean in group Men 
# -0.1718089           0.2026058 
# p-value < 2.2e-16
t.test(residual_age~sex, data = plot.df[plot.df$set == "Test set", ], var.equal = TRUE)
# mean in group Women   mean in group Men 
# -0.1626376           0.1915562 
# p-value < 2.2e-16

violins.sex.gen2 <- vector(mode = "list", length = length(organ.proteins))
names(violins.sex.gen2) <- names(organ.proteins)

for(k in 1:length(violins.sex.gen2)){
  if(!(names(organ.proteins)[k] %in% c("Bladder"))){
    
    plot.df <- data.frame(
      eid = olink_bd_annotation_0$eid,
      chronological_age = olink_bd_annotation_0$age_first_visit, 
      residual_age = predicted.ages.1$gen2$residual[[k]],
      sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men")),
      set = factor(ifelse((olink_bd_annotation_0$eid %in% training.eids), "Training set", "Test set"), levels = c("Training set", "Test set"))
    )
    
    height <- (max(plot.df$residual_age, na.rm = TRUE) - min(plot.df$residual_age, na.rm = TRUE))*0.9+min(plot.df$residual_age, na.rm = TRUE)
    
    significance.df <- data.frame(x = rep(1.5, 2), 
                                  y = rep(height, 2), 
                                  set = factor(c("Training set", "Test set"), levels = c("Training set", "Test set")),
                                  pval = c(NA, NA),
                                  stars = c("n.s.", "n.s."),
                                  size = c(4.5, 4.5),
                                  nudge_y = c(1, 1)
    )
    
    significance.df[significance.df$set == "Training set", ]$pval <- t.test(residual_age~sex, data = plot.df[plot.df$set == "Training set", ])$p.value
    significance.df[significance.df$set == "Test set", ]$pval <- t.test(residual_age~sex, data = plot.df[plot.df$set == "Test set", ])$p.value
    
    significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
    significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
    significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
    significance.df$stars[which(significance.df$pval < 0.001)] <- "***"
    
    significance.df[(significance.df$stars != "n.s."), "size"] <- 10
    significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- 0.13
    
    violins.sex.gen2[[k]] <- ggplot_gtable(ggplot_build(
      ggplot(plot.df, aes(x = sex, y = residual_age)) +
        geom_violin(aes(color = sex, fill = sex)) + # aes(color = Diet)
        facet_wrap(~set) +
        scale_color_manual(values=c("magenta", "#4169E1")) +
        scale_fill_manual(values=c("magenta", "#4169E1")) +
        # facet_wrap(~protein, nrow = 1, strip.position = "bottom") +
        # geom_point(cex = 1, aes(x = jitter(as.numeric(Diet), factor = 1.5), color = Diet)) + # cex = 1.75
        
        # Too many points!
        # geom_quasirandom(cex = 1, aes(x = as.numeric(sex), color = sex, alpha = 1)) +
        
        geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
        
        geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
        
        geom_segment(aes(x = 1.1, y = height, xend = 1.9, yend = height)) +
        
        guides(alpha = "none") + # no legend for "alpha"
        # coord_cartesian(ylim = c(-10, 10)) +
        ggtitle("Mortality-based model") +
        xlab("") + # No label on the horizontal axis
        ylab(paste0("Log(mortality hazard) deviation")) +
        theme_bw() +
        theme(# plot.title = element_text(face = "bold", hjust = 0.5), # centered and bold title
          legend.title = element_text(size = 12), # legend title size
          legend.text = element_text(size = 12), # legend text size
          legend.position = "none", # No legend
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size=12),
          # strip.text.x = element_text(size=12), # the size of the facet labels
          axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 12),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    ))
    
    tmp <- data.frame(
      Set = rep(c("Training set", "Test set"), each = 2*length(xlsx.quantiles)),
      Sex = rep(c(rep("Women", each = length(xlsx.quantiles)), rep("Men", each = length(xlsx.quantiles))), 2),
      Quantile = rep(paste0(100*xlsx.quantiles, "%"), times = 4)
    )
    
    tmp[, "Log(mortality hazard) deviation"] <- c(
      round(quantile(plot.df[(plot.df$set == "Training set") & (plot.df$sex == "Women"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$set == "Training set") & (plot.df$sex == "Men"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$set == "Test set") & (plot.df$sex == "Women"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$set == "Test set") & (plot.df$sex == "Men"), "residual_age"], xlsx.quantiles), 2)
    )
    
    attr(violins.sex.gen2[[k]], "data") <- tmp
    
  }
}


plot(violins.sex.gen2$Conventional)




