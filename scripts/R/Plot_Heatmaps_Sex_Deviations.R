
### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(file = paste0(rds.dir, "organ_proteins.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

### 1. 1st-generation model ###

df.sex.deviations.gen1 <- data.frame(
  organ = names(organ.proteins),
  `DELTA age deviation` = rep(NA, length(organ.proteins)),
  pval = rep(NA, length(organ.proteins)), 
  check.names = FALSE
)

for(k in 1:length(organ.proteins)){
  if(names(organ.proteins)[k] != "Bladder"){
  
  df <- data.frame(
    eid = olink_bd_annotation_0$eid,
    chronological_age = olink_bd_annotation_0$age_first_visit, 
    residual_age = predicted.ages.oof$gen1$residual[[k]],
    sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men"))
  )
  
  ### Statistics ###
  t.stat <- t.test(residual_age~sex, data = df, var.equal = TRUE)
  
  df.sex.deviations.gen1[k, "DELTA age deviation"] <- t.stat$estimate[2]-t.stat$estimate[1]
  df.sex.deviations.gen1[k, "pval"] <- t.stat$p.value
  
  }
}






### 2. mortality-based model ###

df.sex.deviations.gen2 <- data.frame(
  organ = names(organ.proteins),
  `DELTA log(mortality hazard) \ndeviation` = rep(NA, length(organ.proteins)),
  pval = rep(NA, length(organ.proteins)), 
  check.names = FALSE
)

for(k in 1:length(organ.proteins)){
  if(names(organ.proteins)[k] != "Bladder"){
    
    df <- data.frame(
      eid = olink_bd_annotation_0$eid,
      chronological_age = olink_bd_annotation_0$age_first_visit, 
      residual_age = predicted.ages.oof$gen2$residual[[k]],
      sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men"))
    )
    
    ### Statistics ###
    t.stat <- t.test(residual_age~sex, data = df, var.equal = TRUE)
    
    df.sex.deviations.gen2[k, "DELTA log(mortality hazard) \ndeviation"] <- t.stat$estimate[2]-t.stat$estimate[1]
    df.sex.deviations.gen2[k, "pval"] <- t.stat$p.value
    
  }
}

### 3. Make the plot for the 1st-generation model ###

plot.df <- df.sex.deviations.gen1[df.sex.deviations.gen1$organ %in% names(organ.proteins.selected),]
plot.df$qval <- p.adjust(plot.df$pval, method = "hommel")

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

Heatmap_palette <- c("magenta", "white","#4169E1")

range(plot.df$`DELTA age deviation`)
lims <- c(-2, 2)

plot.df$x <- factor(1)
plot.df$organ <- factor(plot.df$organ, levels = rev(names(organ.proteins.selected)))

heatmap.sex.deviations.gen1 <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = x, y = organ, label = stars)) + 
    geom_tile(aes(fill = `DELTA age deviation`)) + # size = minuslog10p , color = p.adjust
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    labs(fill="DELTA age deviation") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    geom_text(size = 4.5, nudge_y = -0.13, na.rm = TRUE) +
    xlab("") +
    ylab("") +
    ggtitle("") +
    # geom_vline(xintercept = 0) +
    # facet_grid(antibody~diet, scales = "free_y") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          plot.title = element_blank(), # element_text(face = "bold"), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.line.y = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_blank(), # element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold",size = 12),
          strip.background = element_blank(),
          plot.margin = unit(c(0, 0, 1, 0), # EXCEPTIONAL: increase bottom margin
                             "inches"))
))

attr(heatmap.sex.deviations.gen1, "data") <- plot.df

plot(heatmap.sex.deviations.gen1)

### 4. Make the plot for the mortality-based model ###

plot.df <- df.sex.deviations.gen2[df.sex.deviations.gen2$organ %in% names(organ.proteins.selected),]
plot.df$qval <- p.adjust(plot.df$pval, method = "hommel")

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

Heatmap_palette <- c("magenta", "white","#4169E1")

range(plot.df$`DELTA log(mortality hazard) \ndeviation`)
lims <- c(-0.5, 0.5)

plot.df$x <- factor(1)
plot.df$organ <- factor(plot.df$organ, levels = rev(names(organ.proteins.selected)))

heatmap.sex.deviations.gen2 <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = x, y = organ, label = stars)) + 
    geom_tile(aes(fill = `DELTA log(mortality hazard) \ndeviation`)) + # size = minuslog10p , color = p.adjust
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    labs(fill="DELTA log(mortality hazard) \ndeviation") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    geom_text(size = 4.5, nudge_y = -0.13, na.rm = TRUE) +
    xlab("") +
    ylab("") +
    ggtitle("") +
    # geom_vline(xintercept = 0) +
    # facet_grid(antibody~diet, scales = "free_y") +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          plot.title = element_blank(), # element_text(face = "bold"), 
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.line.y = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_blank(), # element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold",size = 12),
          strip.background = element_blank(),
          plot.margin = unit(c(0, 0, 1, 0), # EXCEPTIONAL: increase bottom margin
                             "inches"))
))

tmp <- plot.df
colnames(tmp)[2] <- "DELTA log(mortality hazard) deviation"
attr(heatmap.sex.deviations.gen2, "data") <- tmp

plot(heatmap.sex.deviations.gen2)


