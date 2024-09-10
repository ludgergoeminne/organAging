
# From 1_Prepare_Data.R
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

# From 3_Process_GTEx.R
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

# From 4_calculate_predicted_resid_ages.R
predicted.ages.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_1.rds"))

### Plot the correlations with specific diseases ###

df <- data.frame(
  eid = olink_bd_annotation_0$eid,
  set = factor(ifelse((olink_bd_annotation_0$eid %in% training.eids), "Training set", "Test set"), levels = c("Training set", "Test set")), 
  `Chronological age` = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
  `Biological age` = predicted.ages.1$gen1$predicted[["Conventional"]],
  `Age deviation` = predicted.ages.1$gen1$residual[["Conventional"]],
  `Mother's age at death` = olink_bd_annotation_0$mothers_age_at_death,
  `Father's age at death` = olink_bd_annotation_0$fathers_age_at_death,
  `Grip strength strongest hand` = olink_bd_annotation_0$grip_strength_strongest_hand,
  `Grip strength weakest hand` = olink_bd_annotation_0$grip_strength_weakest_hand,
  check.names = FALSE
)

### Correlation heatmap ###

correlation.traits <- colnames(df[, -c(1, 2)], )

res.traits <- data.frame(trait.x = factor(rep(correlation.traits, each = length(correlation.traits)), levels = rev(correlation.traits)),
                      trait.y = factor(rep(correlation.traits, times = length(correlation.traits)), levels = correlation.traits),
                      set = "Diagonal",
                      Estimate = NA,
                      `Std. Error` = NA,
                      `t value` = NA, 
                      `Pr(>|t|)` = NA, 
                      r = NA, 
                      check.names = FALSE)

sum.diag <- length(correlation.traits)+1

res.traits[(as.numeric(res.traits$trait.x)+as.numeric(res.traits$trait.y)) < sum.diag, "set"] <- "Training set"
res.traits[(as.numeric(res.traits$trait.x)+as.numeric(res.traits$trait.y)) > sum.diag, "set"] <- "Test set"

for(i in 1:nrow(res.traits)){
  
  # Lower triangle and diagonal
  if(res.traits$set[i] %in% c("Diagonal", "Training set")){
    tmp <- df[df$set == "Training set", ]
    # Upper triangle
  } else  if(res.traits$set[i] == c("Test set")){
    tmp <- df[df$set == "Test set", ]
  }
  
  summary <- summary(lm(formula(get(as.character(res.traits$trait.x[i]))~get(as.character(res.traits$trait.y[i]))), data = tmp))
  res.traits[i, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")] <- summary$coef["get(as.character(res.traits$trait.y[i]))", ]
  res.traits[i, c("r")] <- sign(summary$coef["get(as.character(res.traits$trait.y[i]))", "Estimate"])*sqrt(summary$r.squared)
}

# Add q-value per group:
res.traits <- res.traits %>% group_by(set) %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))

# saveRDS(res.traits, file = paste0(rds.dir, "res_traits.rds"))
res.traits <- readRDS(file = paste0(rds.dir, "res_traits.rds"))

plot.df <- res.traits

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

# We also want to print the r values
plot.df$corr <- sprintf("%.2f", round(plot.df$r, 2))

plot.df$color <- "black"
plot.df$color[abs(plot.df$r) > 0.8] <- "white"

lims <- c(-1,1)

range(plot.df$r, na.rm = TRUE)
# -0.1904868  1.0000000

c(quantile(plot.df$r, 0.05, na.rm = TRUE), quantile(plot.df$r, 0.95, na.rm = TRUE))
# -0.1851483  1.000000 

# Only here change to factor!
# plot.df$trait.x = factor(plot.df$trait.x, levels = rev(unique(plot.df$trait.x)))
# plot.df$trait.y = factor(plot.df$trait.y, levels = unique(plot.df$trait.y))

# Turn the plot into a lower triangle
# plot.df <- plot.df[nrow(plot.df):1, ]
# plot.df <- plot.df[!duplicated(apply(apply(plot.df[, c("trait.x", "trait.y")], 1, function(x){sort(x)}), 2, function(x){return(paste0(x, collapse = ""))})),]

Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

heatmap.traits <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = trait.x, y = trait.y, label = stars)) + 
    geom_tile(aes(fill = r)) + # size = minuslog10p , color = p.adjust
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    guides(color = "none") + # no legend for color
    labs(fill="r") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    
    # geom_text(size = 4.5, nudge_y = -0.13, na.rm = TRUE) +
    # Exceptional because also text, not only stars:
    geom_text(color = plot.df$color, size = 4.5, nudge_y = 0, na.rm = TRUE) +
    geom_text(aes(label = corr), color = plot.df$color, size = 3, nudge_y = -0.175, na.rm = TRUE) +
    
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
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold",size = 12),
          strip.background = element_blank())
))

attr(heatmap.traits, "data") <- plot.df

# Lower triangle: training set, upper triangle: test set
plot(heatmap.traits)




