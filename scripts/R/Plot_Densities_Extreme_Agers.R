

### This comes from 4_Calculate_Predicted_Residual_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

### Based on mortality-based model ###

pvals <- 2*pnorm(-abs(predicted.ages.oof$gen2$residual[["Conventional"]]), 0, sd(predicted.ages.oof$gen2$residual[["Conventional"]]))
hist(pvals, breaks = 100)
qvals <- p.adjust(pvals, method = "BH")
range(predicted.ages.oof$gen2$residual[["Conventional"]][qvals < 0.05])
# For residuals full dataset:
# 2.680542 5.602474 # The extreme cases all have positive residuals

plot.df <- data.frame(eid = olink_bd_annotation_0$eid,
                      n.diseases = rowSums(!is.na(olink_bd_annotation_0[, grepl("^p13", colnames(olink_bd_annotation_0))])),
                      Group = ifelse(qvals < 0.05, "Extreme", "Rest")
                      )

### Statistics ###
mu <- plot.df %>% 
  group_by(Group) %>%
  summarise(grp.mean = mean(n.diseases)) %>%
  as.data.frame
mu
# Group grp.mean
# 1 Extreme 70.77882
# 2    Rest 39.10410

table(plot.df$Group)
# Extreme    Rest 
# 321   44631 

mu$txt.mean <- sprintf("%.1f", round(mu$grp.mean, 2))
mu$nudge_y <- rep(c(-0.007, 0.003))

significance.df <- data.frame(x = mu %>% dplyr::summarize(mean = mean(grp.mean, na.rm=TRUE)) %>%  `$`("mean"), 
                              y = rep(0.020, 2), 
                              pval = c(NA, NA),
                              stars = c("n.s.", "n.s."),
                              size = c(5, 5),
                              nudge_y = c(0, 0)
)

significance.df$pval <- t.test(n.diseases ~ Group, data = plot.df, var.equal = TRUE)$p.value
significance.df$pval
# 6.95807e-139 6.95807e-139

significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
significance.df$stars[which(significance.df$pval < 0.001)] <- "***"

data <- ggplot_build(
  ggplot(plot.df, aes(x = n.diseases)) + geom_density(aes(fill = Group), alpha = 0.4) +
    geom_vline(aes(xintercept = grp.mean, color = Group), data = mu, linetype = "dashed") +
    scale_color_manual(values = c("#EFC000FF", "#868686FF"))+
    scale_fill_manual(values = c("#EFC000FF", "#868686FF")) +
    guides(color = "none") + # No legend for color, but keep legend for fill
    # geom_vline(xintercept = 0) +
    # facet_grid(antibody~diet, scales = "free_y") +
    
    geom_text(aes(x = grp.mean, y = 0.020, label = txt.mean), nudge_y = mu$nudge_y, size = 4, data = mu, na.rm = TRUE, colour = "black") +
    geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
    # geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
    
    xlab("Number of disease events") +
    ylab("Density") +
    theme_bw() +
    theme(legend.position = c(0.9, 0.5), # The coordinates for legend.position are x- and y- offsets from the bottom-left of the plot, ranging from 0 - 1
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          plot.title = element_blank(), # element_text(face = "bold"), 
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"), #, angle = 45, hjust = 1, vjust = 1
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold",size = 12),
          strip.background = element_blank())
)
density.plot.gen2 <- ggplot_gtable(data)

attr(density.plot.gen2, "data") <- data$data[[1]]

plot(density.plot.gen2)



