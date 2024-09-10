
### TODO: add script Alex here

# saveRDS(correlation.table.genes, file = paste0(rds.dir, "Correlation_table_genes.rds"))
correlation.table.genes <- readRDS(file = paste0(rds.dir, "Correlation_table_genes.rds"))

names(correlation.table.genes) <- c("r", "pval", "qval")

for(i in 1:length(correlation.table.genes)){
  tmp <- as.data.frame(correlation.table.genes[[i]])
  tmp$trait.x <- rownames(tmp)
  correlation.table.genes[[i]] <- pivot_longer(tmp, cols = 1:(ncol(tmp)-1), names_to = "trait.y", values_to = names(correlation.table.genes)[i])
  if(i != 1){
    correlation.table.genes[[i]]$trait.x <- NULL
    correlation.table.genes[[i]]$trait.y <- NULL
  }
}

# saveRDS(correlation.table.functions, file = paste0(rds.dir, "Correlation_table_functions.rds"))
correlation.table.functions <- readRDS(file = paste0(rds.dir, "Correlation_table_functions.rds"))

names(correlation.table.functions) <- c("r", "pval", "qval")

for(i in 1:length(correlation.table.functions)){
  tmp <- as.data.frame(correlation.table.functions[[i]])
  tmp$trait.x <- rownames(tmp)
  correlation.table.functions[[i]] <- pivot_longer(tmp, cols = 1:(ncol(tmp)-1), names_to = "trait.y", values_to = names(correlation.table.functions)[i])
  if(i != 1){
    correlation.table.functions[[i]]$trait.x <- NULL
    correlation.table.functions[[i]]$trait.y <- NULL
  }
}

### 1. Plot correlation table genes ###

plot.df <- do.call("cbind", unname(correlation.table.genes))
plot.df$trait.x <- gsub("ChronoAge", "chronological age", plot.df$trait.x)
plot.df$trait.x <- gsub("Hazard adjusted", "mortality", plot.df$trait.x)

plot.df$trait.y <- gsub("ChronoAge", "chronological age", plot.df$trait.y)
plot.df$trait.y <- gsub("Hazard adjusted", "mortality", plot.df$trait.y)

# Spearman's rho
plot.df$qval[plot.df$r == 1] <- 0
plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

# We also want to print the r values
plot.df$corr <- sprintf("%.2f", round(plot.df$r, 2))

plot.df$color <- "black"
plot.df$color[abs(plot.df$r) > 0.6] <- "white"

lims <- c(-1,1)

range(plot.df$r, na.rm = TRUE)
# 0.03814799  1.0000000

c(quantile(plot.df$r, 0.05, na.rm = TRUE), quantile(plot.df$r, 0.95, na.rm = TRUE))
# 0.04874184 1.00000000 

# Only here change to factor!
plot.df$trait.x = factor(plot.df$trait.x, levels = rev(unique(plot.df$trait.x)))
plot.df$trait.y = factor(plot.df$trait.y, levels = unique(plot.df$trait.y))

# Turn the plot into a lower triangle
plot.df <- plot.df[nrow(plot.df):1, ]
plot.df <- plot.df[!duplicated(apply(apply(plot.df[, c("trait.x", "trait.y")], 1, function(x){sort(x)}), 2, function(x){return(paste0(x, collapse = ""))})),]

Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

heatmap.genes <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = trait.x, y = trait.y, label = stars)) + 
    geom_tile(aes(fill = r)) + # size = minuslog10p , color = p.adjust
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    guides(color = "none") + # no legend for color
    labs(fill="Spearman's rho") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    
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
    theme(legend.position = c(-1, -0.5), # c(1, 0.8), # The coordinates for legend.position are x- and y- offsets from the bottom-left of the plot, ranging from 0 - 1
          axis.line = element_line(colour = "black"),
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

attr(heatmap.genes, "data") <- plot.df

plot(heatmap.genes)

### 2. Plot correlation table functions ###

plot.df <- do.call("cbind", unname(correlation.table.functions))
plot.df$trait.x <- gsub("ChronoAge", "chronological age", plot.df$trait.x)
plot.df$trait.x <- gsub("Hazard adjusted", "mortality", plot.df$trait.x)

plot.df$trait.y <- gsub("ChronoAge", "chronological age", plot.df$trait.y)
plot.df$trait.y <- gsub("Hazard adjusted", "mortality", plot.df$trait.y)

# Spearman's rho
plot.df$qval[plot.df$r == 1] <- 0
plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

# We also want to print the r values
plot.df$corr <- sprintf("%.2f", round(plot.df$r, 2))

plot.df$color <- "black"
plot.df$color[abs(plot.df$r) > 0.6] <- "white"

lims <- c(-1,1)

range(plot.df$r, na.rm = TRUE)
# 0.1072699 1.0000000

c(quantile(plot.df$r, 0.05, na.rm = TRUE), quantile(plot.df$r, 0.95, na.rm = TRUE))
# 0.1205661 1.0000000 

# Only here change to factor!
plot.df$trait.x = factor(plot.df$trait.x, levels = rev(unique(plot.df$trait.x)))
plot.df$trait.y = factor(plot.df$trait.y, levels = unique(plot.df$trait.y))

# Turn the plot into a lower triangle
plot.df <- plot.df[nrow(plot.df):1, ]
plot.df <- plot.df[!duplicated(apply(apply(plot.df[, c("trait.x", "trait.y")], 1, function(x){sort(x)}), 2, function(x){return(paste0(x, collapse = ""))})),]

Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

heatmap.functions <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = trait.x, y = trait.y, label = stars)) + 
    geom_tile(aes(fill = r)) + # size = minuslog10p , color = p.adjust
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    guides(color = "none") + # no legend for color
    labs(fill="Spearman's rho") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    
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
    theme(legend.position = c(-1, -0.5), # c(1, 0.8), # The coordinates for legend.position are x- and y- offsets from the bottom-left of the plot, ranging from 0 - 1
          axis.line = element_line(colour = "black"),
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

attr(heatmap.functions, "data") <- plot.df

plot(heatmap.functions)


