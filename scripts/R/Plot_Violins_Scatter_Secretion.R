
### From Plot_Volvanos_GSEA.R
res.age <- readRDS(file = paste0(rds.dir, "res_age.rds"))
### From Plot_Volvanos_GSEA.R
res.mortality <- readRDS(file = paste0(rds.dir, "res_mortality.rds"))
### From 3_Process_GTEx.R
Uniprot <- readRDS(file = paste0(rds.dir, "Uniprot.rds"))

### 1. Determine the localization of each protein, and save it as localization.df ###

gene.names <- strsplit(Uniprot$`Gene Names`, " ")

localization.indices <- vector(mode = "list", length = length(res.age$protein))
names(localization.indices) <- res.age$protein
proteins <- res.age$protein
proteins <- unlist(lapply(strsplit(proteins, "_"), function(x){return(x[1])}))
proteins[proteins == "ANP32C"] <- "ANP32CP"

for(i in 1:length(proteins)){
  localization.indices[[i]] <- which(unlist(lapply(gene.names, function(x){return(proteins[i] %in% x)})))
}

n <- unlist(lapply(localization.indices, length))
table(n)

names(localization.indices)[n == 0]
names(localization.indices)[n == 2]
Uniprot[localization.indices[["ADAM23"]],]

protein.categories <- vector(mode = "list", length = length(res.age$protein))
names(protein.categories) <- res.age$protein

# i <- 36
for(i in 1:length(localization.indices)){
  
  if(names(localization.indices)[i] != "NTproBNP"){
    
    if((any(grepl(": Secreted", Uniprot[localization.indices[[i]], "Subcellular location [CC]"], fixed = TRUE))) | (any(grepl("extracellular space ", Uniprot[localization.indices[[i]], "Gene Ontology (cellular component)"], fixed = TRUE)))){
      protein.categories[[i]] <- c(protein.categories[[i]], "secreted")
    }
    
    if((any(grepl(": Cell membrane", Uniprot[localization.indices[[i]], "Subcellular location [CC]"], fixed = TRUE))) | (any(grepl("plasma membrane ", Uniprot[localization.indices[[i]], "Gene Ontology (cellular component)"], fixed = TRUE)))){
      protein.categories[[i]] <- c(protein.categories[[i]], "cell membrane")
    }
    
  }
  
}

secreted <- unlist(lapply(protein.categories, function(x){
  return("secreted" %in% x)
}))

external <- unlist(lapply(protein.categories, function(x){
  return(("cell membrane" %in% x) & !("secreted" %in% x))
}))

internal <- unlist(lapply(protein.categories, function(x){
  return(is.null(x))
}))

all(res.mortality$protein == res.age$protein)
# TRUE # Important!!!

localization.df <- data.frame(protein = res.age$protein,
                              protein.name = res.age$protein.name,
                              zval.age = limma::zscoreT(res.age$`t value`, df = res.age$df, approx=FALSE, method = "bailey"), # qnorm(pt(abs(res.age$`t value`), df = res.age$df, lower.tail = FALSE, log.p = TRUE), lower.tail = FALSE, log.p = TRUE) * sign(res.age$`t value`)
                              zval.mortality = res.mortality$z,
                              secreted = secreted,
                              external = external,
                              internal = internal,
                              Significance = "not significant")

localization.df$Significance[res.age$qval < 0.05] <- "only age"
localization.df$Significance[res.mortality$qval < 0.05] <- "only mortality"
localization.df$Significance[(res.age$qval < 0.05) & (res.mortality$qval < 0.05)] <- "both"
localization.df$Significance <- factor(localization.df$Significance, levels = c("not significant", "only age", "only mortality", "both"))

# saveRDS(localization.df, file = paste0(rds.dir, "localization_df.rds"))
localization.df <- readRDS(file = paste0(rds.dir, "localization_df.rds"))

### 2. Make violin plot for chronological age adjusted for sex ###

table(rowSums(localization.df[, c("secreted", "external", "internal")]))
# 1 
# 2923

summary(lm(zval.age~secreted, data = localization.df))

summary(lm(zval.age~external, data = localization.df))

summary(lm(zval.age~internal, data = localization.df))

plot.df <- data.frame(
  protein = c(localization.df$protein[localization.df$secreted], localization.df$protein[localization.df$external], localization.df$protein[localization.df$internal]),
  protein.name = c(localization.df$protein.name[localization.df$secreted], localization.df$protein.name[localization.df$external], localization.df$protein.name[localization.df$internal]),
  zvalue = c(localization.df$zval.age[localization.df$secreted], localization.df$zval.age[localization.df$external], localization.df$zval.age[localization.df$internal]),
  localization = c(rep("Secreted", sum(localization.df$secreted)), rep("Extracellular membrane", sum(localization.df$external)), rep("Intracellular", sum(localization.df$internal)))
)

plot.df$localization <- factor(plot.df$localization, levels = c("Intracellular", "Extracellular membrane", "Secreted"))

summary <- summary(glht(lm(zvalue~localization, data = plot.df), linfct = mcp(localization = c("`Extracellular membrane` - Intracellular = 0", "Secreted - Intracellular = 0", "Secreted - `Extracellular membrane` = 0"))))
summary
# Simultaneous Tests for General Linear Hypotheses
# Multiple Comparisons of Means: User-defined Contrasts
# Fit: lm(formula = zvalue ~ localization, data = plot.df)
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
#   `Extracellular membrane` - Intracellular == 0   6.4388     0.9783   6.582 1.41e-10 ***
#   Secreted - Intracellular == 0                  14.6418     0.9323  15.704  < 1e-10 ***
#   Secreted - `Extracellular membrane` == 0        8.2030     1.0073   8.144  < 1e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

height1 <- (max(plot.df$zvalue, na.rm = TRUE) - min(plot.df$zvalue, na.rm = TRUE))*0.8+min(plot.df$zvalue, na.rm = TRUE)
height2 <- (max(plot.df$zvalue, na.rm = TRUE) - min(plot.df$zvalue, na.rm = TRUE))*0.97+min(plot.df$zvalue, na.rm = TRUE)

significance.df <- data.frame(x = c(1.5, 2, 2.5), 
                              xend = c(1.9, 2.9, 2.9), 
                              y = c(height1, height2, height1),
                              pval = summary$test$pvalues,
                              stars = c("n.s.", "n.s.", "n.s."),
                              size = c(4.5, 4.5, 4.5),
                              nudge_y = c(1, 1, 1)
)

significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
significance.df$stars[which(significance.df$pval < 0.001)] <- "***"

significance.df[(significance.df$stars != "n.s."), "size"] <- 10
significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- 5

violin.secreted.aging <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = localization, y = zvalue)) +
    geom_violin(aes(color = localization, fill = localization)) + # aes(color = Diet)
    
    scale_color_manual(values=c("#556219", "#418979", "#96410e")) +
    scale_fill_manual(values=c("#556219", "#418979", "#96410e")) +
    
    # facet_wrap(~protein, nrow = 1, strip.position = "bottom") +
    # geom_point(cex = 1, aes(x = jitter(as.numeric(Diet), factor = 1.5), color = Diet)) + # cex = 1.75
    
    # Too many points!
    # geom_quasirandom(cex = 1, aes(x = as.numeric(sex), color = sex, alpha = 1)) +
    
    geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
    guides(alpha = "none") + # no legend for "alpha"
    
    geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
    
    geom_segment(aes(x = significance.df$x[1]-0.4, y = significance.df$y[1], xend = significance.df$xend[1], yend = significance.df$y[1])) +
    geom_segment(aes(x = significance.df$x[2]-0.9, y = significance.df$y[2], xend = significance.df$xend[2], yend = significance.df$y[2])) +
    geom_segment(aes(x = significance.df$x[3]-0.4, y = significance.df$y[3], xend = significance.df$xend[3], yend = significance.df$y[3])) +
    
    # coord_cartesian(ylim = c(-10, 10)) +
    
    ggtitle("Chronological age adjusted for sex") +
    xlab("") + # No label on the horizontal axis
    ylab(paste0("Z value")) +
    theme_bw() +
    theme(# plot.title = element_text(face = "bold", hjust = 0.5), # centered and bold title
      legend.title = element_text(size = 12), # legend title size
      legend.text = element_text(size = 12), # legend text size
      legend.position = "none", # No legend
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size=12),
      strip.text.x = element_text(size=12), # the size of the facet labels
      axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
      axis.title.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line.x = element_line(color="black", linewidth = 0.5),
      axis.line.y = element_line(color="black", linewidth = 0.5),
      axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
      axis.ticks.y = element_line(color="black", linewidth = 0.5),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
))

attr(violin.secreted.aging, "data") <- plot.df

plot(violin.secreted.aging)


### 3. Make violin plot for mortality adjusted for age x sex interaction ###

table(rowSums(localization.df[, c("secreted", "external", "internal")]))
# 1 
# 2923

summary(lm(zval.mortality~secreted, data = localization.df))

summary(lm(zval.mortality~external, data = localization.df))

summary(lm(zval.mortality~internal, data = localization.df))

plot.df <- data.frame(
  protein = c(localization.df$protein[localization.df$secreted], localization.df$protein[localization.df$external], localization.df$protein[localization.df$internal]),
  protein.name = c(localization.df$protein.name[localization.df$secreted], localization.df$protein.name[localization.df$external], localization.df$protein.name[localization.df$internal]),
  zvalue = c(localization.df$zval.mortality[localization.df$secreted], localization.df$zval.mortality[localization.df$external], localization.df$zval.mortality[localization.df$internal]),
  localization = c(rep("Secreted", sum(localization.df$secreted)), rep("Extracellular membrane", sum(localization.df$external)), rep("Intracellular", sum(localization.df$internal)))
)

plot.df$localization <- factor(plot.df$localization, levels = c("Intracellular", "Extracellular membrane", "Secreted"))

summary(lm(zvalue~localization, data = plot.df))

summary <- summary(glht(lm(zvalue~localization, data = plot.df), linfct = mcp(localization = c("`Extracellular membrane` - Intracellular = 0", "Secreted - Intracellular = 0", "Secreted - `Extracellular membrane` = 0"))))
summary
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: lm(formula = zvalue ~ localization, data = plot.df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# `Extracellular membrane` - Intracellular == 0 -0.06709    0.08042  -0.834    0.682
# Secreted - Intracellular == 0                 -0.01841    0.07665  -0.240    0.969
# Secreted - `Extracellular membrane` == 0       0.04868    0.08280   0.588    0.826
# (Adjusted p values reported -- single-step method)

height1 <- (max(plot.df$zvalue, na.rm = TRUE) - min(plot.df$zvalue, na.rm = TRUE))*0.8+min(plot.df$zvalue, na.rm = TRUE)
height2 <- (max(plot.df$zvalue, na.rm = TRUE) - min(plot.df$zvalue, na.rm = TRUE))*0.97+min(plot.df$zvalue, na.rm = TRUE)

significance.df <- data.frame(x = c(1.5, 2, 2.5), 
                              xend = c(1.9, 2.9, 2.9), 
                              y = c(height1, height2, height1),
                              pval = summary$test$pvalues,
                              stars = c("n.s.", "n.s.", "n.s."),
                              size = c(4.5, 4.5, 4.5),
                              nudge_y = c(1, 1, 1)
)

significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
significance.df$stars[which(significance.df$pval < 0.001)] <- "***"

significance.df[(significance.df$stars != "n.s."), "size"] <- 10
significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- 5


violin.secreted.mortality <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = localization, y = zvalue)) +
    geom_violin(aes(color = localization, fill = localization)) + # aes(color = Diet)
    
    scale_color_manual(values=c("#556219", "#418979", "#96410e")) +
    scale_fill_manual(values=c("#556219", "#418979", "#96410e")) +
    
    geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
    guides(alpha = "none") + # no legend for "alpha"
    
    geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
    
    geom_segment(aes(x = significance.df$x[1]-0.4, y = significance.df$y[1], xend = significance.df$xend[1], yend = significance.df$y[1])) +
    geom_segment(aes(x = significance.df$x[2]-0.9, y = significance.df$y[2], xend = significance.df$xend[2], yend = significance.df$y[2])) +
    geom_segment(aes(x = significance.df$x[3]-0.4, y = significance.df$y[3], xend = significance.df$xend[3], yend = significance.df$y[3])) +

    ggtitle("Mortality adjusted for chronological age and sex") +
    xlab("") + # No label on the horizontal axis
    ylab(paste0("Z value")) +
    theme_bw() +
    theme(# plot.title = element_text(face = "bold", hjust = 0.5), # centered and bold title
      legend.title = element_text(size = 12), # legend title size
      legend.text = element_text(size = 12), # legend text size
      legend.position = "none", # No legend
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size=12),
      strip.text.x = element_text(size=12), # the size of the facet labels
      axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
      axis.title.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line.x = element_line(color="black", linewidth = 0.5),
      axis.line.y = element_line(color="black", linewidth = 0.5),
      axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
      axis.ticks.y = element_line(color="black", linewidth = 0.5),
      # strip.text = element_text(face = "black", linewidth = 12),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
))

attr(violin.secreted.mortality, "data") <- plot.df

plot(violin.secreted.mortality)


### 4. Make scatterplots for res.age and res.mortality for secreted and non-secreted proteins ###

all(res.age$protein == res.mortality$protein)
# TRUE

### 1. Non-secreted proteins ###

plot.df <- localization.df
plot.df <- plot.df[!plot.df$secreted,]

corrplot.nonsecreted <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = zval.age, y = zval.mortality, color = Significance)) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("black", "#04a3bd", "#f0be3d", "#247d3f"), 0.5)) +
    ggtitle(paste0("Non-secreted proteins (r² = ", round(summary(lm(plot.df[!plot.df$secreted,"zval.age"] ~ plot.df[!plot.df$secreted,"zval.mortality"]))$r.squared, 2), ")")) +
    xlab("Z value chronological age adj. for sex") +
    ylab("Z value mortality adj. for age x sex") +
    xlim(c(-80, 152)) +
    ylim(c(-5, 8.2)) +
    theme_bw() +
    theme(
      legend.position = "none", # No legend
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size=12),
      axis.title.x = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line.x = element_line(color="black", linewidth = 0.5),
      axis.line.y = element_line(color="black", linewidth = 0.5),
      axis.ticks.x = element_line(color="black", linewidth = 0.5),
      axis.ticks.y = element_line(color="black", linewidth = 0.5),
      strip.text = element_text(face = "black", size = 12),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
))

attr(corrplot.nonsecreted, "data") <- plot.df

plot(corrplot.nonsecreted)

sqrt(summary(lm(localization.df[localization.df$secreted,"zval.age"] ~ localization.df[localization.df$secreted,"zval.mortality"]))$r.squared)
# r = 0.52
sqrt(summary(lm(localization.df[!localization.df$secreted,"zval.age"] ~ localization.df[!localization.df$secreted,"zval.mortality"]))$r.squared)
# r = 0.34

### 2. Secreted proteins ###

plot.df <- localization.df
plot.df <- plot.df[plot.df$secreted,]

corrplot.secreted <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = zval.age, y = zval.mortality, color = Significance)) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("black", "#04a3bd", "#f0be3d", "#247d3f"), 0.5)) +
    ggtitle(paste0("Secreted proteins (r² = ", round(summary(lm(plot.df[plot.df$secreted,"zval.age"] ~ plot.df[plot.df$secreted,"zval.mortality"]))$r.squared, 2), ")")) +
    xlab("Z value chronological age adj. for sex") +
    ylab("Z value mortality adj. for age x sex") +
    theme_bw() +
    theme(
      legend.position = "none", # No legend
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size=12),
      axis.title.x = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.line.x = element_line(color="black", linewidth = 0.5),
      axis.line.y = element_line(color="black", linewidth = 0.5),
      axis.ticks.x = element_line(color="black", linewidth = 0.5),
      axis.ticks.y = element_line(color="black", linewidth = 0.5),
      strip.text = element_text(face = "black", size = 12),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
))

attr(corrplot.secreted, "data") <- plot.df

plot(corrplot.secreted)



