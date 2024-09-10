
# From 1_Prepare_Data.R
coding143 <- readRDS(paste0(rds.dir, "coding143.rds"))
olink_data_unfiltered_0 <- readRDS(file = paste0(rds.dir, "olink_data_unfiltered_0.rds"))
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

# From 3_Process_GTEx.R
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### Dimensions ###
dim(olink_data_unfiltered_0)
# 53014  2923

# Here, we want to do the analysis for all proteins, so no filtering here

res.age <- data.frame(
  protein.code = colnames(olink_data_unfiltered_0),
  Estimate = NA,
  `Std. Error` = NA,
  df = NA,
  `t value` = NA,
  `Pr(>|t|)` = NA, 
  ll = NA, 
  ul = NA, 
  r = NA,
  check.names = FALSE)

df <- cbind(data.frame(
  eid = bd[rownames(olink_data_unfiltered_0),]$eid,
  chronological_age = bd[rownames(olink_data_unfiltered_0),]$age_first_visit, # == df.test$age_first_visit,
  sex = bd[rownames(olink_data_unfiltered_0),]$p31,
  COPD = as.numeric(!is.na(bd[rownames(olink_data_unfiltered_0),]$p131492))
), olink_data_unfiltered_0)

### Statistics associations SCGB1A1 - COPD ###

summary(lm(COPD ~  SCGB1A1 + chronological_age * sex, data = df))
summary(lm(COPD ~  SCGB1A1 + sex, data = df))
summary(lm(COPD ~  SCGB1A1, data = df))
summary(lm(chronological_age ~  SCGB1A1 + sex, data = df))

for(i in 1:length(res.age$protein.code)){
  model <- lm(chronological_age ~  get(res.age$protein.code[i]) + sex, data = df)
  summary <- summary(model)
  res.age[i, c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary$coef["get(res.age$protein.code[i])", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], confint(model, parm = "get(res.age$protein.code[i])"), sqrt(summary$r.squared)*sign(model$coefficients["get(res.age$protein.code[i])"]))
  res.age[i, c("df")] <- model$df.residual
}

res.age <- dplyr::right_join(coding143, res.age, by = c("protein" = "protein.code"))

# Add q-value per group:
res.age <- res.age %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))

# saveRDS(res.age, file = paste0(rds.dir, "res_age.rds"))
res.age <- readRDS(file = paste0(rds.dir, "res_age.rds"))

### Volcano plot chronological age ###

dim(res.age[which((res.age$Estimate > 0) & (res.age$qval < 0.05)),])
# 1659   12
dim(res.age[which((res.age$Estimate < 0) & (res.age$qval < 0.05)),])
# 647  12

plot.df <- res.age
plot.df$Significance <- ifelse(plot.df$qval < 0.05, "FDR < 5%", "FDR >= 5%")
plot.df <- na.omit(plot.df)

plot.df$zval <- limma::zscoreT(plot.df$`t value`, df = plot.df$df, approx=FALSE, method = "bailey")

plot.df$color <- 0
plot.df$color[(plot.df$zval < 0) & (plot.df$qval < 0.05)] <- 1
plot.df$color[(plot.df$zval > 0) & (plot.df$qval < 0.05)] <- 2
plot.df$color <- as.factor(plot.df$color)

# Those are the same!
# 2*pt(abs(3.730877e-03), df = 51274, lower.tail = FALSE)
# 10^((log(2)+pt(abs(3.730877e-03), df = 51274, lower.tail = FALSE, log.p = TRUE))/log(10))
### Add log10p to avoid 0 values! ###
plot.df$minuslog10p <- -((log(2)+pt(abs(plot.df$`t value`), df = plot.df$df, lower.tail = FALSE, log.p = TRUE))/log(10))

volcano.age <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = Estimate, y = minuslog10p, color = color)) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("black", "#ff0000ff", "#00be64ff"), 0.5)) +
    ggtitle("") + # ggtitle("Chronological age adjusted for sex") +
    xlab("Estimated effect on chronological age") +
    ylab("-log10(p-value)") +
    theme_bw() +
    theme(# plot.title = element_text(face = "bold", hjust = 0.5), # centered and bold title
      # legend.title = element_text(size = 12), # legend title size
      # legend.text = element_text(size = 12), # legend text size
      legend.position = "none", # No legend
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size=12),
      axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
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

attr(volcano.age, "data") <- plot.df

plot(volcano.age)



### fgsea analysis ###

species <- "Homo sapiens"
m_df_bp <- msigdbr(species = species, category = "C5",subcategory = "BP")
m_df_kegg <- msigdbr(species = species, category = "C2",subcategory = "CP:KEGG")
m_df_reactome <- msigdbr(species = species, category = "C2",subcategory = "CP:REACTOME")
m_df_hall <- msigdbr(species = species, category = "H")
m_df <- rbind(rbind(m_df_bp,m_df_kegg),rbind(m_df_reactome,m_df_hall))
functions_list <- levels(factor(m_df$gs_name))
geneset_annotation <- list()
for (i in functions_list){
  temp_data <- m_df[m_df$gs_name==i,]
  temp_data <- as.character(temp_data$gene_symbol)
  geneset_annotation[[i]] <- temp_data
}

# saveRDS(geneset_annotation, file = paste0(rds.dir, "geneset_annotation.rds"))
geneset_annotation <- readRDS(file = paste0(rds.dir, "geneset_annotation.rds"))

### 1. Msigdbr analysis for chronological age ###

zvalues <- limma::zscoreT(res.age$`t value`, df = res.age$df, approx=FALSE, method = "bailey")
names(zvalues) <- res.age$protein
zvalues <- na.omit(zvalues)
zvalues[((is.infinite(zvalues)) & zvalues > 0)] <- 1e99
zvalues[((is.infinite(zvalues)) & zvalues < 0)] <- -1e99
# nPermSimple = 100000 to avoid pathways for which p-values are not properly calculated
# Now three such pathways for which log2err = NA
fgsea.age <- fgsea(geneset_annotation, zvalues, minSize = 0, maxSize = Inf, eps = 0, nPermSimple = 100000) # if not parallellisation: , nproc = 1
fgsea.age$zval <- sign(fgsea.age$NES)*qnorm(fgsea.age$pval/2, lower.tail = FALSE)
fgsea.age <- fgsea.age %>% arrange(pval)

# saveRDS(fgsea.age, file = paste0(rds.dir, "fgsea_age.rds"))
fgsea.age <- readRDS(file = paste0(rds.dir, "fgsea_age.rds"))

plot.df <- fgsea.age
plot.df <- plot.df %>% arrange(pval)

top.pathways <- c(plot.df[plot.df$NES < 0, ]$pathway[1:5], 
                  plot.df[plot.df$NES > 0, ]$pathway[1:5])

plot.df <- plot.df[plot.df$pathway %in% top.pathways, ]

plot.df$stars <- ""
plot.df$stars[plot.df$padj < 0.1] <- "'"
plot.df$stars[plot.df$padj < 0.05] <- "*"
plot.df$stars[plot.df$padj < 0.01] <- "**"
plot.df$stars[plot.df$padj < 0.001] <- "***"

plot.df$GO_description <- as.character(plot.df$pathway)
plot.df$GO_description <- paste(toupper(substr(plot.df$GO_description, 1, 1)), tolower(substr(plot.df$GO_description, 2, nchar(plot.df$GO_description))), sep="") # tolower(plot.df$GO_description)
plot.df$GO_description <- gsub("_", " ", plot.df$GO_description)
# plot.df$GO_description <- paste0(gsub("%GO[a-z]{2}%GO:", " (GO:", plot.df$GO_description), ")")
plot.df$GO_description <- gsub("MRNA ", "mRNA ", plot.df$GO_description)
plot.df$GO_description <- gsub("kegg", "KEGG", plot.df$GO_description)
plot.df$GO_description <- gsub("Kegg ", "KEGG: ", plot.df$GO_description)
plot.df$GO_description <- gsub("Gobp ", "GObp: ", plot.df$GO_description)
plot.df$GO_description <- gsub("Reactome ", "Reactome: ", plot.df$GO_description)
plot.df$GO_description <- gsub("p450", "P450", plot.df$GO_description)
plot.df$GO_description <- gsub(" system kks", " system (KKS)", plot.df$GO_description)
plot.df$GO_description <- gsub(" system cas", " system (CAS)", plot.df$GO_description)
plot.df$GO_description <- gsub("signaling by kit in", "signaling by KIT in", plot.df$GO_description)
plot.df$GO_description <- gsub("dna ", "DNA ", plot.df$GO_description)
plot.df$GO_description <- gsub("rhoj gtpase", "RHOJ GTPase", plot.df$GO_description)

# Sort by p-value (most significant first)
plot.df <- plot.df %>% arrange(pval)
# Then sort by sign; is not the same as sorting on z-value, which will put lowest z-values last!
plot.df <- plot.df %>% arrange(-sign(NES))
# To keep order in plot!
plot.df$GO_description <- factor(plot.df$GO_description, levels = rev(plot.df$GO_description))
# ### Overrule, take custom order: ###
# plot.df$GO_description <- factor(plot.df$GO_description, levels = rev(c("Mitochondrion (GO:0005739)",
#                                                                   "Lipolysis in adipocytes (KEGG: mmu04923)",
#                                                                   "Structural constituent of muscle (GO:0008307)",
#                                                                   "mRNA splicing, via spliceosome (GO:0000398)")))

hjust <- rep(NA, nrow(plot.df))
hjust[plot.df$NES > 0] <- 1.6
hjust[plot.df$NES < 0] <- -0.7

x_pos_stars <- rep(NA, nrow(plot.df))
x_pos_stars[plot.df$NES > 0] <- 2.7
x_pos_stars[plot.df$NES < 0] <- -2.7

rescale_colors <- function(x){
  x <- c(0,x,1)
  out <- scales::rescale(log(x+0.0001))
  return(out[-c(1,length(out))])
}
labels_scale <- c(0, 0.001, 0.01, 0.05, 0.1, 1)

# To make sure that a zero after the comma is also printed out
numbers <- unlist(rapply(list(plot.df$NES), sprintf, fmt = "%0.1f", how = "replace"))

main = paste0("")
xlim = c(-3, 3)

barplot.fgsea.age <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = NES, y = GO_description)) + # 1:10
    ggtitle(main) +
    ylab("") + # "Gene set"
    xlab("Normalized Enrichment Score") +
    geom_bar(stat="identity", aes(fill = rescale_colors(padj))) + # fill = colors
    scale_fill_viridis(option = "plasma", name = "q-value", limits = c(0, 1), 
                       values = rescale_colors(labels_scale), 
                       breaks = rescale_colors(labels_scale), 
                       labels = labels_scale,
                       direction = -1,
                       guide = guide_colorbar(reverse = TRUE)) +
    geom_text(aes(label = numbers), hjust = hjust, color="white", size=3, fontface = "bold") + # Not useful, all q-values of 0.01
    geom_text(aes(label = stars), x = x_pos_stars, color="black", size=5) + # , fontface = "bold"
    geom_vline(xintercept = 0) +
    xlim(xlim[1], xlim[2]) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size=12),
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12, face="bold"), # smaller title, bold
          legend.title = element_text(size = 10, face = "bold")
    )
))

attr(barplot.fgsea.age, "data") <- as.data.frame(plot.df)

plot(barplot.fgsea.age)

### Volcano plot statistics mortality ###

res.mortality <- data.frame(
  protein.code = colnames(olink_data_allproteins_0),
  `exp(coef)` = NA,
  coef = NA,
  z = NA,
  `Pr(>|z|)` = NA, 
  ll = NA, 
  ul = NA, 
  check.names = FALSE)

df <- data.frame(
  eid = olink_bd_annotation_0$eid,
  status_os = ifelse(is.na(olink_bd_annotation_0[, "mortality"]), 0, 1),
  time = olink_bd_annotation_0[, "mortality"],
  chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
  sex = olink_bd_annotation_0$p31
)

df$time <- as.numeric(as.Date(df$time, na.rm = TRUE))-as.numeric(as.Date(olink_bd_annotation_0[,"p53_i0"]))
df$time[df$time < 0] <- 0 # The ones that already had dementia, we put at 0
df$time[is.na(df$time)] <- (latest.date[i]-as.numeric(as.Date(olink_bd_annotation_0[,"p53_i0"])))[is.na(df$time)]
df <- cbind(df, olink_data_allproteins_0)

for(i in 1:length(res.mortality$protein.code)){
  cox.model <- coxph(Surv(time = time, event = status_os, type='right') ~ chronological_age * sex + get(res.mortality$protein.code[i]), data = df)
  # summary(cox.model)
  res.mortality[i, c("exp(coef)", "coef", "z", "Pr(>|z|)", "ll", "ul")] <- c(summary(cox.model)$coef["get(res.mortality$protein.code[i])", c("exp(coef)", "coef", "z", "Pr(>|z|)")], exp(confint(cox.model, parm = "get(res.mortality$protein.code[i])")))
}

res.mortality <- dplyr::right_join(coding143, res.mortality, by = c("protein" = "protein.code"))

# Add q-value per group:
res.mortality <- res.mortality %>% mutate(qval = p.adjust(`Pr(>|z|)`, method = "BH"))

# saveRDS(res.mortality, file = paste0(rds.dir, "res_mortality.rds"))
res.mortality <- readRDS(file = paste0(rds.dir, "res_mortality.rds"))


nlabels <- 10

plot.df <- res.mortality
plot.df$Significance <- ifelse(plot.df$qval < 0.05, "FDR < 5%", "FDR >= 5%")
plot.df <- na.omit(plot.df)

plot.df$color <- 0
plot.df$color[(plot.df$z < 0) & (plot.df$qval < 0.05)] <- 1
plot.df$color[(plot.df$z > 0) & (plot.df$qval < 0.05)] <- 2
plot.df$color <- as.factor(plot.df$color)

volcano.mortality <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = coef, y = -log10(`Pr(>|z|)`), color = color)) +
    geom_point(cex = 2.5) +
    scale_color_manual(values = alpha(c("black", "#ff0000ff", "#00be64ff"), 0.5)) +
    ggtitle("") + # ggtitle("Mortality") +
    xlab("Estimated log(mortality hazard ratio)") +
    ylab("-log10(p-value)") +
    theme_bw() +
    theme(# plot.title = element_text(face = "bold", hjust = 0.5), # centered and bold title
      # legend.title = element_text(size = 12), # legend title size
      # legend.text = element_text(size = 12), # legend text size
      legend.position = "none", # No legend
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size=12),
      axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
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

attr(volcano.mortality, "data") <- plot.df

plot(volcano.mortality)

### 1. Msigdbr analysis for mortality ###

zvalues <-  res.mortality$z
names(zvalues) <- res.mortality$protein
zvalues <- na.omit(zvalues)
zvalues[((is.infinite(zvalues)) & zvalues > 0)] <- 1e99
zvalues[((is.infinite(zvalues)) & zvalues < 0)] <- -1e99
fgsea.mortality <- fgsea(geneset_annotation, zvalues, minSize = 0, maxSize = Inf, eps = 0) # if not parallellisation: , nproc = 1
fgsea.mortality$zval <- sign(fgsea.mortality$NES)*qnorm(fgsea.mortality$pval/2, lower.tail = FALSE)
fgsea.mortality <- fgsea.mortality %>% arrange(pval)

# saveRDS(fgsea.mortality, file = paste0(rds.dir, "fgsea_mortality.rds"))
fgsea.mortality <- readRDS(file = paste0(rds.dir, "fgsea_mortality.rds"))



plot.df <- fgsea.mortality
plot.df <- plot.df %>% arrange(pval)

top.pathways <- c(plot.df[plot.df$NES < 0, ]$pathway[1:5], 
                  plot.df[plot.df$NES > 0, ]$pathway[1:5])

plot.df <- plot.df[plot.df$pathway %in% top.pathways, ]

plot.df$stars <- ""
plot.df$stars[plot.df$padj < 0.1] <- "'"
plot.df$stars[plot.df$padj < 0.05] <- "*"
plot.df$stars[plot.df$padj < 0.01] <- "**"
plot.df$stars[plot.df$padj < 0.001] <- "***"

plot.df$GO_description <- as.character(plot.df$pathway)
plot.df$GO_description <- paste(toupper(substr(plot.df$GO_description, 1, 1)), tolower(substr(plot.df$GO_description, 2, nchar(plot.df$GO_description))), sep="") # tolower(plot.df$GO_description)
plot.df$GO_description <- gsub("_", " ", plot.df$GO_description)
# plot.df$GO_description <- paste0(gsub("%GO[a-z]{2}%GO:", " (GO:", plot.df$GO_description), ")")
plot.df$GO_description <- gsub("MRNA ", "mRNA ", plot.df$GO_description)
plot.df$GO_description <- gsub("kegg", "KEGG", plot.df$GO_description)
plot.df$GO_description <- gsub("Kegg ", "KEGG: ", plot.df$GO_description)
plot.df$GO_description <- gsub("Gobp ", "GObp: ", plot.df$GO_description)
plot.df$GO_description <- gsub("Reactome ", "Reactome: ", plot.df$GO_description)
plot.df$GO_description <- gsub("p450", "P450", plot.df$GO_description)
plot.df$GO_description <- gsub(" system kks", " system (KKS)", plot.df$GO_description)
plot.df$GO_description <- gsub(" system cas", " system (CAS)", plot.df$GO_description)

# Insert line break here
plot.df$GO_description <- gsub(" system \\(CAS\\)", " system (CAS)\n", plot.df$GO_description)

# Sort by p-value (most significant first)
plot.df <- plot.df %>% arrange(pval)
# Then sort by sign; is not the same as sorting on z-value, which will put lowest z-values last!
plot.df <- plot.df %>% arrange(-sign(NES))
# To keep order in plot!
plot.df$GO_description <- factor(plot.df$GO_description, levels = rev(plot.df$GO_description))
# ### Overrule, take custom order: ###
# plot.df$GO_description <- factor(plot.df$GO_description, levels = rev(c("Mitochondrion (GO:0005739)",
#                                                                   "Lipolysis in adipocytes (KEGG: mmu04923)",
#                                                                   "Structural constituent of muscle (GO:0008307)",
#                                                                   "mRNA splicing, via spliceosome (GO:0000398)")))

hjust <- rep(NA, nrow(plot.df))
hjust[plot.df$NES > 0] <- 1.6
hjust[plot.df$NES < 0] <- -0.7

x_pos_stars <- rep(NA, nrow(plot.df))
x_pos_stars[plot.df$NES > 0] <- 2.7
x_pos_stars[plot.df$NES < 0] <- -2.7

rescale_colors <- function(x){
  x <- c(0,x,1)
  out <- scales::rescale(log(x+0.0001))
  return(out[-c(1,length(out))])
}
labels_scale <- c(0, 0.001, 0.01, 0.05, 0.1, 1)

# To make sure that a zero after the comma is also printed out
numbers <- unlist(rapply(list(plot.df$NES), sprintf, fmt = "%0.1f", how = "replace"))

main <- paste0("")
xlim <- c(-3, 3)

barplot.fgsea.mortality <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = NES, y = GO_description)) + # 1:10
    ggtitle(main) +
    ylab("") + # "Gene set"
    xlab("Normalized Enrichment Score") +
    geom_bar(stat="identity", aes(fill = rescale_colors(padj))) + # fill = colors
    scale_fill_viridis(option = "plasma", name = "q-value", limits = c(0, 1), 
                       values = rescale_colors(labels_scale), 
                       breaks = rescale_colors(labels_scale), 
                       labels = labels_scale,
                       direction = -1,
                       guide = guide_colorbar(reverse = TRUE)) +
    geom_text(aes(label = numbers), hjust = hjust, color="white", size=3, fontface = "bold") + # Not useful, all q-values of 0.01
    geom_text(aes(label = stars), x = x_pos_stars, color="black", size=5) + # , fontface = "bold"
    geom_vline(xintercept = 0) +
    xlim(xlim[1], xlim[2]) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_line(color="black", linewidth = 0.5),
          axis.line.y = element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size=12),
          strip.background = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12, face="bold"), # smaller title, bold
          legend.title = element_text(size = 10, face = "bold")
    )
))

attr(barplot.fgsea.mortality, "data") <- as.data.frame(plot.df)

plot(barplot.fgsea.mortality)

### Results association SCGB1A1 with age and mortality

res.age[res.age$protein == "SCGB1A1", ]
# p < 1e-16
res.mortality[res.mortality$protein == "SCGB1A1", ]
# p = 0.01



