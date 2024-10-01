
### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
coefficients.longitudinal <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

dict <- read.table(paste0(input.GSE169148.dir, "dict.tsv"), sep = "\t", header = TRUE, comment.char = "")

filenames <- paste0(input.GSE169148.dir, "GSE169148_series_matrix.txt.gz")
samples <- geograbi.get.samples(filename = filenames)
metadata <- geograbi::geograbi.extract.characteristics(samples)

dat <- exprs(GEOquery::getGEO(filename = filenames))
any(is.na(dat))
# FALSE
# No filtering or imputation needed here!

tmp.list <- vector(mode = "list", length = length(coefficients.longitudinal$gen2))
names(tmp.list) <- names(coefficients.longitudinal$gen2)

for(k in 1:length(tmp.list)){
  if(!(names(tmp.list)[k] %in% c("Bladder", "Thyroid"))){
    tmp.list[[k]] <- data.frame(organ = names(coefficients.longitudinal$gen2)[k], 
                                protein = names(coefficients.longitudinal$gen2[[k]]),
                                value = coefficients.longitudinal$gen2[[k]])
  }
}

organ.df <- do.call("rbind", tmp.list)
rownames(organ.df) <- NULL
organ.df <- organ.df[organ.df$value != 0, ]

maped_organs = organ.df |>
  inner_join(dict, by = c("protein" = "Assay"))

res <- dat |>
  as_tibble(rownames = "ID") |>
  mutate(ID = as.numeric(ID)) |>
  inner_join(maped_organs, by = "ID") |>
  mutate(across(starts_with("GSM"), ~ .x * value)) |>
  group_by(organ) |>
  summarize(across(starts_with("GSM"), sum))

# Combine the data and metadata
t_res <- res |>
  as.data.frame() |>
  tibble::column_to_rownames(var = "organ") |>
  t() |>
  as_tibble()

data <- cbind(t_res, metadata) |>
  as_tibble()

# Reshape the data from wide to long format
full.df <- melt(
  data,
  id.vars = c("subject status", "treatment", "tissue"),
  variable.name = "organ",
  value.name = "biological_age"
)

full.df$condition <- paste0(full.df$subject, " ", full.df$treatment)

# full.df$condition <- full.df$group
full.df$condition[full.df$condition == "healthy control none"] <- "Healthy control"
full.df$condition[full.df$condition == "Sarcoidosis patient none"] <- "Sarcoidosis"
full.df$condition[full.df$condition == "Sarcoidosis patient tofacitinib"] <- "Sarcoidosis with tofacitinib"

table(full.df$condition[full.df$organ == "Conventional"])
# Healthy control               Sarcoidosis Sarcoidosis with tofacitinib 
# 11                            9                           11 

pvals.full <- rep(NA, length(organ.proteins))
names(pvals.full) <- names(organ.proteins)
adjps.full <- rep(NA, length(organ.proteins))
names(adjps.full) <- names(organ.proteins)
adjps.selected <- rep(NA, length(organ.proteins.selected))
names(adjps.selected) <- names(organ.proteins.selected)

for(k in 1:length(organ.proteins)){
  
  if(!(names(organ.proteins)[k] %in% c("Bladder", "Thyroid"))){
    df <- full.df[full.df$organ == names(organ.proteins)[k],]
    df$condition <- factor(df$condition, levels = c("Healthy control", "Sarcoidosis", "Sarcoidosis with tofacitinib"))
    
    model <- lm(biological_age ~ condition, data = df)
    
    if(!any(is.na(summary(model)$coef[, "Pr(>|t|)"]))){
      summary <- tryCatch(summary(glht(model, 
                                       linfct = mcp(condition = c("(Sarcoidosis + `Sarcoidosis with tofacitinib`)/2 - `Healthy control` = 0")))), 
                          error = function(e){return(NA)})
    }
    pvals.full[k] <- summary$test$pvalues
  }
}

adjps.full <- p.adjust(pvals.full, method = "hommel")
adjps.selected <- p.adjust(pvals.full[names(organ.proteins.selected)], method = "hommel")

adjps.selected
# Conventional        Brain       Artery        Liver    Intestine       Immune       Kidney         Skin         Lung 
# 0.077952537  0.935134158  0.378797765  0.588807934  0.144336660  0.003258034  0.935134158  0.935134158  0.005079618 

### Be careful, we only correct for organ.proteins.selected!
### See line if(names(organ.proteins)[k] %in% organ.proteins.selected)

violin.plots.GSE169148 <- vector(mode = "list", length = length(organ.proteins))
names(violin.plots.GSE169148) <- names(organ.proteins)

for(k in 1:length(violin.plots.GSE169148)){
  
  if(!(names(violin.plots.GSE169148)[k] %in% c("Bladder", "Thyroid"))){
    
    plot.df <- full.df[full.df$organ == names(organ.proteins)[k],]
    
    plot.df$grouping <- plot.df$condition
    plot.df$grouping[plot.df$grouping == "Sarcoidosis with tofacitinib"] <- "Sarcoidosis"
    plot.df$grouping <- factor(plot.df$grouping, levels = c("Sarcoidosis", "Healthy control"))
    
    
    plot.df$condition[plot.df$condition == "Sarcoidosis with tofacitinib"] <- "Sarcoidosis \nwith tofacitinib"
    plot.df$condition <- factor(plot.df$condition, levels = c("Healthy control", "Sarcoidosis", "Sarcoidosis \nwith tofacitinib"))
    
    plot.df$x.numeric <- as.numeric(plot.df$condition)
    plot.df$x.numeric[plot.df$x.numeric == 2] <- 2.0
    plot.df$x.numeric[plot.df$x.numeric == 3] <- 2.6
    
    height1 <- (max(plot.df$biological_age, na.rm = TRUE) - min(plot.df$biological_age, na.rm = TRUE))*0.95+min(plot.df$biological_age, na.rm = TRUE)
    
    significance.df <- data.frame(x = 1.5, 
                                  xend = 1.75, 
                                  y = height1,
                                  adjp = adjps.full[k], 
                                  stars = "n.s.",
                                  size = 4.5,
                                  nudge_y = (max(plot.df$biological_age, na.rm = TRUE) - min(plot.df$biological_age, na.rm = TRUE))/13
    )
    
    if(names(organ.proteins)[k] %in% names(organ.proteins.selected)){
      significance.df$adjp <- adjps.selected[names(organ.proteins)[k]]
    }
    
    significance.df$stars[which(significance.df$adjp < 0.1)] <- "'"
    significance.df$stars[which(significance.df$adjp < 0.05)] <- "*"
    significance.df$stars[which(significance.df$adjp < 0.01)] <- "**"
    significance.df$stars[which(significance.df$adjp < 0.001)] <- "***"
    
    significance.df[(significance.df$stars != "n.s."), "size"] <- 10
    significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- (max(plot.df$biological_age, na.rm = TRUE) - min(plot.df$biological_age, na.rm = TRUE))/30
    
    violin.plots.GSE169148[[k]] <- ggplot_gtable(ggplot_build(
      ggplot(plot.df, aes(x = x.numeric, y = biological_age)) +
        
        geom_quasirandom(cex = 1.5, aes(x = x.numeric, color = condition, alpha = 1)) +
        
        scale_color_manual(values = c("#5087B6", "#E9322D", "#CB9018")) +
        scale_fill_manual(values = c("#5087B6", "#E9322D", "#CB9018")) +
        
        geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA, aes(group = condition)) + # fill = Diet, alpha = 1
        guides(alpha = "none") + # no legend for "alpha"
        scale_x_continuous(breaks = sort(unique(plot.df$x.numeric)), labels = levels(plot.df$condition)) + 
        
        geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
        geom_segment(aes(x = significance.df$x[1]-0.25, y = significance.df$y[1], xend = significance.df$xend[1], yend = significance.df$y[1])) +
        
        ggtitle(paste0(gsub("_", "-", paste(toupper(substr(names(violin.plots.GSE169148)[k], 1, 1)), substr(names(violin.plots.GSE169148)[k], 2, nchar(names(violin.plots.GSE169148)[k])), sep="")), " model")) +
        xlab("") + # No label on the horizontal axis
        ylab("Predicted rel. log(mortality hazard)") +
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
          axis.ticks.x = element_blank(), 
          axis.ticks.y = element_line(color="black", linewidth = 0.5),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
    ))
    
    tmp <- plot.df
    colnames(tmp)[colnames(tmp) == "biological_age"] <- "Predicted rel. log(mortality hazard)"
    attr(violin.plots.GSE169148[[k]], "data") <- tmp
    
  }
}


