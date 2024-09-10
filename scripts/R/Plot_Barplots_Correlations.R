
# From 1_Prepare_Data.R
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

# From 3_Process_GTEx.R
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### "Male", "Female", "Organismal", "Multi-organ" are not being discussed in the manuscript
organ.proteins.1B <- organ.proteins[!(names(organ.proteins) %in% c("Male", "Female", "Organismal", "Multi-organ"))]

df.r.list <- vector(mode = "list", length = 2)
names(df.r.list) <- c("gen1", "gen2")

df.r.list <- lapply(df.r.list, function(x){
  return(
    data.frame(set = rep(c("training dataset", "test dataset"), times = length(organ.proteins.1B)),
               organ = rep(names(organ.proteins.1B), each = 2),
               r = NA,
               r2 = NA)
  )
})

for(g in 1:length(df.r.list)){
  for(i in 1:2){
    for(k in 1:length(organ.proteins.1B)){
      if(!(names(organ.proteins.1B)[k] %in% c("Bladder"))){
        
        if(i == 1){
          df <- data.frame(
            chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% training.eids)],
            biological_age = predicted.ages.1[[g]]$predicted[[names(organ.proteins.1B)[k]]][training.eids],
            sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men"))[(olink_bd_annotation_0$eid %in% training.eids)]
          )
        } else if(i == 2){
          df <- data.frame(
            chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% test.eids)],
            biological_age = predicted.ages.1[[g]]$predicted[[names(organ.proteins.1B)[k]]][test.eids],
            sex = factor(ifelse((olink_bd_annotation_0$p31) == "Female", "Women", "Men"), levels = c("Women", "Men"))[(olink_bd_annotation_0$eid %in% test.eids)]
          )
        }
        
        summary <- summary(lm(biological_age ~ chronological_age, data = df))

        df.r.list[[g]][(df.r.list[[g]]$organ == names(organ.proteins.1B)[k]) & (df.r.list[[g]]$set == c("training dataset", "test dataset")[i]), "r"] <- sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared)
        df.r.list[[g]][(df.r.list[[g]]$organ == names(organ.proteins.1B)[k]) & (df.r.list[[g]]$set == c("training dataset", "test dataset")[i]), "r2"] <- summary$r.squared
        
      }
    }
  } 
}

# saveRDS(df.r.list, file = paste0(rds.dir, "df_r_list.rds"))
df.r.list <- readRDS(file = paste0(rds.dir, "df_r_list.rds"))

### Selection of organs to keep for further analysis, order based on gen2 training dataset ###
ordered.organs <- df.r.list[["gen2"]][which((df.r.list[["gen2"]]$set == "training dataset")),] %>% arrange(-r2) %>% `[`("organ") %>% unlist
ordered.organs <- ordered.organs[ordered.organs != "Bladder"]
ordered.organs <- paste(toupper(substr(ordered.organs, 1, 1)), substr(ordered.organs, 2, nchar(ordered.organs)), sep="")

barplots.age.r <- vector(mode = "list", length = 2)
names(barplots.age.r) <- c("gen1", "gen2")

titles <- c("1st-generation model", "Mortality-based model")

for(g in 1:length(df.r.list)){
  plot.df <- df.r.list[[g]]
  plot.df$color <- plot.df$r > 0.3
  plot.df$organ <- paste(toupper(substr(plot.df$organ, 1, 1)), substr(plot.df$organ, 2, nchar(plot.df$organ)), sep="")
  plot.df <- plot.df[plot.df$organ %in% ordered.organs,]
  plot.df$organ <- factor(plot.df$organ, levels = rev(ordered.organs))
  plot.df$set <- factor(plot.df$set, levels = c("training dataset", "test dataset"))
  
  barplots.age.r[[g]] <- ggplot_gtable(ggplot_build(
    ggplot(data = plot.df, aes(x = r, y = organ)) +
      geom_bar(aes(fill = color), stat='identity', position="dodge") +
      facet_grid(~set, scales = "free_y", space="free") +
      geom_vline(xintercept = 0.3, colour = "darkgrey") +
      ggtitle(paste0(titles[g])) +
      scale_fill_manual(values = c("#f83e31ff", "#006bc4ff")) +
      xlab("Correlation coefficient (r)") +
      ylab("") +
      theme_bw() +
      theme(legend.position = "none", # No legend
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 12),
        strip.text.x = element_text(size= 12), # the size of the facet labels
        axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        axis.ticks.x = element_line(color="black", linewidth = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.5),
        # strip.text = element_text(face = "black", size = 12),
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  ))
  attr(barplots.age.r[[g]], "data") <- plot.df
}

# r consistently larger than 0.3

plot(barplots.age.r[["gen2"]])
plot(barplots.age.r[["gen1"]])

# Save selected organ-specific models

organ.proteins.selected <- organ.proteins[c("Conventional", "Brain", "Artery", "Liver", "Intestine", "Immune", "Kidney", "Skin", "Lung")]

# saveRDS(organ.proteins.selected, file = paste0(rds.dir, "organ_proteins_selected.rds"))
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))


