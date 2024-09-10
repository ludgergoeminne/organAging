
### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

# Participants were asked to select a job group by navigating down a three-level tree. 
# This field records the code they selected at the final level of the tree. 

# Let's go one level up, hence pick only selectable, then take parent ID
# Then check what parent ID means in the non-selectable

coding497 <- read.table(paste0(input.coding.dir, "coding497.tsv"), quote = "", comment.char = "", sep = "\t", header = TRUE)
coding497$coding <- as.character(coding497$coding)
coding497$meaning <- as.character(coding497$meaning)
coding497$node_id <- as.character(coding497$node_id)
coding497$parent_id <- as.character(coding497$parent_id)
coding497$selectable <- as.character(coding497$selectable)

lowest.level <- coding497[coding497$selectable == "Y",]
level.1 <- coding497[coding497$node_id %in% unique(lowest.level$parent_id),]

jobs.level.1 <- vector(mode = "list", length = nrow(olink_bd_annotation_0))

for(i in 1:length(jobs.level.1)){
  tmp <- as.character(na.omit(unlist(olink_bd_annotation_0[i, paste0("p22601_a", 0:39)])))
  jobs.level.1[[i]] <- lowest.level[lowest.level$coding %in% tmp, "parent_id"] # Should already be unique
}

all.jobs.lvl1 <- vector(mode = "list", length = length(level.1$node_id))
names(all.jobs.lvl1) <- level.1$node_id

for(i in 1:length(all.jobs.lvl1)){
  all.jobs.lvl1[[i]] <- unlist(lapply(jobs.level.1, function(x){
    level.1$node_id[i] %in% x
  }))
}

table(all.jobs.lvl1[[1]])

res.list.jobs.gen2 <- vector(mode = "list", length = length(organ.proteins))
names(res.list.jobs.gen2) <- names(organ.proteins)

res.list.jobs.gen2 <- lapply(res.list.jobs.gen2, function(x){
  return(data.frame(job = level.1$node_id,
                    meaning = level.1$meaning,
                    Estimate = NA,
                    `Std. Error` = NA,
                    `t value` = NA, 
                    `Pr(>|t|)` = NA, 
                    ll = NA, 
                    ul = NA, 
                    r = NA, 
                    check.names = FALSE))
})

# We only retain the job codings for which there actually are participants
indices.jobs <- which(unlist(lapply(jobs.level.1, length) != 0))

for(k in 1:length(res.list.jobs.gen2)){
  
  if(names(organ.proteins)[k] != "Bladder"){
    
    df <- data.frame(
      eid = olink_bd_annotation_0$eid,
      job = NA,
      chronological_age = olink_bd_annotation_0$age_first_visit,
      biological_age = predicted.ages.oof$gen2$predicted[[k]],
      sex = factor(olink_bd_annotation_0$p31, levels = c("Female", "Male"), ordered = FALSE), 
      deprivation = olink_bd_annotation_0$p22189, # Townsend deprivation index
      centre = factor(olink_bd_annotation_0$p54_i0),
      IPAQ_activity = factor(olink_bd_annotation_0$p22032_i0, level = c("low", "moderate", "high"), ordered = FALSE),
      smoking_status = factor(olink_bd_annotation_0$p20116_i0, levels = c("Prefer not to answer", "Never", "Previous", "Current"), ordered = FALSE)
    )
    
    df <- df[indices.jobs,]
    
    for(i in 1:length(all.jobs.lvl1)){
      
      df$job <- all.jobs.lvl1[[i]][indices.jobs]
      
          model <- lm(biological_age ~ chronological_age * sex + job + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          # summary(model)
          summary.vals <- tryCatch(summary(model)$coef["jobTRUE", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
          res.list.jobs.gen2[[k]][(res.list.jobs.gen2[[k]]$job == names(all.jobs.lvl1)[i]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model, parm = "jobTRUE"), sign(model$coef["jobTRUE"])*sqrt(summary(model)$r.squared))
 
    }
    
    # Add q-value per group:
    res.list.jobs.gen2[[k]] <- res.list.jobs.gen2[[k]] %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))
    
  }
  
}

for(k in 1:length(res.list.jobs.gen2)){
  res.list.jobs.gen2[[k]]$organ <- names(res.list.jobs.gen2)[k]
}

res.list.jobs.gen2[["Bladder"]] <- NULL
res.df.jobs.gen2 <- do.call("rbind", res.list.jobs.gen2)

# saveRDS(res.df.jobs.gen2, file = paste0(rds.dir, "res_df_jobs_gen2.rds"))
res.df.jobs.gen2 <- readRDS(file = paste0(rds.dir, "res_df_jobs_gen2.rds"))


### Heatmap top jobs ###

res.df.jobs.gen2 <- res.df.jobs.gen2[res.df.jobs.gen2$organ %in% names(organ.proteins.selected), ]

tmp <- res.df.jobs.gen2[(res.df.jobs.gen2$organ == "Conventional"),] %>% arrange(`Pr(>|t|)`)

top.jobs <- c(tmp[tmp$Estimate < 0,] %>% head(10) %>% `[`("meaning") %>% unlist, 
              tmp[tmp$Estimate > 0,] %>% head(10) %>% `[`("meaning") %>% unlist)

plot.df <- res.df.jobs.gen2[(res.df.jobs.gen2$meaning %in% top.jobs),]

# Add q-value per group:
plot.df <- plot.df %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))

plot.df$meaning <- paste(toupper(substr(plot.df$meaning, 1, 1)), substr(plot.df$meaning, 2, nchar(plot.df$meaning)), sep="")
plot.df$meaning <- gsub("human resources", "HR", plot.df$meaning)
plot.df$meaning <- gsub(" and ", " & ", plot.df$meaning)
plot.df$meaning <- gsub(" such as ", ", e.g. ", plot.df$meaning)
plot.df$meaning <- gsub("Manager or senior", "Manager/senior", plot.df$meaning)
plot.df$meaning <- gsub("; property/housing manager", ", property/housing", plot.df$meaning)
plot.df$meaning <- gsub(" in sales ", ": sales ", plot.df$meaning)
plot.df$meaning <- gsub(" in the production of ", " in producing ", plot.df$meaning)

plot.df$organ <- gsub("_", " ", plot.df$organ)
plot.df$organ <- paste(toupper(substr(plot.df$organ, 1, 1)), substr(plot.df$organ, 2, nchar(plot.df$organ)), sep="")

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

lims <- c(-0.6, 0.6)

range(plot.df$Estimate, na.rm = TRUE)
# -0.5964480  0.5859368

c(quantile(plot.df$Estimate, 0.05, na.rm = TRUE), quantile(plot.df$Estimate, 0.95, na.rm = TRUE))
# -0.1741207  0.2375542

# Remove Bladder
plot.df <- plot.df[plot.df$organ != "Bladder",]

# Only here change to factor!
plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))
meanings <- plot.df[(plot.df$organ == "Conventional"),] %>% arrange(`t value`) %>% `[`("meaning") %>% unlist
plot.df$meaning = factor(plot.df$meaning, levels = rev(meanings))

Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

heatmap.top.jobs.gen2 <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = organ, y = meaning, label = stars)) + 
    geom_tile(aes(fill = Estimate)) + # size = minuslog10p , color = p.adjust
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    guides(color = "none") + # no legend for color
    labs(fill="Effect size") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
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

attr(heatmap.top.jobs.gen2, "data") <- plot.df

plot(heatmap.top.jobs.gen2)


