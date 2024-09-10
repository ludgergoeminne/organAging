
### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

total_weight_by_food_group_yesterday <- read.table(paste0(input.coding.dir, "total_weight_by_food_group_yesterday.tsv"), quote = "", comment.char = "", sep = "\t", header = TRUE)
total_weight_by_food_group_yesterday$Field.ID <- as.character(total_weight_by_food_group_yesterday$Field.ID)
total_weight_by_food_group_yesterday$Description <- as.character(total_weight_by_food_group_yesterday$Description)

# Category 100020 say if diet yesterday was typical!!!
# Data-Field 100026:	Daily dietary data not credible!!!

res.list.foods.gen2 <- vector(mode = "list", length = length(organ.proteins))
names(res.list.foods.gen2) <- names(organ.proteins)

res.list.foods.gen2 <- lapply(res.list.foods.gen2, function(x){
  return(data.frame(food = total_weight_by_food_group_yesterday$Field.ID,
                    meaning = total_weight_by_food_group_yesterday$Description,
                    Estimate = NA,
                    `Std. Error` = NA,
                    `t value` = NA, 
                    df = NA,
                    `Pr(>|t|)` = NA, 
                    ll = NA, 
                    ul = NA, 
                    r = NA, 
                    check.names = FALSE))
})

tmp_bd_annotation_0 <- olink_bd_annotation_0

for(k in 1:length(res.list.foods.gen2)){
  
  if(names(organ.proteins)[k] != "Bladder"){
    
    df <- data.frame(
      eid = olink_bd_annotation_0$eid,
      food = NA, # olink_bd_annotation_0[, paste0("p", total_weight_by_food_group_yesterday$Field.ID[i], "_i0")],
      chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
      biological_age = predicted.ages.oof$gen2$predicted[[k]],
      sex = factor(olink_bd_annotation_0$p31, levels = c("Female", "Male"), ordered = FALSE), 
      deprivation = olink_bd_annotation_0$p22189, # Townsend deprivation index
      centre = factor(olink_bd_annotation_0$p54_i0),
      IPAQ_activity = factor(olink_bd_annotation_0$p22032_i0, level = c("low", "moderate", "high"), ordered = FALSE),
      smoking_status = factor(olink_bd_annotation_0$p20116_i0, levels = c("Prefer not to answer", "Never", "Previous", "Current"), ordered = FALSE)
    )
    
    for(i in 1:nrow(total_weight_by_food_group_yesterday)){
      
      ### Exclude people with non-credible or non-typical diets ###
      tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i0 == "No" | tmp_bd_annotation_0$p100026_i0 == "No"), paste0("p", total_weight_by_food_group_yesterday$Field.ID[i], "_i0")] <- NA
      tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i1 == "No" | tmp_bd_annotation_0$p100026_i1 == "No"), paste0("p", total_weight_by_food_group_yesterday$Field.ID[i], "_i1")] <- NA
      tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i2 == "No" | tmp_bd_annotation_0$p100026_i2 == "No"), paste0("p", total_weight_by_food_group_yesterday$Field.ID[i], "_i2")] <- NA
      tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i3 == "No" | tmp_bd_annotation_0$p100026_i3 == "No"), paste0("p", total_weight_by_food_group_yesterday$Field.ID[i], "_i3")] <- NA
      tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i4 == "No" | tmp_bd_annotation_0$p100026_i4 == "No"), paste0("p", total_weight_by_food_group_yesterday$Field.ID[i], "_i4")] <- NA
      
      # All food questionnaires were sent out after the first and before the third visit!
      df$food <- rowMeans(tmp_bd_annotation_0[, paste0("p", total_weight_by_food_group_yesterday$Field.ID[i], "_i", 0:4)], na.rm = TRUE)
      
      # !!! Keep only the ones for which food info is available! Otherwise, all foods are younger !!!
      # df <- df[unlist(lapply(foods.level.1, length) != 0),]

          model <- lm(biological_age ~ chronological_age * sex + food + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          # summary(model)
          summary.vals <- tryCatch(summary(model)$coef["food", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
          res.list.foods.gen2[[k]][res.list.foods.gen2[[k]]$food == total_weight_by_food_group_yesterday$Field.ID[i], c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model, parm = "food"), sign(model$coef["food"])*sqrt(summary(model)$r.squared))
          res.list.foods.gen2[[k]][res.list.foods.gen2[[k]]$food == total_weight_by_food_group_yesterday$Field.ID[i], "df"] <- model$df.residual
    }
    
    # Add q-value per group:
    res.list.foods.gen2[[k]] <- res.list.foods.gen2[[k]] %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))
    
  }
  
}
rm(tmp_bd_annotation_0)

for(k in 1:length(res.list.foods.gen2)){
  res.list.foods.gen2[[k]]$organ <- names(res.list.foods.gen2)[k]
}

res.list.foods.gen2[["Bladder"]] <- NULL
res.df.foods.gen2 <- do.call("rbind", res.list.foods.gen2)

# saveRDS(res.df.foods.gen2, file = paste0(rds.dir, "res_df_foods_gen2.rds"))
res.df.foods.gen2 <- readRDS(paste0(rds.dir, "res_df_foods_gen2.rds"))

### Heatmap top foods ###

res.df.foods.gen2 <- res.df.foods.gen2[res.df.foods.gen2$organ %in% names(organ.proteins.selected), ]

tmp <- res.df.foods.gen2[(res.df.foods.gen2$organ == "Conventional"),] %>% arrange(`Pr(>|t|)`)
# tmp <- res.df.foods.gen2 %>% arrange(`Pr(>|t|)`)
# tmp <- tmp[!duplicated(tmp$food),]

top.foods <- c(tmp[tmp$Estimate < 0,] %>% head(10) %>% `[`("meaning") %>% unlist, 
               tmp[tmp$Estimate > 0,] %>% head(10) %>% `[`("meaning") %>% unlist %>% rev)

# saveRDS(top.foods, file = "/PHShome/lg002/UKBB/data/top_foods.rds")
top.foods <- readRDS(file = "/PHShome/lg002/UKBB/data/top_foods.rds")

plot.df <- res.df.foods.gen2

# Add q-value per group:
plot.df <- plot.df %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))

# This should come after calculating the q-values!
plot.df <- plot.df[(plot.df$meaning %in% top.foods),]

# Arrange data frame by top foods (to turn into a factor later)
plot.df <- plot.df %>% arrange(match(meaning, top.foods))
all(unique(plot.df$meaning) == top.foods)
# TRUE
plot.df$meaning <- paste(toupper(substr(plot.df$meaning, 1, 1)), substr(plot.df$meaning, 2, nchar(plot.df$meaning)), sep="")
plot.df$meaning <- gsub(" and ", " & ", plot.df$meaning)

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

lims <- c(-0.007, 0.007)

range(plot.df$Estimate, na.rm = TRUE)
# -0.007176112  0.012899285

c(quantile(plot.df$Estimate, 0.05, na.rm = TRUE), quantile(plot.df$Estimate, 0.95, na.rm = TRUE))
# -0.001341260  0.001467783 

# Remove Bladder
plot.df <- plot.df[plot.df$organ != "Bladder",]

# Only here change to factor!
plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))
plot.df$meaning = factor(plot.df$meaning, levels = rev(unique(plot.df$meaning))) # rev(meanings)

Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

heatmap.top.foods.gen2 <- ggplot_gtable(ggplot_build(
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

attr(heatmap.top.foods.gen2, "data") <- plot.df

plot(heatmap.top.foods.gen2)



