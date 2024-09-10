
### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

### This comes from Plot_Heatmaps_Medications.R ###
medicine.indices <- readRDS(file = paste0(rds.dir, "medicine_indices.rds"))

### 1. Analyze the effect of metformin after correcting for diabetes ###

df <- data.frame(
  eid = olink_bd_annotation_0$eid,
  metformin = FALSE, # apply(olink_bd_annotation_0[, paste0("p20003_i0_a", 1:47)], 1, function(x){return(any(uk_biobank_codings %in% x))}),
  chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
  biological_age = NA,
  sex = factor(olink_bd_annotation_0$p31, levels = c("Female", "Male"), ordered = FALSE), 
  deprivation = olink_bd_annotation_0$p22189, # Townsend deprivation index
  centre = factor(olink_bd_annotation_0$p54_i0),
  IPAQ_activity = factor(olink_bd_annotation_0$p22032_i0, level = c("low", "moderate", "high"), ordered = FALSE),
  smoking_status = factor(olink_bd_annotation_0$p20116_i0, levels = c("Prefer not to answer", "Never", "Previous", "Current"), ordered = FALSE)
)

df$diabetes <- NA

df$diabetes[(olink_bd_annotation_0$status_type_2_diabetes == 0)] <- "No"
df$diabetes[(olink_bd_annotation_0$time_type_2_diabetes <= 0) & (olink_bd_annotation_0$status_type_2_diabetes == 1)] <- "Past"
df$diabetes[(olink_bd_annotation_0$time_type_2_diabetes > 0) & (olink_bd_annotation_0$status_type_2_diabetes == 1)] <- "Future"

df$diabetes <- factor(df$diabetes, levels = c("No", "Past", "Future"))

# Equal to: medications.lvls[[2]] (see Plot_Heatmaps_Medications.R)
names(medicine.indices[[2]]$full)[13]
# [13] "LIPID MODIFYING AGENTS"

df$lipid_modifying_agents <- 0
df$lipid_modifying_agents[medicine.indices[[2]]$full$`LIPID MODIFYING AGENTS`] <- 1


### 2. Analyze the effect of metformin after correcting for diabetes and simvastatin after correcting for other lipid modifying agents ###

res.meds.corrected <- vector(mode = "list", length = 2)
names(res.meds.corrected) <- c("metformin", "simvastatin")

res.meds.corrected <- lapply(res.meds.corrected, function(x){
  x <- data.frame(organ = names(organ.proteins),
                  Estimate = NA,
                  `Std. Error` = NA,
                  `t value` = NA, 
                  df = NA,
                  `Pr(>|t|)` = NA, 
                  ll = NA, 
                  ul = NA, 
                  r = NA, 
                  check.names = FALSE)
  return(x)
})

for(k in 1:nrow(res.meds.corrected[[1]])){
  # if(names(organ.proteins)[k] != "Bladder"){
  if(names(organ.proteins)[k] %in% names(organ.proteins.selected)){
    
    df$biological_age <- predicted.ages.oof$gen2$predicted[[k]]
    
    for(i in 1:length(res.meds.corrected)){
      
      if(names(res.meds.corrected)[i] == "metformin"){
        
        df$med <- FALSE
        df$med[medicine.indices[[5]][["full"]][["metformin"]]] <- TRUE
        
        model <- lm(biological_age ~ chronological_age * sex + med + diabetes + deprivation + centre + IPAQ_activity + smoking_status, data = df)
        summary.vals <- tryCatch(summary(model)$coef["medTRUE", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
        res.meds.corrected[[i]][k, c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model, parm = "medTRUE"), sign(model$coef["medTRUE"])*sqrt(summary(model)$r.squared))
        res.meds.corrected[[i]][k, c("df")] <- model$df.residual
        
        
      } else if(names(res.meds.corrected)[i] == "simvastatin"){
        
        df$med <- FALSE
        df$med[medicine.indices[[5]][["full"]][["simvastatin"]]] <- TRUE
        
        model <- lm(biological_age ~ chronological_age * sex + med + lipid_modifying_agents + centre + IPAQ_activity + smoking_status, data = df)
        summary.vals <- tryCatch(summary(model)$coef["medTRUE", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
        res.meds.corrected[[i]][k, c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model, parm = "medTRUE"), sign(model$coef["medTRUE"])*sqrt(summary(model)$r.squared))
        res.meds.corrected[[i]][k, c("df")] <- model$df.residual
        
      }
      
    }
  }
}

for(i in 1:length(res.meds.corrected)){
  res.meds.corrected[[i]]$med <- names(res.meds.corrected)[i]
}

# saveRDS(res.meds.corrected, file = paste0(rds.dir, "res_meds_corrected.rds"))
res.meds.corrected <- readRDS(file = paste0(rds.dir, "res_meds_corrected.rds"))

### 2. Make the heatmap ###

plot.df <- do.call("rbind", res.meds.corrected)
plot.df <- plot.df[plot.df$organ %in% names(organ.proteins.selected),]

# Add q-value per group:
plot.df <- plot.df %>% group_by(med) %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

lims <- c(-1.5, 1.5)

range(plot.df$Estimate, na.rm = TRUE)
# -0.05419004  0.48101512

c(quantile(plot.df$Estimate, 0.05, na.rm = TRUE), quantile(plot.df$Estimate, 0.95, na.rm = TRUE))
# -0.03799033  0.34103266

# Only here change to factor!

plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))

plot.df$med <- paste(toupper(substr(plot.df$med, 1, 1)), substr(plot.df$med, 2, nchar(plot.df$med)), sep="")
plot.df$med <- factor(plot.df$med, levels = rev(c("Simvastatin", "Metformin")))

Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

heatmap.metformin.simvastatin <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = med, y = organ, label = stars)) + 
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

attr(heatmap.metformin.simvastatin, "data") <- as.data.frame(plot.df)

plot(heatmap.metformin.simvastatin)


