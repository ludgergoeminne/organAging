
### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from 4_calculate_predicted_resid_ages.R ###
coefficients.longitudinal <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal.rds"))
predicted.ages.longitudinal.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_longitudinal_1.rds"))
predicted.ages.longitudinal.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_longitudinal_oof.rds"))
longitudinal.training.eids <- readRDS(file = paste0(rds.dir, "longitudinal_training_eids.rds"))
longitudinal.test.eids <- readRDS(file = paste0(rds.dir, "longitudinal_test_eids.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### Import the Parkinson's Disease / MSA dataset (Dammer et al., 2022) ###
### This comes from 6_Prepare_Parkinson_MSA_Data.R ###
PD_annotation_list <- readRDS(file = paste0(rds.dir, "PD_annotation_list.rds"))
predicted.ages.PD <- readRDS(file = paste0(rds.dir, "predicted_ages_PD.rds"))


### Import the COVID dataset (Filbin et al., 2021) ###
### This comes from Plots_Violins_COVID.R ###
df.COVID.list <- readRDS(file = paste0(rds.dir, "df_COVID_list.rds"))

### Put the data together in a data frame ###

sets <- c("Training set", "Test set", "Out of fold", "Third visit", "Fourth visit", "Dammer et al. (2022) - first visit", "Dammer et al. (2022) - last visit", "Filbin et al. (2021)")

summary.list <- vector(mode = "list", length = 2)
names(summary.list) <- c("1st-generation models", "Mortality-based models")

summary.list <- lapply(summary.list, function(x){
  x <- vector(mode = "list", length = 5)
  names(x) <- c("r", "r²", "MAE", "MSE", "pval")
  
  x <- lapply(x, function(x){
    y <- data.frame(
      organ = rep(names(organ.proteins), each = length(sets)),
      set = rep(sets, times = length(organ.proteins)),
      value = NA
    )
    return(y)
  })
  
  return(x)
})

for(g in 1:length(summary.list)){
  
  for(k in 1:length(organ.proteins)){
    if(!(names(organ.proteins)[k] %in% c("Bladder", "Thyroid"))){
      
      for(i in 1:length(sets)){
        
        if(sets[i] == "Training set"){
          
          df <- data.frame(
            chronological_age = olink_bd_annotation_list[["instance_0"]]$age_first_visit[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.training.eids)],
            biological_age = predicted.ages.longitudinal.1[["instance_0"]][[g]]$predicted[[k]][longitudinal.training.eids])
          
        } else if(sets[i] == "Test set"){
          
          df <- data.frame(
            chronological_age = olink_bd_annotation_list[["instance_0"]]$age_first_visit[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.test.eids)],
            biological_age = predicted.ages.longitudinal.1[["instance_0"]][[g]]$predicted[[k]][longitudinal.test.eids])
          
        } else if(sets[i] == "Out of fold"){
          
          df <- data.frame(
            chronological_age = olink_bd_annotation_list[["instance_0"]]$age_first_visit,
            biological_age = predicted.ages.longitudinal.oof[[g]]$predicted[[k]])
          
        } else if(sets[i] == "Third visit"){
          
          df <- data.frame(
            chronological_age = olink_bd_annotation_list[["instance_2"]]$age_third_visit,
            biological_age = predicted.ages.longitudinal.1[["instance_2"]][[g]]$predicted[[k]])
          
        } else if(sets[i] == "Fourth visit"){
          
          df <- data.frame(
            chronological_age = olink_bd_annotation_list[["instance_3"]]$age_fourth_visit,
            biological_age = predicted.ages.longitudinal.1[["instance_3"]][[g]]$predicted[[k]])
          
        } else if(sets[i] == "Dammer et al. (2022) - first visit"){
          
          df <- cbind(data.frame(
            biological_age = predicted.ages.PD[["plasma"]][[g]][["predicted"]][[k]]), 
            PD_annotation_list[["plasma"]]
          )
          
          ### Only select the last visit here ###
          df <- df[df$visit_month == 0, ]
          df <- df[!duplicated(df$participant_id),] # see prepare_CSF_models.R
          # Remove empty age values
          df <- df[!is.na(df$age_current),]
          # Rename to align with everything else
          df$chronological_age <- df$age_current
          
          df <- df[, c("biological_age", "chronological_age")]
          
        } else if(sets[i] == "Dammer et al. (2022) - last visit"){
          
          df <- cbind(data.frame(
            biological_age = predicted.ages.PD[["plasma"]][[g]][["predicted"]][[k]]), 
            PD_annotation_list[["plasma"]]
          )
          
          ### Only select the last visit here ###
          df <- df %>% arrange(-visit_month)
          df <- df[!duplicated(df$participant_id),] # see prepare_CSF_models.R
          # Remove empty age values
          df <- df[!is.na(df$age_current),]
          # Rename to align with everything else
          df$chronological_age <- df$age_current
          
          df <- df[, c("biological_age", "chronological_age")]
          
        } else if(sets[i] == "Filbin et al. (2021)"){
          
          df.0 <- df.COVID.list[[k]][df.COVID.list[[k]]$Day == 0,]
          df.0$chronological_age <- df.0[, "age.numeric"]
          df.0$biological_age <- df.0[, "predicted.ages"]
          
          df <- df.0[, c("biological_age", "chronological_age")]
        }
        
        model <- lm(biological_age ~ chronological_age, data = df)
        summary <- summary(model)
       
          for(j in 1:length(summary.list[[g]])){
            
            if(names(summary.list[[g]])[j] == "r"){
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared)
            } else if(names(summary.list[[g]])[j] == "r²"){
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- summary$r.squared
            } else if(names(summary.list[[g]])[j] == "MAE"){
              resid <- resid(model)
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- mean(abs(resid))
            } else if(names(summary.list[[g]])[j] == "MSE"){
              resid <- resid(model)
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- mean(resid^2)
            } else if(names(summary.list[[g]])[j] == "pval"){
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- summary$coefficients["chronological_age", "Pr(>|t|)"]
            }
            
          }
      
    }
  }
  
  }
}

summary.df.list.longitudinal <- lapply(summary.list, function(x){
  x <- do.call("cbind", x)
  x <- x[, c(1, 2, seq(3, ncol(x), by = 3))]
  colnames(x)[c(1, 2)] <- gsub("^.+?\\.", "", colnames(x)[c(1, 2)])
  colnames(x)[-c(1, 2)] <- gsub("\\..+$", "", colnames(x)[-c(1, 2)])
  return(x)
})

# Needed for Export_Supp_Tables.R
# saveRDS(summary.df.list.longitudinal, file = paste0(rds.dir, "summary_df_list_longitudinal.rds"))
summary.df.list.longitudinal <- readRDS(file = paste0(rds.dir, "summary_df_list_longitudinal.rds"))

### Plot the correlation heatmaps for the non-longitudinal models ###

heatmaps.r.longitudinal <- vector(mode = "list", length = length(summary.df.list.longitudinal))
names(heatmaps.r.longitudinal) <- names(summary.df.list.longitudinal)

for(g in 1:length(heatmaps.r.longitudinal)){
  
  plot.df <- summary.df.list.longitudinal[[g]]
  plot.df <- plot.df[plot.df$organ %in% names(organ.proteins.selected),]
  
  # Add q-value per group:
  plot.df <- plot.df %>% mutate(qval = p.adjust(`pval`, method = "BH"))
  
  plot.df$stars <- ""
  plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
  plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
  plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
  plot.df$stars[which(plot.df$qval < 0.001)] <- "***"
  
  lims <- c(-1, 1)
  
  # Only here change to factor!
  plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))
  plot.df$set = factor(plot.df$set, levels = rev(unique(plot.df$set)))
  
  Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  
  heatmaps.r.longitudinal[[g]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = organ, y = set, label = stars)) + 
      geom_tile(aes(fill = r)) + # size = minuslog10p , color = p.adjust
      # scale_shape_manual(values=c(22, 24, 25)) +
      # facet_grid(pathway~diet, scales = "free_y", space="free") +
      scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
      guides(color = "none") + # no legend for color
      labs(fill="r") + 
      geom_text(size = 4.5, nudge_y = -0.13, na.rm = TRUE) +
      xlab("") +
      ylab("") +
      ggtitle("") +
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
  attr(heatmaps.r.longitudinal[[g]], "data") <- plot.df
}

plot(heatmaps.r.longitudinal[[1]])
plot(heatmaps.r.longitudinal[[2]])

### Plot the MAE heatmaps for the non-longitudinal models ###

heatmaps.MAE.longitudinal <- vector(mode = "list", length = length(summary.df.list.longitudinal))
names(heatmaps.MAE.longitudinal) <- names(summary.df.list.longitudinal)

lims.max <- c(3.5, 0.6)

for(g in 1:length(heatmaps.MAE.longitudinal)){
  
  plot.df <- summary.df.list.longitudinal[[g]]
  plot.df <- plot.df[plot.df$organ %in% names(organ.proteins.selected),]
  
  lims <- c(0, lims.max[g])
  
  # Only here change to factor!
  plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))
  plot.df$set = factor(plot.df$set, levels = rev(unique(plot.df$set)))
  
  Heatmap_palette <- c("white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  
  heatmaps.MAE.longitudinal[[g]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = organ, y = set)) + 
      geom_tile(aes(fill = MAE)) + # size = minuslog10p , color = p.adjust
      scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
      guides(color = "none") + # no legend for color
      labs(fill="MAE") + 
      xlab("") +
      ylab("") +
      ggtitle("") +
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
  attr(heatmaps.MAE.longitudinal[[g]], "data") <- plot.df
}

plot(heatmaps.MAE.longitudinal[[1]])
plot(heatmaps.MAE.longitudinal[[2]])


