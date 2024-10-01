
# From 1_Prepare_Data.R
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

# From 3_Process_GTEx.R
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

# From 4_calculate_predicted_resid_ages.R
coefficients <- readRDS(paste0(rds.dir, "coefficients.rds"))
predicted.ages.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_1.rds"))
training.eids <- readRDS(file = paste0(rds.dir, "training_eids.rds"))
test.eids <- readRDS(file = paste0(rds.dir, "test_eids.rds"))

# To export the summarized data to Excel:
organs.save.xlsx <- c("Conventional", "Artery", "Brain")
# In case we don't need to export the data to Excel (faster):
# organs.save.xlsx <- NULL

### 1. Plot the predictions in training and test datasets for the 1st-generation models ###

prediction.plots.gen1 <- vector(mode = "list", length = 2)
names(prediction.plots.gen1) <- c("training dataset", "test dataset")

for(i in 1:length(prediction.plots.gen1)){
  prediction.plots.gen1[[i]] <- vector(mode = "list", length = length(organ.proteins))
  names(prediction.plots.gen1[[i]]) <- names(organ.proteins)
}

names <- paste(toupper(substr(names(organ.proteins), 1, 1)), substr(names(organ.proteins), 2, nchar(names(organ.proteins))), sep="")
names <- gsub("_", "-", names)

for(i in 1:length(prediction.plots.gen1)){
  for(k in 1:length(prediction.plots.gen1[[i]])){
    if(!(names(organ.proteins)[k] %in% c("Bladder"))){
      
      if(names(prediction.plots.gen1)[i] == "training dataset"){
        plot.df <- data.frame(
          chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% training.eids)],
          biological_age = predicted.ages.1$gen1$predicted[[k]][training.eids],
          sex = olink_bd_annotation_0$p31[(olink_bd_annotation_0$eid %in% training.eids)])
      } else if(names(prediction.plots.gen1)[i] == "test dataset"){
        plot.df <- data.frame(
          chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% test.eids)],
          biological_age = predicted.ages.1$gen1$predicted[[k]][test.eids],
          sex = olink_bd_annotation_0$p31[(olink_bd_annotation_0$eid %in% test.eids)])
      }
      
      summary <- summary(lm(biological_age ~ chronological_age, data = plot.df))
      resid <- resid(lm(biological_age ~ chronological_age, data = plot.df))
      metrics <- paste0("r = ", sprintf("%.2f", round(sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared), 2)), ", r² = ", sprintf("%.2f", round(summary$r.squared, 2)), ", MAE = ", sprintf("%.2f", round(mean(abs(resid)), 2)))
      
      prediction.plots.gen1[[i]][[k]] <- ggplot_gtable(ggplot_build(
        ggplot(data = plot.df, aes(x = chronological_age, y = biological_age)) +
          geom_point(size = 1, aes(color = sex, shape = sex, alpha = 0.05)) + 
          scale_color_manual(values=c("magenta", "#4169E1")) +
          geom_smooth(method = MASS::rlm , color="turquoise", fill="turquoise", se=TRUE) + 
          ggtitle(paste0(names[k], " 1st-generation model: \n", names(prediction.plots.gen1)[i], ": ", sum(coefficients[["gen1"]][[k]][-1] != 0), " / ", length(organ.proteins[[k]]), " proteins", "\n", metrics)) +
          xlab("Chronological age") +
          ylab("Predicted age") +
          theme_bw() +
          theme(
            legend.position = "none", # No legend
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
            panel.grid.minor = element_blank())
      ))
      
      if(!is.null(organs.save.xlsx)){
        if(names(prediction.plots.gen1[[i]])[k] %in% organs.save.xlsx){
          
          tmp.male <- kMeansEqual(kdat = plot.df[plot.df$sex == "Male", ], 
                                  k = floor(nrow(plot.df[plot.df$sex == "Male", ])/10),
                                  columns = c("chronological_age", "biological_age"),
                                  add.centers = FALSE
          )
          tmp.female <- kMeansEqual(kdat = plot.df[plot.df$sex == "Female", ], 
                                    k = floor(nrow(plot.df[plot.df$sex == "Female", ])/10),
                                    columns = c("chronological_age", "biological_age"),
                                    add.centers = FALSE
          )
          
          tmp.male <- tmp.male %>% group_by(assigned, sex) %>% summarize(
            `N participants` = n(),
            `Avg. chronological age` = mean(chronological_age, na.rm = TRUE), 
            `Avg. biological age` = mean(biological_age, na.rm = TRUE)
          ) 
          
          tmp.female <- tmp.female %>% group_by(assigned, sex) %>% summarize(
            `N participants` = n(),
            `Avg. chronological age` = mean(chronological_age, na.rm = TRUE), 
            `Avg. biological age` = mean(biological_age, na.rm = TRUE)
          ) 
          
          tmp <- rbind(tmp.male, tmp.female) %>% arrange(`Avg. chronological age`)
          tmp <- tmp[, c("N participants", "sex", "Avg. chronological age", "Avg. biological age")]
          colnames(tmp)[2] <- "Sex"
          
          attr(prediction.plots.gen1[[i]][[k]], "data") <- as.data.frame(tmp)
          rm(tmp.male, tmp.female, tmp)
        }
      }
      
    }
  }
}

plot(prediction.plots.gen1[["training dataset"]]$Conventional)
plot(prediction.plots.gen1[["test dataset"]]$Conventional)

### 2. Plot the predictions in training and test datasets for the mortality-based models ###

prediction.plots.gen2 <- vector(mode = "list", length = 2)
names(prediction.plots.gen2) <- c("training dataset", "test dataset")

for(i in 1:length(prediction.plots.gen2)){
  prediction.plots.gen2[[i]] <- vector(mode = "list", length = length(organ.proteins))
  names(prediction.plots.gen2[[i]]) <- names(organ.proteins)
}

for(i in 1:length(prediction.plots.gen2)){
  for(k in 1:length(prediction.plots.gen2[[i]])){
    if(!(names(organ.proteins)[k] %in% c("Bladder"))){ # , "Thyroid"
      
      if(names(prediction.plots.gen2)[i] == "training dataset"){
        plot.df <- data.frame(
          chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% training.eids)],
          biological_age = predicted.ages.1$gen2$predicted[[k]][training.eids],
          sex = olink_bd_annotation_0$p31[(olink_bd_annotation_0$eid %in% training.eids)])
      } else if(names(prediction.plots.gen2)[i] == "test dataset"){
        plot.df <- data.frame(
          chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% test.eids)],
          biological_age = predicted.ages.1$gen2$predicted[[k]][test.eids],
          sex = olink_bd_annotation_0$p31[(olink_bd_annotation_0$eid %in% test.eids)])
      }
      
      summary <- summary(lm(biological_age ~ chronological_age, data = plot.df))
      resid <- resid(lm(biological_age ~ chronological_age, data = plot.df))
      metrics <- paste0("r = ", sprintf("%.2f", round(sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared), 2)), ", r² = ", sprintf("%.2f", round(summary$r.squared, 2)), ", MAE = ", sprintf("%.2f", round(mean(abs(resid)), 2)))
      
      prediction.plots.gen2[[i]][[k]] <- ggplot_gtable(ggplot_build(
        ggplot(data = plot.df, aes(x = chronological_age, y = biological_age)) +
          # geom_boxplot(aes(color = diet)) + # fill = diet
          geom_point(size = 1, aes(color = sex, shape = sex, alpha = 0.05)) + # aes(color = Gender, shape = Sequence)
          scale_color_manual(values=c("magenta", "#4169E1")) +
          geom_smooth(method = MASS::rlm , color="turquoise", fill="turquoise", se=TRUE) + 
          ggtitle(paste0(names(organ.proteins)[k], " mortality-based model: \n", names(prediction.plots.gen2)[i], ": ", sum(coefficients[["gen2"]][[k]] != 0), " / ", length(organ.proteins[[k]]), " proteins", "\n", metrics)) +
          xlab("Chronological age") +
          ylab("Predicted rel. log(mortality hazard)") +
          theme_bw() +
          theme(
            legend.position = "none", # No legend
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
      
      if(!is.null(organs.save.xlsx)){
        if(names(prediction.plots.gen2[[i]])[k] %in% organs.save.xlsx){
          
          tmp.male <- kMeansEqual(kdat = plot.df[plot.df$sex == "Male", ], 
                                           k = floor(nrow(plot.df[plot.df$sex == "Male", ])/10),
                                           columns = c("chronological_age", "biological_age"),
                                           add.centers = FALSE
                               )
          tmp.female <- kMeansEqual(kdat = plot.df[plot.df$sex == "Female", ], 
                                           k = floor(nrow(plot.df[plot.df$sex == "Female", ])/10),
                                           columns = c("chronological_age", "biological_age"),
                                           add.centers = FALSE
                               )
          
          tmp.male <- tmp.male %>% group_by(assigned, sex) %>% summarize(
            `N participants` = n(),
            `Avg. chronological age` = mean(chronological_age, na.rm = TRUE), 
            `Avg. biological age` = mean(biological_age, na.rm = TRUE)
            ) 
          
          tmp.female <- tmp.female %>% group_by(assigned, sex) %>% summarize(
            `N participants` = n(),
            `Avg. chronological age` = mean(chronological_age, na.rm = TRUE), 
            `Avg. biological age` = mean(biological_age, na.rm = TRUE)
          ) 
          
          tmp <- rbind(tmp.male, tmp.female) %>% arrange(`Avg. chronological age`)
          tmp <- tmp[, c("N participants", "sex", "Avg. chronological age", "Avg. biological age")]
          colnames(tmp)[2] <- "Sex"
          
          attr(prediction.plots.gen2[[i]][[k]], "data") <- as.data.frame(tmp)
          rm(tmp.male, tmp.female, tmp)
        }
      }
      
    }
  }
}

plot(prediction.plots.gen2[["training dataset"]]$Conventional)
plot(prediction.plots.gen2[["test dataset"]]$Conventional)



