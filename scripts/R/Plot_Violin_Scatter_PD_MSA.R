
### This comes from 6_Prepare_Parkinson_MSA_Data.R ###
PD_annotation_list <- readRDS(file = paste0(rds.dir, "PD_annotation_list.rds"))
PD_olink_data_list_imputed <- readRDS(file = paste0(rds.dir, "cleanDat_Olink_imputed_list.rds"))
coefficients.longitudinal.Parkinson <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal_Parkinson.rds"))
predicted.ages.PD <- readRDS(file = paste0(rds.dir, "predicted_ages_PD.rds"))

# Only needed for export to Excel
xlsx.quantiles <- c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99, 0.999, 1)

### 1. Plot correlations with chronological age ###

prediction.plots.Parkinson <- vector(mode = "list", length = 2)
names(prediction.plots.Parkinson) <- c("plasma", "CSF")

for(m in 1:length(prediction.plots.Parkinson)){
  prediction.plots.Parkinson[[m]] <- vector(mode = "list", length = length(coefficients.longitudinal.Parkinson[[m]][["gen2"]]))
  names(prediction.plots.Parkinson[[m]]) <- names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])
}

for(m in 1:length(prediction.plots.Parkinson)){
  for(k in 1:length(prediction.plots.Parkinson[[m]])){
    if(!(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k] %in% c("Bladder", "Thyroid"))){
      
      plot.df <- cbind(data.frame(
        predicted.ages = predicted.ages.PD[[m]][["gen2"]][["predicted"]][[k]]
      ), PD_annotation_list[[m]]
      )
      
      ### Only select the last visit here ###
      plot.df <- plot.df %>% arrange(-visit_month)
      plot.df <- plot.df[!duplicated(plot.df$participant_id),] # see prepare_CSF_models.R
      # Remove empty age values
      plot.df <- plot.df[!is.na(plot.df$age_current),]
      # Turn sex into a factor
      plot.df$sex <- factor(plot.df[, "sex.text"])
      
      summary <- summary(lm(predicted.ages ~ age_current, data = plot.df))
      resid <- resid(lm(predicted.ages ~ age_current, data = plot.df))
      metrics <- paste0("r = ", round(sign(summary$coefficients["age_current", "Estimate"])*sqrt(summary$r.squared), 2), ", r² = ", round(summary$r.squared, 2), ", MAE = ", round(mean(abs(resid)), 2))
      
      prediction.plots.Parkinson[[m]][[k]] <- ggplot_gtable(ggplot_build(
        ggplot(data = plot.df, aes(x = age_current, y = predicted.ages)) +
          geom_point(size = 1, aes(color = sex, shape = sex, alpha = 0.05)) + # aes(color = Gender, shape = Sequence)
          scale_color_manual(values=c("magenta", "#4169E1")) +
        geom_smooth(method = MASS::rlm , color="turquoise", fill="turquoise", se=TRUE) + 
          ggtitle(paste0(paste(toupper(substr(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k], 1, 1)), substr(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k], 2, nchar(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k])), sep=""), " model in ", names(prediction.plots.Parkinson)[m], ": ", sum(coefficients.longitudinal.Parkinson[[m]][["gen2"]][[k]] != 0), " / ", length(coefficients.longitudinal.Parkinson[[m]][["gen2"]][[k]]), " proteins", "\n", metrics)) +
          xlab(paste0("Chronological age at last visit")) +
          ylab("Predicted rel. log(mortality hazard)") +
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
            strip.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
      ))
      
      tmp.male <- kMeansEqual(kdat = plot.df[plot.df$sex == "Male", ], 
                              k = floor(nrow(plot.df[plot.df$sex == "Male", ])/5),
                              columns = c("age_current", "predicted.ages"),
                              add.centers = FALSE
      )
      tmp.female <- kMeansEqual(kdat = plot.df[plot.df$sex == "Female", ], 
                                k = floor(nrow(plot.df[plot.df$sex == "Female", ])/5),
                                columns = c("age_current", "predicted.ages"),
                                add.centers = FALSE
      )
      
      tmp.male <- tmp.male %>% group_by(assigned, sex) %>% summarize(
        `N participants` = n(),
        `Avg. chronological age at last visit` = mean(age_current, na.rm = TRUE), 
        `Avg. predicted rel. log(mortality hazard)` = mean(predicted.ages, na.rm = TRUE)
      ) 
      
      tmp.female <- tmp.female %>% group_by(assigned, sex) %>% summarize(
        `N participants` = n(),
        `Avg. chronological age at last visit` = mean(age_current, na.rm = TRUE), 
        `Avg. predicted rel. log(mortality hazard)` = mean(predicted.ages, na.rm = TRUE)
      ) 
      
      tmp <- rbind(tmp.male, tmp.female) %>% arrange(`Avg. chronological age at last visit`)
      tmp <- tmp[, c("N participants", "sex", "Avg. chronological age at last visit", "Avg. predicted rel. log(mortality hazard)")]
      colnames(tmp)[2] <- "Sex"
      
      attr(prediction.plots.Parkinson[[m]][[k]], "data") <- as.data.frame(tmp)
      rm(tmp.male, tmp.female, tmp)
      
    }
  }
}

plot(prediction.plots.Parkinson[["plasma"]]$Conventional)
plot(prediction.plots.Parkinson[["CSF"]]$Conventional)
plot(prediction.plots.Parkinson[["plasma"]]$Brain)
plot(prediction.plots.Parkinson[["CSF"]]$Brain)

### Get dimensions ###
df <- PD_annotation_list[[m]]
### Only select the last visit here ###
df <- df %>% arrange(-visit_month)
df <- df[!duplicated(df$participant_id),] # see prepare_CSF_models.R
# Remove empty age values
df <- df[!is.na(df$age_current),]
dim(df)
# 212  31 # 212 participants at last visit, for both plasma and CSF

### Statistics ###
plot.df <- cbind(data.frame(
  predicted.ages = predicted.ages.PD[["plasma"]][["gen2"]][["predicted"]][[24]]
), PD_annotation_list[["plasma"]]
)
plot.df <- PD_annotation_list[["plasma"]] %>% arrange(-visit_month)
# df$delta_age <- df$visit_month/12
plot.df <- plot.df[!duplicated(plot.df$participant_id),]
plot.df <- df[!is.na(plot.df$age_current),]
plot.df$sex <- factor(plot.df[, "sex.text"])
length(plot.df$participant_id)
# 212 # same as CSF

plot.df <- cbind(data.frame(
  predicted.ages = predicted.ages.PD[["CSF"]][["gen2"]][["predicted"]][[24]]
), PD_annotation_list[["CSF"]]
)
plot.df <- PD_annotation_list[["CSF"]] %>% arrange(-visit_month)
# df$delta_age <- df$visit_month/12
plot.df <- plot.df[!duplicated(plot.df$participant_id),]
plot.df <- df[!is.na(plot.df$age_current),]
plot.df$sex <- factor(plot.df[, "sex.text"])
length(plot.df$participant_id)
# 212 # same as plasma


### 2. Make violin plots for associations with Parkinson's disease and MSA ###

violin.Parkinson.model <- vector(mode = "list", length = 2)
names(violin.Parkinson.model) <- c("plasma", "CSF")

for(m in 1:length(violin.Parkinson.model)){
  violin.Parkinson.model[[m]] <- vector(mode = "list", length = length(coefficients.longitudinal.Parkinson[[m]][["gen2"]]))
  names(violin.Parkinson.model[[m]]) <- names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])
}

for(m in 1:length(violin.Parkinson.model)){
  for(k in 1:length(violin.Parkinson.model[[m]])){
    if(!(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k] %in% c("Bladder", "Thyroid"))){
      
      plot.df <- cbind(data.frame(
        predicted.ages = predicted.ages.PD[[m]][["gen2"]][["predicted"]][[k]],
        resid.ages = predicted.ages.PD[[m]][["gen2"]][["residual"]][[k]]
      ), PD_annotation_list[[m]]
      )
      
      ### Only select the last visit here ###
      plot.df <- plot.df %>% arrange(-visit_month)
      # df$delta_age <- df$visit_month/12
      plot.df <- plot.df[!duplicated(plot.df$participant_id),] # see also Prepare_CSF_Models.R
      plot.df <- plot.df[!is.na(plot.df$age_current),]
      plot.df$sex <- factor(plot.df[, "sex.text"])
      
      plot.df$diagnosis <- plot.df[, "diagnosis_latest"]
      plot.df$diagnosis[plot.df$diagnosis %in% c("Idiopathic PD", "Parkinson's Disease")] <- "Parkinson's Disease"
      plot.df$diagnosis[plot.df$diagnosis == "No PD Nor Other Neurological Disorder"] <- "Control"
      
      plot.df$diagnosis[plot.df$diagnosis == "Essential Tremor"] <- "ET"
      plot.df$diagnosis[plot.df$diagnosis == "Parkinson's Disease"] <- "PD"
      plot.df$diagnosis[plot.df$diagnosis == "Multiple System Atrophy"] <- "MSA"
      
      plot.df$diagnosis <- factor(plot.df$diagnosis, levels = c("Control", "ET", "PD", "MSA"))
      
      model <- lm(predicted.ages ~ age_at_visit*sex + diagnosis, data = plot.df)
      
      # Catch the case in which all predicted relative log(mortality hazards) are 0
      if(all(is.na(summary(model)$coef[, "Pr(>|t|)"]))){
        pvals <- NA
      } else{
        pvals <- summary(glht(model, linfct = mcp(diagnosis = c("ET - Control = 0", "PD - Control = 0", "MSA - Control = 0"))))$test$pvalues
      }
      
      height1 <- (max(plot.df$resid.ages, na.rm = TRUE) - min(plot.df$resid.ages, na.rm = TRUE))*0.98+min(plot.df$resid.ages, na.rm = TRUE)

      significance.df <- data.frame(x = c(2, 3, 4), 
                                    y = c(height1, height1, height1),
                                    pval = pvals,
                                    stars = c("n.s.", "n.s.", "n.s."),
                                    size = c(4.5, 4.5, 4.5),
                                    nudge_y = (max(plot.df$resid.ages, na.rm = TRUE) - min(plot.df$resid.ages, na.rm = TRUE))/13
      )
      
      significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
      significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
      significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
      significance.df$stars[which(significance.df$pval < 0.001)] <- "***"
      
      significance.df[(significance.df$stars != "n.s."), "size"] <- 10
      significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- (max(plot.df$resid.ages, na.rm = TRUE) - min(plot.df$resid.ages, na.rm = TRUE))/30
      
      violin.Parkinson.model[[m]][[k]] <- ggplot_gtable(ggplot_build(
        ggplot(plot.df, aes(x = diagnosis, y = resid.ages)) +
          
          geom_quasirandom(cex = 1.5, aes(x = diagnosis, color = diagnosis, alpha = 1)) +
          # If too many points:
          # geom_violin(aes(color = diagnosis, fill = diagnosis)) + # aes(color = Diet)
          
          scale_color_manual(values=c("#469d76",  "#c08e39", "#c93f55", "#924099")) +
          scale_fill_manual(values=c("#469d76",  "#c08e39", "#c93f55", "#924099")) +
          
          geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
          guides(alpha = "none") + # no legend for "alpha"
          
          geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
          
          ggtitle(paste0(paste(toupper(substr(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k], 1, 1)), substr(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k], 2, nchar(names(coefficients.longitudinal.Parkinson[[m]][["gen2"]])[k])), sep=""), " model in ", names(prediction.plots.Parkinson)[m], ": \n", sum(coefficients.longitudinal.Parkinson[[m]][["gen2"]][[k]] != 0), " / ", length(coefficients.longitudinal.Parkinson[[m]][["gen2"]][[k]]), " proteins")) +
          xlab("Diagnosis") + # No label on the horizontal axis
          ylab(paste0("Log(mortality hazard) deviation")) +
          theme_bw() +
          theme(legend.title = element_text(size = 12), # legend title size
            legend.text = element_text(size = 12), # legend text size
            legend.position = "none", # No legend
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(size = 12),
            strip.text.x = element_text(size = 12), # the size of the facet labels
            axis.title.x = element_text(size = 12, color = "black"), 
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
      
      tmp <- data.frame(
        diagnosis = c(rep("Control", length(xlsx.quantiles)), "ET", rep("PD", length(xlsx.quantiles)), "MSA"),
        Quantile = c(paste0(100*xlsx.quantiles, "%"), "50%", paste0(100*xlsx.quantiles, "%"), "50%")
      )
      
      tmp[, "Log(mortality hazard) deviation"] <- c(
        round(quantile(plot.df[(plot.df$diagnosis == "Control"), "resid.ages"], xlsx.quantiles), 2), 
        round(quantile(plot.df[(plot.df$diagnosis == "ET"), "resid.ages"], 0.5), 2), 
        round(quantile(plot.df[(plot.df$diagnosis == "PD"), "resid.ages"], xlsx.quantiles), 2), 
        round(quantile(plot.df[(plot.df$diagnosis == "MSA"), "resid.ages"], 0.5), 2)
      )
      
      attr(violin.Parkinson.model[[m]][[k]], "data") <- tmp
      
    }
  }
}

plot(violin.Parkinson.model[["plasma"]][["Conventional"]])
plot(violin.Parkinson.model[["plasma"]][["Brain"]])
plot(violin.Parkinson.model[["CSF"]][["Conventional"]])
plot(violin.Parkinson.model[["CSF"]][["Brain"]])

### Statistics ###

table(plot.df$diagnosis)
# Control  ET    PD      MSA 
# 90       2     118       2 

m <- "plasma"
m <- "CSF"
k <- "Conventional"
k <- "Brain"

plot.df <- cbind(data.frame(
  predicted.ages = predicted.ages.PD[[m]][["gen2"]][["predicted"]][[k]],
  resid.ages = predicted.ages.PD[[m]][["gen2"]][["residual"]][[k]]
), PD_annotation_list[[m]]
)

### Only select the last visit here ###
plot.df <- plot.df %>% arrange(-visit_month)
# df$delta_age <- df$visit_month/12
plot.df <- plot.df[!duplicated(plot.df$participant_id),] # see also Prepare_CSF_Models.R
plot.df <- plot.df[!is.na(plot.df$age_current),]
plot.df$sex <- factor(plot.df[, "sex.text"])

plot.df$diagnosis <- plot.df[, "diagnosis_latest"]
plot.df$diagnosis[plot.df$diagnosis %in% c("Idiopathic PD", "Parkinson's Disease")] <- "Parkinson's Disease"
plot.df$diagnosis[plot.df$diagnosis == "No PD Nor Other Neurological Disorder"] <- "Control"

plot.df$diagnosis[plot.df$diagnosis == "Essential Tremor"] <- "ET"
plot.df$diagnosis[plot.df$diagnosis == "Parkinson's Disease"] <- "PD"
plot.df$diagnosis[plot.df$diagnosis == "Multiple System Atrophy"] <- "MSA"

plot.df$diagnosis <- factor(plot.df$diagnosis, levels = c("Control", "ET", "PD", "MSA"))### plasma - conventional ###

summary(glht(lm(predicted.ages ~ age_at_visit*sex + diagnosis, data = plot.df), linfct = mcp(diagnosis = c("ET - Control = 0", "PD - Control = 0", "MSA - Control = 0"))))
### plasma - conventional ###
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: lm(formula = predicted.ages ~ age_at_visit * sex + diagnosis, 
#         data = plot.df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)  
# ET - Control == 0  -0.05187    0.54045  -0.096   0.9995  
# PD - Control == 0   0.24267    0.10706   2.267   0.0712 .
# MSA - Control == 0  0.45022    0.54410   0.827   0.7910  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

### plasma - brain ###
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: lm(formula = predicted.ages ~ age_at_visit * sex + diagnosis, 
#         data = plot.df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# ET - Control == 0   -0.3974     0.6012  -0.661    0.880
# PD - Control == 0    0.1092     0.1191   0.917    0.735
# MSA - Control == 0  -0.2134     0.6053  -0.353    0.979
# (Adjusted p values reported -- single-step method)

### CSF - conventional ###
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: lm(formula = predicted.ages ~ age_at_visit * sex + diagnosis, 
#         data = plot.df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)  
# ET - Control == 0    0.2199     0.5872   0.375   0.9748  
# PD - Control == 0    0.1013     0.1163   0.871   0.7645  
# MSA - Control == 0   1.2685     0.5911   2.146   0.0953 .
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

### CSF - brain ###
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: lm(formula = predicted.ages ~ age_at_visit * sex + diagnosis, 
#         data = plot.df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)   
# ET - Control == 0  -0.03229    0.59803  -0.054  0.99992   
# PD - Control == 0   0.14408    0.11843   1.217  0.53161   
# MSA - Control == 0  1.87292    0.60203   3.111  0.00637 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

summary(glht(lm(predicted.ages ~ age_at_visit*sex + diagnosis, data = plot.df), linfct = mcp(diagnosis = c("ET - Control = 0", "PD - Control = 0", "MSA - Control = 0"))))$test$pvalues

