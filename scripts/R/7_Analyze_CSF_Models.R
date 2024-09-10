
### This comes from 6_Prepare_Parkinson_MSA_Data.R
PD_olink_data_list_raw <- readRDS(file = paste0(rds.dir, "PD_olink_data_list_raw.rds"))
PD_annotation_list <- readRDS(file = paste0(rds.dir, "PD_annotation_list.rds"))

### This comes from Prepare_CSF_Models.R ###
organ.proteins.CSF <- readRDS(paste0(rds.dir, "organ_proteins_CSF.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

# Only needed for export to Excel
xlsx.quantiles <- c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99, 0.999, 1)

# For these ones, we need the leave-one-out-predictions
df.loo.imputed.data.CSF <- read_csv_fast(paste0(dir.loo.imputed.data.CSF, "test_data_LOO_imputed.csv"), header = TRUE, sep = ",", quote = "\"")

# An indicator to show if an observation corresponds to a first visit
first.visit <- PD_annotation_list[["CSF"]]$GUID %in% rownames(df.loo.imputed.data.CSF)
# An indicator to show if an observation corresponds to the last visit

### Only select the last visit here ###
tmp <- PD_annotation_list[["CSF"]] %>% arrange(-visit_month)
tmp <- tmp[!duplicated(tmp$participant_id),]
tmp <- tmp[!is.na(tmp$age_current),]
last.visit <- PD_annotation_list[["CSF"]]$GUID %in% tmp$GUID

### Sanity check ###
all(first.visit == (PD_annotation_list[["CSF"]]$visit_month == 0))
# TRUE

# impute_olink was sourced in 1_Prepare_Data.R (Note: do NOT use "last.visit" here!)
df.other.imputed.data.CSF <- impute_olink(PD_olink_data_list_raw[["CSF"]][which(!first.visit),])

coefficients.CSF <- vector(mode = "list", length = length(organ.proteins.CSF))
names(coefficients.CSF) <- names(organ.proteins.CSF)

predicted.ages.CSF <- vector(mode = "list", length = 2)
names(predicted.ages.CSF) <- c("predicted", "residual")
predicted.ages.CSF <- lapply(predicted.ages.CSF, function(y){
  y <- vector(mode = "list", length = length(organ.proteins.CSF))
  names(y) <- names(organ.proteins.CSF)
  y <- lapply(y, function(z){
    z <- rep(NA, nrow(PD_annotation_list[["CSF"]]))
    names(z) <- PD_annotation_list[["CSF"]]$GUID
    return(z)
  })
  return(y)
})

for(k in 1:length(organ.proteins.CSF)){
  if(!(names(organ.proteins.CSF)[k] %in% c("Bladder", "Thyroid"))){
    
    ### For the first visit, we use leave-one-out predictions ###
    coefficients.CSF[[k]] <- read.csv(paste0(dir.gen1.models.CSF, names(organ.proteins.CSF)[k], "_coefs.csv"), header = TRUE, check.names = FALSE)
    coef.tmp <- coefficients.CSF[[k]][, -1, drop = FALSE]
    intercept <- coefficients.CSF[[k]][, "Intercept"]
    
    predicted.ages.CSF[["predicted"]][[k]][which(first.visit)] <- intercept + rowSums(df.loo.imputed.data.CSF[, colnames(coef.tmp), drop = FALSE] * coef.tmp)
    predicted.ages.CSF[["residual"]][[k]][which(first.visit)] <- resid(lm(predicted.ages.CSF[["predicted"]][[k]][which(first.visit)] ~ PD_annotation_list[["CSF"]]$age_at_visit[which(first.visit)]))
    
    ### For the last visit, we use the first leave-one-out cross-validation model ###
    coef.tmp <- unlist(coefficients.CSF[[k]][1, -1, drop = FALSE])
    intercept <- coefficients.CSF[[k]][1, "Intercept"]
    
    predicted.ages.CSF[["predicted"]][[k]][which(!first.visit)] <- intercept + rowSums(sweep(as.matrix(df.other.imputed.data.CSF[, names(coef.tmp), drop = FALSE]), MARGIN=2, coef.tmp, `*`))
    predicted.ages.CSF[["residual"]][[k]][which(!first.visit)] <- resid(lm(predicted.ages.CSF[["predicted"]][[k]][which(!first.visit)] ~ PD_annotation_list[["CSF"]]$age_at_visit[which(!first.visit)]))
    
  }
}

# saveRDS(predicted.ages.CSF, file = paste0(rds.dir, "predicted_ages_CSF.rds"))
predicted.ages.CSF <- readRDS(file = paste0(rds.dir, "predicted_ages_CSF.rds"))

# saveRDS(coefficients.CSF, file = paste0(rds.dir, "coefficients_CSF.rds"))
coefficients.CSF <- readRDS(file = paste0(rds.dir, "coefficients_CSF.rds"))

### Correlations with chronological age ###

prediction.plots.CSF <- vector(mode = "list", length = 2)
names(prediction.plots.CSF) <- c("first visit", "last visit") # c("training dataset", "test dataset", "third visit", "fourth visit")

for(i in 1:length(prediction.plots.CSF)){
  prediction.plots.CSF[[i]] <- vector(mode = "list", length = length(organ.proteins.CSF))
  names(prediction.plots.CSF[[i]]) <- names(organ.proteins.CSF)
}

age.visit <- c("first", "last")

for(i in 1:length(age.visit)){
  for(k in 1:length(prediction.plots.CSF[[i]])){
    if(!(names(organ.proteins.CSF)[k] %in% c("Bladder", "Thyroid"))){
      
      if(names(prediction.plots.CSF)[i] == "first visit"){
        plot.df <- data.frame(GUID = PD_annotation_list[["CSF"]][which(first.visit), "GUID"], 
                         biological_age = predicted.ages.CSF[["predicted"]][[k]][which(first.visit)],
                         age_current = PD_annotation_list[["CSF"]][which(first.visit), "age_current"],
                         sex = factor(PD_annotation_list[["CSF"]][which(first.visit), "sex.text"])
        )
      } else if(names(prediction.plots.CSF)[i] == "last visit"){
        plot.df <- data.frame(GUID = PD_annotation_list[["CSF"]][which(last.visit), "GUID"], 
                         biological_age = predicted.ages.CSF[["predicted"]][[k]][which(last.visit)],
                         age_current = PD_annotation_list[["CSF"]][which(last.visit), "age_current"],
                         sex = factor(PD_annotation_list[["CSF"]][which(last.visit), "sex.text"])
        )
      }

      plot.df <- plot.df[!is.na(plot.df$age_current),]
      
      summary <- summary(lm(biological_age ~ age_current, data = plot.df))
      resid <- resid(lm(biological_age ~ age_current, data = plot.df))
      metrics <- paste0("r = ", sprintf("%.2f", round(sign(summary$coefficients["age_current", "Estimate"])*sqrt(summary$r.squared), 2)), ", r² = ", sprintf("%.2f", round(summary$r.squared, 2)), ", MAE = ", sprintf("%.2f", round(mean(abs(resid)), 2)))
      
      prediction.plots.CSF[[i]][[k]] <- ggplot_gtable(ggplot_build(
        ggplot(data = plot.df, aes(x = age_current, y = biological_age)) +
          
          geom_point(size = 1, aes(color = sex, shape = sex, alpha = 0.05)) + 
          scale_color_manual(values=c("magenta", "#4169E1")) +

          geom_smooth(method = MASS::rlm , color="turquoise", fill="turquoise", se=TRUE) + 
          ggtitle(paste0(paste(toupper(substr(names(organ.proteins.CSF)[k], 1, 1)), substr(names(organ.proteins.CSF)[k], 2, nchar(names(organ.proteins.CSF)[k])), sep=""), " CSF model: ", names(prediction.plots.CSF)[i], ": ", sum(coefficients.CSF[[k]][1, -1] != 0), " / ", length(organ.proteins.CSF[[k]]), " proteins", "\n", metrics)) +
          xlab(paste0("Chronological age at ", age.visit[i], " visit")) +
          ylab("Predicted age") +
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
                              columns = c("age_current", "biological_age"),
                              add.centers = FALSE
      )
      tmp.female <- kMeansEqual(kdat = plot.df[plot.df$sex == "Female", ], 
                                k = floor(nrow(plot.df[plot.df$sex == "Female", ])/5),
                                columns = c("age_current", "biological_age"),
                                add.centers = FALSE
      )
      
      tmp.male <- tmp.male %>% group_by(assigned, sex) %>% summarize(
        `N participants` = n(),
        `Avg. chronological age at TMP visit` = mean(age_current, na.rm = TRUE), 
        `Avg. predicted age` = mean(biological_age, na.rm = TRUE)
      ) 
      
      tmp.female <- tmp.female %>% group_by(assigned, sex) %>% summarize(
        `N participants` = n(),
        `Avg. chronological age at TMP visit` = mean(age_current, na.rm = TRUE), 
        `Avg. predicted age` = mean(biological_age, na.rm = TRUE)
      ) 
      
      tmp <- rbind(tmp.male, tmp.female) %>% arrange(`Avg. chronological age at TMP visit`)
      tmp <- tmp[, c("N participants", "sex", "Avg. chronological age at TMP visit", "Avg. predicted age")]
      colnames(tmp)[2] <- "Sex"
      colnames(tmp)[3] <- paste0("Chronological age at ", age.visit[i], " visit")
      
      attr(prediction.plots.CSF[[i]][[k]], "data") <- as.data.frame(tmp)
      rm(tmp.male, tmp.female, tmp)
      
    }
  }
}

plot(prediction.plots.CSF[["first visit"]]$Conventional)
plot(prediction.plots.CSF[["last visit"]]$Conventional)

plot(prediction.plots.CSF[["first visit"]]$Brain)
plot(prediction.plots.CSF[["last visit"]]$Brain)

violin.model.CSF <- vector(mode = "list", length = length(organ.proteins.CSF))
names(violin.model.CSF) <- names(organ.proteins.CSF)

for(k in 1:length(violin.model.CSF)){
  if(!(names(organ.proteins.CSF)[k] %in% c("Bladder", "Thyroid"))){
    
    plot.df <- data.frame(
               GUID = PD_annotation_list[["CSF"]][which(last.visit), "GUID"], 
               diagnosis = PD_annotation_list[["CSF"]][which(last.visit), "diagnosis_latest"], 
               biological_age = predicted.ages.CSF[["predicted"]][[k]][which(last.visit)],
               residual_age = predicted.ages.CSF[["residual"]][[k]][which(last.visit)],
               age_current = PD_annotation_list[["CSF"]][which(last.visit), "age_current"],
               sex = factor(PD_annotation_list[["CSF"]][which(last.visit), "sex.text"])
    )
    
    plot.df$diagnosis[plot.df$diagnosis %in% c("Idiopathic PD", "Parkinson's Disease")] <- "Parkinson's Disease"
    plot.df$diagnosis[plot.df$diagnosis == "No PD Nor Other Neurological Disorder"] <- "Control"
    
    plot.df$diagnosis[plot.df$diagnosis == "Essential Tremor"] <- "ET"
    plot.df$diagnosis[plot.df$diagnosis == "Parkinson's Disease"] <- "PD"
    plot.df$diagnosis[plot.df$diagnosis == "Multiple System Atrophy"] <- "MSA"
    
    plot.df$diagnosis <- factor(plot.df$diagnosis, levels = c("Control", "ET", "PD", "MSA"))
    
    summary <- summary(glht(lm(biological_age ~ age_current*sex + diagnosis, data = plot.df), linfct = mcp(diagnosis = c("ET - Control = 0", "PD - Control = 0", "MSA - Control = 0"))))
    
    height1 <- (max(plot.df$residual_age, na.rm = TRUE) - min(plot.df$residual_age, na.rm = TRUE))*0.98+min(plot.df$residual_age, na.rm = TRUE)
    
    significance.df <- data.frame(x = c(2, 3, 4), 
                                  y = c(height1, height1, height1),
                                  pval = summary$test$pvalues,
                                  stars = c("n.s.", "n.s.", "n.s."),
                                  size = c(4.5, 4.5, 4.5),
                                  nudge_y = (max(plot.df$residual_age, na.rm = TRUE) - min(plot.df$residual_age, na.rm = TRUE))/13
    )
    
    significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
    significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
    significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
    significance.df$stars[which(significance.df$pval < 0.001)] <- "***"
    
    significance.df[(significance.df$stars != "n.s."), "size"] <- 10
    significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- (max(plot.df$residual_age, na.rm = TRUE) - min(plot.df$residual_age, na.rm = TRUE))/30
    
    ### Statistics ###
    # k <- 24
    # summary(glht(lm(residual_age~diagnosis, data = plot.df), linfct = mcp(diagnosis = c("`ET` - Control = 0", "PD - Control = 0", "MSA - Control = 0"))))
    # 0.992 0.066 . 0.145  
    
    # k <- 5
    # summary(glht(lm(residual_age~diagnosis, data = plot.df), linfct = mcp(diagnosis = c("`ET` - Control = 0", "PD - Control = 0", "MSA - Control = 0"))))
    # 0.83812  0.01084 * 0.00179 **
    
    violin.model.CSF[[k]] <- ggplot_gtable(ggplot_build(
      ggplot(plot.df, aes(x = diagnosis, y = residual_age)) +
        
        geom_quasirandom(cex = 1.5, aes(x = diagnosis, color = diagnosis, alpha = 1)) +
        scale_color_manual(values=c("#469d76",  "#c08e39", "#c93f55", "#924099")) +
        scale_fill_manual(values=c("#469d76",  "#c08e39", "#c93f55", "#924099")) +
        
        geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
        guides(alpha = "none") + # no legend for "alpha"
        
        geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
        
        ggtitle(paste0(paste(toupper(substr(names(organ.proteins.CSF)[k], 1, 1)), substr(names(organ.proteins.CSF)[k], 2, nchar(names(organ.proteins.CSF)[k])), sep=""), " CSF model: last visit: ", sum(coefficients.CSF[[k]][1, -1] != 0), " / ", length(organ.proteins.CSF[[k]]), " proteins")) +
        xlab("Diagnosis") + # No label on the horizontal axis
        ylab(paste0("Age deviation")) +
        theme_bw() +
        theme(legend.title = element_text(size = 12), # legend title size
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
    
    tmp <- data.frame(
      diagnosis = c(rep("Control", length(xlsx.quantiles)), "ET", rep("PD", length(xlsx.quantiles)), "MSA"),
      Quantile = c(paste0(100*xlsx.quantiles, "%"), "50%", paste0(100*xlsx.quantiles, "%"), "50%")
    )
    
    tmp[, "Age deviation"] <- c(
      round(quantile(plot.df[(plot.df$diagnosis == "Control"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$diagnosis == "ET"), "residual_age"], 0.5), 2), 
      round(quantile(plot.df[(plot.df$diagnosis == "PD"), "residual_age"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$diagnosis == "MSA"), "residual_age"], 0.5), 2)
    )
    
    attr(violin.model.CSF[[k]], "data") <- tmp
    
  }
}


plot(violin.model.CSF[["Conventional"]])
plot(violin.model.CSF[["Brain"]])

plot(prediction.plots.CSF[["first visit"]]$Conventional)
plot(prediction.plots.CSF[["last visit"]]$Conventional)


k <- 24
# k <- 5
plot.df <- data.frame(
  GUID = PD_annotation_list[["CSF"]][which(last.visit), "GUID"], 
  diagnosis = PD_annotation_list[["CSF"]][which(last.visit), "diagnosis_latest"], 
  biological_age = predicted.ages.CSF[["predicted"]][[k]][which(last.visit)],
  residual_age = predicted.ages.CSF[["residual"]][[k]][which(last.visit)],
  age_current = PD_annotation_list[["CSF"]][which(last.visit), "age_current"],
  sex = factor(PD_annotation_list[["CSF"]][which(last.visit), "sex.text"])
)

plot.df$diagnosis[plot.df$diagnosis %in% c("Idiopathic PD", "Parkinson's Disease")] <- "Parkinson's Disease"
plot.df$diagnosis[plot.df$diagnosis == "No PD Nor Other Neurological Disorder"] <- "Control"

plot.df$diagnosis[plot.df$diagnosis == "Essential Tremor"] <- "ET"
plot.df$diagnosis[plot.df$diagnosis == "Parkinson's Disease"] <- "PD"
plot.df$diagnosis[plot.df$diagnosis == "Multiple System Atrophy"] <- "MSA"

plot.df$diagnosis <- factor(plot.df$diagnosis, levels = c("Control", "ET", "PD", "MSA"))

summary(glht(lm(biological_age ~ age_current*sex + diagnosis, data = plot.df), linfct = mcp(diagnosis = c("ET - Control = 0", "PD - Control = 0", "MSA - Control = 0"))))

### For conventional ###
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: lm(formula = biological_age ~ age_current * sex + diagnosis, 
#         data = plot.df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)  
# ET - Control == 0   -1.5643     4.3968  -0.356   0.9782  
# PD - Control == 0    1.9102     0.8707   2.194   0.0851 .
# MSA - Control == 0   4.6853     4.4262   1.059   0.6405  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

### For brain ###
# Simultaneous Tests for General Linear Hypotheses
# 
# Multiple Comparisons of Means: User-defined Contrasts
# 
# 
# Fit: lm(formula = biological_age ~ age_current * sex + diagnosis, 
#         data = plot.df)
# 
# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)    
# ET - Control == 0     1.276      5.856   0.218   0.9948    
# PD - Control == 0     2.769      1.160   2.388   0.0524 .  
# MSA - Control == 0   27.911      5.895   4.734   <1e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)
