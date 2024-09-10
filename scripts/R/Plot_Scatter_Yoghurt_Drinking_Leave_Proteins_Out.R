### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
df.test.list <- readRDS(file = paste0(rds.dir, "df_test_list.rds"))

### 1. We are only interested in the associations of these organs with these outcomes ###

outcomes <- c("Full fat yoghurt", "Drinking")

coefs <- vector(mode = "list", length = 2)
names(coefs) <- c("Artery", "Intestine")

coefs[["Artery"]] <- read.csv(paste0(dir.gen2.models, "Artery", "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
coefs[["Intestine"]] <- read.csv(paste0(dir.gen2.models, "Intestine", "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)

### 2. Perform the analysis ###

leave.protein.out.yoghurt.drinking.list <- vector(mode = "list", length = length(organ.proteins))
names(leave.protein.out.yoghurt.drinking.list) <- names(organ.proteins)

for(k in 1:length(leave.protein.out.yoghurt.drinking.list)){
  if(names(organ.proteins)[k] != "Bladder"){
    leave.protein.out.yoghurt.drinking.list[[k]] <- vector(mode = "list", length = length(outcomes))
    names(leave.protein.out.yoghurt.drinking.list[[k]]) <- outcomes
    
    leave.protein.out.yoghurt.drinking.list[[k]] <- lapply(leave.protein.out.yoghurt.drinking.list[[k]], function(x){
      x <- data.frame(protein = organ.proteins[[k]],
                      Estimate = NA,
                      `Std. Error` = NA,
                      `t value` = NA, 
                      df = NA,
                      `Pr(>|t|)` = NA, 
                      ll = NA, 
                      ul = NA, 
                      LRTstat = NA, 
                      pLRTA = NA, 
                      pLRTB = NA, 
                      pLRTAB = NA, 
                      check.names = FALSE)
      return(x)
    })
  }
}

tmp_bd_annotation_0 <- olink_bd_annotation_0

for(k in 1:length(organ.proteins)){
  # if(names(organ.proteins)[k] != "Bladder"){
  if(names(organ.proteins)[k] %in% c("Artery", "Intestine")){
    coef.tmp <- read.csv(paste0(dir.gen2.models, names(organ.proteins)[k], "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    
    for(j in 1:length(coef.tmp)){
      
      tmp <- vector(mode = "list", length = 5)
      for(m in 1:5){
        tmp[[m]] <- rowSums(sweep(as.matrix(df.test.list[[m]][, colnames(coef.tmp)[-j], drop = FALSE]), MARGIN=2, unlist(coef.tmp[m, -j]), `*`))
      }
      
      predicted.ages.oof.minus.j <- (do.call("c", tmp))[as.character(olink_bd_annotation_0$eid)]
      
      for(i in 1:length(outcomes)){
        
        df <- data.frame(
          eid = olink_bd_annotation_0$eid,
          status_os = olink_bd_annotation_0[, paste0("status_", names(general_hazard_outcomes)[i])],
          time = olink_bd_annotation_0[, paste0("time_", names(general_hazard_outcomes)[i])],
          chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
          biological_age_full = predicted.ages.oof[["gen2"]][["predicted"]][[k]],
          biological_age_noJ = predicted.ages.oof.minus.j,
          sex = olink_bd_annotation_0$p31,
          drinking = olink_bd_annotation_0$p26030_i0, 
          deprivation = olink_bd_annotation_0$p22189, # Townsend deprivation index
          centre = factor(olink_bd_annotation_0$p54_i0),
          IPAQ_activity = factor(olink_bd_annotation_0$p22032_i0, level = c("low", "moderate", "high"), ordered = FALSE),
          smoking_status = factor(olink_bd_annotation_0$p20116_i0, levels = c("Prefer not to answer", "Never", "Previous", "Current"), ordered = FALSE)
        )
        
        if(outcomes[i] == "Full fat yoghurt"){
          
          tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i0 == "No" | tmp_bd_annotation_0$p100026_i0 == "No"), paste0("p26096_i0")] <- NA
          tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i1 == "No" | tmp_bd_annotation_0$p100026_i1 == "No"), paste0("p26096_i1")] <- NA
          tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i2 == "No" | tmp_bd_annotation_0$p100026_i2 == "No"), paste0("p26096_i2")] <- NA
          tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i3 == "No" | tmp_bd_annotation_0$p100026_i3 == "No"), paste0("p26096_i3")] <- NA
          tmp_bd_annotation_0[which(tmp_bd_annotation_0$p100020_i4 == "No" | tmp_bd_annotation_0$p100026_i4 == "No"), paste0("p26096_i4")] <- NA
          
          # All food questionnaires were sent out after the first and before the third visit!
          df$food <- rowMeans(tmp_bd_annotation_0[, paste0("p26096_i", 0:4)], na.rm = TRUE)
          
          # Both need to be fit with x = TRUE, so that the x matrix is returned in component x!
          # model_full <- lm(biological_age_full ~ chronological_age * sex + food + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          # model_noJ <- lm(biological_age_noJ ~ chronological_age * sex + food + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          # parm <- "food"
          
          model_full <- lm(food ~ chronological_age * sex + biological_age_full + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          model_noJ <- lm(food ~ chronological_age * sex + biological_age_noJ + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          parm <- "biological_age_noJ"
          
        } else if(outcomes[i] == "Drinking"){
          
          # model_full <- lm(biological_age_full ~ chronological_age * sex + drinking + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          # model_noJ <- lm(biological_age_noJ ~ chronological_age * sex + drinking + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          # parm <- "drinking"
          
          model_full <- lm(drinking ~ chronological_age * sex + biological_age_full + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          model_noJ <- lm(drinking ~ chronological_age * sex + biological_age_noJ + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          parm <- "biological_age_noJ"
        }
        
        lr.object <- vuongtest(object1 = model_full, object2 = model_noJ, nested = FALSE)
        summary.vals <- tryCatch(summary(model_noJ)$coef[parm, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
        
        leave.protein.out.yoghurt.drinking.list[[k]][[i]][(leave.protein.out.yoghurt.drinking.list[[k]][[i]]$protein == colnames(coef.tmp)[j]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model_noJ, parm = parm), sign(model_noJ$coef[parm])*sqrt(summary(model_noJ)$r.squared))
        leave.protein.out.yoghurt.drinking.list[[k]][[i]][(leave.protein.out.yoghurt.drinking.list[[k]][[i]]$protein == colnames(coef.tmp)[j]), "df"] <- model_noJ$df.residual
        leave.protein.out.yoghurt.drinking.list[[k]][[i]][(leave.protein.out.yoghurt.drinking.list[[k]][[i]]$protein == colnames(coef.tmp)[j]), c("LRTstat", "pLRTA", "pLRTB", "pLRTAB")] <- c(lr.object$LRTstat, lr.object$p_LRT$A, lr.object$p_LRT$B, 2*min(lr.object$p_LRT$A, lr.object$p_LRT$B))
        
      }
      
      
    }
  }
}

lr.object
# Model 1 
# Class: lm 
# Call: lm(formula = biological_age_full ~ chronological_age * sex + ...
#          
#          Model 2 
#          Class: lm 
#          Call: lm(formula = biological_age_noJ ~ chronological_age * sex + food + ...
#                   
#                   Variance test 
#                   H0: Model 1 and Model 2 are indistinguishable 
#                   H1: Model 1 and Model 2 are distinguishable 
#                   w2 = 0.174,   p = <2e-16
#                   
#                   Non-nested likelihood ratio test 
#                   H0: Model fits are equal for the focal population 
#                   H1A: Model 1 fits better than Model 2 
#                   z = 22.761,   p = <2e-16
#                   H1B: Model 2 fits better than Model 1 
#                   z = 22.761,   p = 1

# pLRTA small and LRTstat positive means that the full model fits better than the model in which the specific protein is removed, as one would expect
# pLRTAB is for the two-sided test
# Hence most appropriate: signed -log10(pLRTAB)

# With Bioage as a predictor:
# saveRDS(leave.protein.out.yoghurt.drinking.list, file = paste0(rds.dir, "leave_protein_out_yoghurt_drinking_list.rds"))
leave.protein.out.yoghurt.drinking.list <- readRDS(file = paste0(rds.dir, "leave_protein_out_yoghurt_drinking_list.rds"))

leave.protein.out.yoghurt.drinking.list[["Artery"]][["Full fat yoghurt"]]
leave.protein.out.yoghurt.drinking.list[["Intestine"]][["Full fat yoghurt"]]
leave.protein.out.yoghurt.drinking.list[["Artery"]][["Drinking"]]
leave.protein.out.yoghurt.drinking.list[["Intestine"]][["Drinking"]]

# Inspect the distributions
hist(-log10(leave.protein.out.yoghurt.drinking.list[["Artery"]][["Full fat yoghurt"]]$pLRTAB)*sign(leave.protein.out.yoghurt.drinking.list[["Artery"]][["Full fat yoghurt"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.yoghurt.drinking.list[["Intestine"]][["Full fat yoghurt"]]$pLRTAB)*sign(leave.protein.out.yoghurt.drinking.list[["Intestine"]][["Full fat yoghurt"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.yoghurt.drinking.list[["Artery"]][["Drinking"]]$pLRTAB)*sign(leave.protein.out.yoghurt.drinking.list[["Artery"]][["Drinking"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.yoghurt.drinking.list[["Intestine"]][["Drinking"]]$pLRTAB)*sign(leave.protein.out.yoghurt.drinking.list[["Intestine"]][["Drinking"]]$LRTstat), breaks = 100)


scatterplot.yoghurt.drinking.leaveout <- vector(mode = "list", length = 4)
names(scatterplot.yoghurt.drinking.leaveout) <- c("Artery model - full fat yoghurt",
                                                "Intestine model - full fat yoghurt", 
                                                "Artery model - drinking", 
                                                "Intestine model - drinking")

coef.index <- c(1, 2, 1, 2)
outcome.index <- c(1, 1, 2, 2)

for(j in 1:length(scatterplot.yoghurt.drinking.leaveout)){
  
  plot.df <- data.frame(
    protein = organ.proteins[[names(coefs)[coef.index[j]]]], 
    `Minimum coefficient` = apply(coefs[[coef.index[j]]], 2, min),
    `Median coefficient` = apply(coefs[[coef.index[j]]], 2, median),
    `Maximum coefficient` = apply(coefs[[coef.index[j]]], 2, max),
    `Signed -log10(p-value) likelihood ratio test` = -log10(leave.protein.out.yoghurt.drinking.list[[names(coefs)[coef.index[j]]]][[outcomes[outcome.index[j]]]]$pLRTAB)*sign(leave.protein.out.yoghurt.drinking.list[[names(coefs)[coef.index[j]]]][[outcomes[outcome.index[j]]]]$LRTstat),
    pval = leave.protein.out.yoghurt.drinking.list[[names(coefs)[coef.index[j]]]][[outcomes[outcome.index[j]]]]$pLRTAB,
    check.names = FALSE
  )
  
  # Exclude proteins that were not included in the model to begin with
  plot.df <- plot.df[!((plot.df$`Minimum coefficient` == 0) & (plot.df$`Maximum coefficient` == 0)), ]
  
  plot.df$qval <- p.adjust(plot.df$pval, method = "BH")
  
  plot.df$Significance <- "Not significant"
  plot.df$Significance[((plot.df$`Signed -log10(p-value) likelihood ratio test` < 0) & (plot.df$qval < 0.05))] <- "Negative effect & q-value < 5%"
  plot.df$Significance[((plot.df$`Signed -log10(p-value) likelihood ratio test` > 0) & (plot.df$qval < 0.05))] <- "Positive effect & q-value < 5%"
  plot.df$Significance <- factor(plot.df$Significance, levels = c("Not significant", "Negative effect & q-value < 5%", "Positive effect & q-value < 5%"))
  
  plot.df <- plot.df %>% arrange(-abs(`Signed -log10(p-value) likelihood ratio test`))
  plot.df$label <- ""
  plot.df$label[1:min(c(10, nrow(plot.df)))] <- plot.df$protein[1:min(c(10, nrow(plot.df)))]
  
  scatterplot.yoghurt.drinking.leaveout[[j]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = `Signed -log10(p-value) likelihood ratio test`, y = `Median coefficient`)) +
      geom_pointrange(aes(ymin = `Minimum coefficient`, ymax = `Maximum coefficient`, color = Significance), size = 0.7)  +
      scale_color_manual(values = alpha(c("black", "#ff0000ff", "#00be64ff"), 0.5), breaks = levels(plot.df$Significance)) +
      ggtitle("") + # ggtitle("Chronological age adjusted for sex") +
      ggtitle(names(scatterplot.yoghurt.drinking.leaveout)[j]) +
      geom_vline(xintercept = 0, colour = "darkgrey") +
      geom_hline(yintercept = 0, colour = "darkgrey") +
      geom_text_repel(aes(label = label)) +
      theme_bw() +
      theme(
        legend.position = "none", # No legend
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=12),
        axis.title.x = element_text(size = 12, color = "black"),
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
  
  attr(scatterplot.yoghurt.drinking.leaveout[[j]], "data") <- plot.df
  
}



plot(scatterplot.yoghurt.drinking.leaveout[[1]])
plot(scatterplot.yoghurt.drinking.leaveout[[2]])
plot(scatterplot.yoghurt.drinking.leaveout[[3]])
plot(scatterplot.yoghurt.drinking.leaveout[[4]])




