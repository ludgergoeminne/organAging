
### This comes from 3_Process_GTEx.R ###
organ.proteins.longitudinal <- readRDS(paste0(rds.dir, "organ_proteins_longitudinal.rds"))

### This comes from 6_Prepare_Parkinson_MSA_Data.R ###
# This one has"plasma" and "CSF" as levels
PD_annotation_list <- readRDS(file = paste0(rds.dir, "PD_annotation_list.rds"))
PD_olink_data_list_imputed <- readRDS(file = paste0(rds.dir, "cleanDat_Olink_imputed_list.rds"))
coefficients.longitudinal.Parkinson <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal_Parkinson.rds"))
predicted.ages.PD <- readRDS(file = paste0(rds.dir, "predicted_ages_PD.rds"))


coefs.brain <- read.csv(paste0(dir.gen2.models.longitudinal, "Brain", "_mortality_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
coefs.brain <- list(coefs.brain[, colnames(coefs.brain) %in% names(coefficients.longitudinal.Parkinson$plasma[["gen2"]]$Brain)], 
                    coefs.brain[, colnames(coefs.brain) %in% names(coefficients.longitudinal.Parkinson$CSF[["gen2"]]$Brain)])
names(coefs.brain) <- names(coefficients.longitudinal.Parkinson)

for(m in 1:length(coefs.brain)){
  coefs.brain[[m]] <- coefs.brain[[m]][, names(coefficients.longitudinal.Parkinson[[m]][["gen2"]]$Brain)]
}

### 1. We are only interested in the associations of these organs with these outcomes ###

outcomes <- c("PD", "MSA") # i
fluids <- c("plasma", "CSF") # m

### 2. Perform the analysis ###

leave.protein.out.PD.MSA.list <- vector(mode = "list", length = length(fluids))
names(leave.protein.out.PD.MSA.list) <- fluids

for(m in 1:length(leave.protein.out.PD.MSA.list)){
  leave.protein.out.PD.MSA.list[[m]] <- vector(mode = "list", length = length(outcomes))
  names(leave.protein.out.PD.MSA.list[[m]]) <- outcomes
  leave.protein.out.PD.MSA.list[[m]] <- lapply(leave.protein.out.PD.MSA.list[[m]], function(y){
    y <- vector(mode = "list", length = length(organ.proteins.longitudinal))
    names(y) <- names(organ.proteins.longitudinal)
    for(k in 1:length(y)){
      if((names(y)[k] %in% c("Brain"))){
        y[[k]] <- data.frame(protein = names(coefs.brain[[m]]),
                             Estimate = NA,
                             `Std. Error` = NA,
                             df = NA,
                             `t value` = NA,
                             `Pr(>|t|)` = NA, 
                             ll = NA, 
                             ul = NA, 
                             r = NA, 
                             LRTstat = NA, 
                             pLRTA = NA, 
                             pLRTB = NA, 
                             pLRTAB = NA, 
                             check.names = FALSE)
      }
    }
    return(y)
  })
}

for(m in 1:length(fluids)){
  
  for(k in 1:length(organ.proteins.longitudinal)){
    # if(names(organ.proteins.longitudinal)[k] != "Bladder"){
    if(names(organ.proteins.longitudinal)[k] %in% c("Brain")){
      coef.tmp <- coefficients.longitudinal.Parkinson[[m]][["gen2"]][[k]]
      
      tmp2 <- PD_olink_data_list_imputed[[m]][ , names(coef.tmp), drop = FALSE]
      
      df.predictions <- cbind(data.frame(
        GUID = rownames(tmp2),
        predicted.ages = predicted.ages.PD[[m]][["gen2"]][["predicted"]][[k]],
        resid.ages = predicted.ages.PD[[m]][["gen2"]][["residual"]][[k]],
        PD_annotation_list[[m]]
      ))
      
      # j: the indicator that loops over leaving every protein out
      for(j in 1:length(coef.tmp)){
        
        df.predictions$predicted.ages.oof.minus.j <- rowSums(sweep(as.matrix(tmp2[, names(coef.tmp)[-j], drop = FALSE]), MARGIN=2, coef.tmp[-j], `*`))
        
        for(i in 1:length(outcomes)){
          
          if(outcomes[i] == "PD"){
            df <- df.predictions[df.predictions$diagnosis %in% c("Control", "Parkinson's Disease"), ]
            df$diagnosis <- factor(df$diagnosis, levels = c("Control", "Parkinson's Disease"))
            
          } else if(outcomes[i] == "MSA"){
            df <- df.predictions[df.predictions$diagnosis %in% c("Control", "Multiple System Atrophy"), ]
            df$diagnosis <- factor(df$diagnosis, levels = c("Control", "Multiple System Atrophy"))
          } else{
            stop("Outcome not defined!")
          }
          
          df$diagnosis <- as.numeric(df$diagnosis)-1
          df$sex.text <- factor(df$sex.text)
          
          model_full <- lm(diagnosis ~ age_at_visit*sex.text + predicted.ages, data = df)
          model_noJ <- lm(diagnosis ~ age_at_visit*sex.text + predicted.ages.oof.minus.j, data = df)
          parm <- "predicted.ages.oof.minus.j"
          
          lr.object <- vuongtest(object1 = model_full, object2 = model_noJ, nested = FALSE)
          summary.vals <- tryCatch(summary(model_noJ)$coef[parm, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
          
          leave.protein.out.PD.MSA.list[[m]][[i]][[k]][(leave.protein.out.PD.MSA.list[[m]][[i]][[k]]$protein == names(coef.tmp)[j]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model_noJ, parm = parm), sign(model_noJ$coef[parm])*sqrt(summary(model_noJ)$r.squared))
          leave.protein.out.PD.MSA.list[[m]][[i]][[k]][(leave.protein.out.PD.MSA.list[[m]][[i]][[k]]$protein == names(coef.tmp)[j]), "df"] <- model_noJ$df.residual
          leave.protein.out.PD.MSA.list[[m]][[i]][[k]][(leave.protein.out.PD.MSA.list[[m]][[i]][[k]]$protein == names(coef.tmp)[j]), c("LRTstat", "pLRTA", "pLRTB", "pLRTAB")] <- c(lr.object$LRTstat, lr.object$p_LRT$A, lr.object$p_LRT$B, 2*min(lr.object$p_LRT$A, lr.object$p_LRT$B))
          
        }
        
        
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
# saveRDS(leave.protein.out.PD.MSA.list, file = paste0(rds.dir, "leave_protein_out_PD_MSA_list.rds"))
leave.protein.out.PD.MSA.list <- readRDS(file = paste0(rds.dir, "leave_protein_out_PD_MSA_list.rds"))

leave.protein.out.PD.MSA.list[["plasma"]][["PD"]][["Brain"]]
leave.protein.out.PD.MSA.list[["CSF"]][["PD"]][["Brain"]]
leave.protein.out.PD.MSA.list[["plasma"]][["MSA"]][["Brain"]]
leave.protein.out.PD.MSA.list[["CSF"]][["MSA"]][["Brain"]]

# Inspect the distributions
hist(-log10(leave.protein.out.PD.MSA.list[["plasma"]][["PD"]][["Brain"]]$pLRTAB)*sign(leave.protein.out.PD.MSA.list[["plasma"]][["PD"]][["Brain"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.PD.MSA.list[["CSF"]][["PD"]][["Brain"]]$pLRTAB)*sign(leave.protein.out.PD.MSA.list[["CSF"]][["PD"]][["Brain"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.PD.MSA.list[["plasma"]][["MSA"]][["Brain"]]$pLRTAB)*sign(leave.protein.out.PD.MSA.list[["plasma"]][["MSA"]][["Brain"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.PD.MSA.list[["CSF"]][["MSA"]][["Brain"]]$pLRTAB)*sign(leave.protein.out.PD.MSA.list[["CSF"]][["MSA"]][["Brain"]]$LRTstat), breaks = 100)


scatterplot.PD.MSA.leaveout <- vector(mode = "list", length = 4)
names(scatterplot.PD.MSA.leaveout) <- c("Plasma - Parkinson's Disease",
                                                  "CSF - Parkinson's Disease", 
                                                  "Plasma - Multiple System Atrophy", 
                                                  "CSF - Multiple System Atrophy")

fluid.index <- c(1, 2, 1, 2)
outcome.index <- c(1, 1, 2, 2)

for(j in 1:length(scatterplot.PD.MSA.leaveout)){
  
  plot.df <- data.frame(
    protein = names(coefs.brain[[names(coefs.brain)[fluid.index[j]]]]), 
    `Minimum coefficient` = apply(coefs.brain[[fluid.index[j]]], 2, min),
    `Median coefficient` = apply(coefs.brain[[fluid.index[j]]], 2, median),
    `Maximum coefficient` = apply(coefs.brain[[fluid.index[j]]], 2, max),
    `Signed -log10(p-value) likelihood ratio test` = -log10(leave.protein.out.PD.MSA.list[[names(coefs.brain)[fluid.index[j]]]][[outcomes[outcome.index[j]]]][["Brain"]]$pLRTAB)*sign(leave.protein.out.PD.MSA.list[[names(coefs.brain)[fluid.index[j]]]][[outcomes[outcome.index[j]]]][["Brain"]]$LRTstat),
    pval = leave.protein.out.PD.MSA.list[[names(coefs.brain)[fluid.index[j]]]][[outcomes[outcome.index[j]]]][["Brain"]]$pLRTAB,
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
  
  scatterplot.PD.MSA.leaveout[[j]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = `Signed -log10(p-value) likelihood ratio test`, y = `Median coefficient`)) +
      geom_pointrange(aes(ymin = `Minimum coefficient`, ymax = `Maximum coefficient`, color = Significance), size = 0.7)  +
      scale_color_manual(values = alpha(c("black", "#ff0000ff", "#00be64ff"), 0.5), breaks = levels(plot.df$Significance)) +
      ggtitle("") + # ggtitle("Chronological age adjusted for sex") +
      ggtitle(names(scatterplot.PD.MSA.leaveout)[j]) +
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
  
  attr(scatterplot.PD.MSA.leaveout[[j]], "data") <- plot.df
  
}

plot(scatterplot.PD.MSA.leaveout[[1]])
plot(scatterplot.PD.MSA.leaveout[[2]])
plot(scatterplot.PD.MSA.leaveout[[3]])
plot(scatterplot.PD.MSA.leaveout[[4]])




