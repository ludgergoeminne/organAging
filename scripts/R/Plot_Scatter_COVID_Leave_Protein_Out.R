
### This comes from 3_Process_GTEx.R ###
organ.proteins.longitudinal <- readRDS(paste0(rds.dir, "organ_proteins_longitudinal.rds"))

### This comes from Plot_Violin_COVID.R ###
coefficients.longitudinal.gen2.COVID <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal_gen2_COVID.rds"))
df.COVID.list <- readRDS(file = paste0(rds.dir, "df_COVID_list.rds"))
Olink_Proteomics.COVID <- readRDS(file = paste0(rds.dir, "Olink_Proteomics_COVID.rds"))

coefs <- vector(mode = "list", length = 2)
names(coefs) <- c("Lung", "Conventional")

coefs[["Lung"]] <- read.csv(paste0(dir.gen2.models.longitudinal, "Lung", "_mortality_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)
coefs[["Conventional"]] <- read.csv(paste0(dir.gen2.models.longitudinal, "Conventional", "_mortality_coefs_GTEx_4x_FC_longitudinal.csv"), header = TRUE, check.names = FALSE)

for(k in 1:length(coefs)){
  coefs[[k]] <- coefs[[k]][, names(coefficients.longitudinal.gen2.COVID[[names(coefs)[k]]])]
}

### Sanity checks ###
all(coefficients.longitudinal.gen2.COVID[["Lung"]] == coefs[["Lung"]][1, ])
# TRUE
all(coefficients.longitudinal.gen2.COVID[["Conventional"]] == coefs[["Conventional"]][1, ])
# TRUE
all(rownames(Olink_Proteomics.COVID) == df.COVID.list[["Conventional"]]$Unique.ID)
# TRUE
all(rownames(Olink_Proteomics.COVID) == df.COVID.list[["Lung"]]$Unique.ID)
# TRUE

### 1. We are only interested in the associations of these organs with these outcomes ###

outcomes <- c("COVID", "Acuity.max") # i

### 2. Perform the analysis ###

leave.protein.out.COVID.list <- vector(mode = "list", length = length(outcomes))
names(leave.protein.out.COVID.list) <- outcomes

leave.protein.out.COVID.list <- lapply(leave.protein.out.COVID.list, function(x){
  x <- vector(mode = "list", length = length(organ.proteins.longitudinal))
    names(x) <- names(organ.proteins.longitudinal)
   for(k in 1:length(x)){
      if((names(x)[k] %in% names(coefs))){
      x[[k]] <- data.frame(protein = names(coefs[[names(x)[k]]]),
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
    return(x)
  })

for(k in 1:length(organ.proteins.longitudinal)){
  # if(names(organ.proteins.longitudinal)[k] != "Bladder"){
  if(names(organ.proteins.longitudinal)[k] %in% names(coefs)){
    
    coef.tmp <- coefficients.longitudinal.gen2.COVID[[k]]
    
    tmp2 <- Olink_Proteomics.COVID[df.COVID.list[[k]]$Day == 0, colnames(Olink_Proteomics.COVID) %in% names(coef.tmp), drop = FALSE]
    df.0 <- df.COVID.list[[k]][df.COVID.list[[k]]$Day == 0, ]
    df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
    df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
    tmp2 <- tmp2[ , names(coef.tmp), drop = FALSE]
    
    # j: the indicator that loops over leaving every protein out
    for(j in 1:length(coef.tmp)){
      
      df.0$predicted.ages.oof.minus.j <- rowSums(sweep(as.matrix(tmp2[, names(coef.tmp)[-j], drop = FALSE]), MARGIN=2, coef.tmp[-j], `*`))
      
      for(i in 1:length(outcomes)){
        
        if(outcomes[i] == "COVID"){
          conversion <- c(1,2)
          names(conversion) <- c("Neg.", "Pos.")
          df.0$outcome <- as.numeric(factor(names(conversion[df.0[, outcomes[i]]]), levels = names(conversion)))-1
 
        } else if(outcomes[i] == "Acuity.max"){
          conversion <- 1:5
          names(conversion) <- c("Discharged", "Hospitalized - O2", "Hospitalized + O2", "Intubated/ventilated", "Death")
          df.0$outcome <- as.numeric(factor(names(conversion[df.0[, outcomes[i]]]), levels = names(conversion)))-1
        } else{
          stop("Outcome not defined!")
        }
        
        model_full <- lm(outcome ~ age.numeric + predicted.ages, data = df.0)
        model_noJ <- lm(outcome ~ age.numeric + predicted.ages.oof.minus.j, data = df.0)
        parm <- "predicted.ages.oof.minus.j"
        
        lr.object <- vuongtest(object1 = model_full, object2 = model_noJ, nested = FALSE)
        summary.vals <- tryCatch(summary(model_noJ)$coef[parm, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
        
        leave.protein.out.COVID.list[[i]][[k]][(leave.protein.out.COVID.list[[i]][[k]]$protein == names(coef.tmp)[j]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model_noJ, parm = parm), sign(model_noJ$coef[parm])*sqrt(summary(model_noJ)$r.squared))
        leave.protein.out.COVID.list[[i]][[k]][(leave.protein.out.COVID.list[[i]][[k]]$protein == names(coef.tmp)[j]), "df"] <- model_noJ$df.residual
        leave.protein.out.COVID.list[[i]][[k]][(leave.protein.out.COVID.list[[i]][[k]]$protein == names(coef.tmp)[j]), c("LRTstat", "pLRTA", "pLRTB", "pLRTAB")] <- c(lr.object$LRTstat, lr.object$p_LRT$A, lr.object$p_LRT$B, 2*min(lr.object$p_LRT$A, lr.object$p_LRT$B))
        
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
# saveRDS(leave.protein.out.COVID.list, file = paste0(rds.dir, "leave_protein_out_COVID_list.rds"))
leave.protein.out.COVID.list <- readRDS(file = paste0(rds.dir, "leave_protein_out_COVID_list.rds"))

leave.protein.out.COVID.list[["COVID"]][["Lung"]]
leave.protein.out.COVID.list[["COVID"]][["Conventional"]]
leave.protein.out.COVID.list[["Acuity.max"]][["Lung"]]
leave.protein.out.COVID.list[["Acuity.max"]][["Conventional"]]

# Inspect the distributions
hist(-log10(leave.protein.out.COVID.list[["COVID"]][["Lung"]]$pLRTAB)*sign(leave.protein.out.COVID.list[["COVID"]][["Lung"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.COVID.list[["COVID"]][["Conventional"]]$pLRTAB)*sign(leave.protein.out.COVID.list[["COVID"]][["Conventional"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.COVID.list[["Acuity.max"]][["Lung"]]$pLRTAB)*sign(leave.protein.out.COVID.list[["Acuity.max"]][["Lung"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.COVID.list[["Acuity.max"]][["Conventional"]]$pLRTAB)*sign(leave.protein.out.COVID.list[["Acuity.max"]][["Conventional"]]$LRTstat), breaks = 100)


scatterplot.COVID.leaveout <- vector(mode = "list", length = 4)
names(scatterplot.COVID.leaveout) <- c("Lung - COVID status",
                                        "Conventional - COVID status", 
                                        "Lung - acuity", 
                                        "Conventional - acuity")

coef.index <- c(1, 2, 1, 2)
outcome.index <- c(1, 1, 2, 2)

for(j in 1:length(scatterplot.COVID.leaveout)){
  
  plot.df <- data.frame(
    protein = names(coefs[[names(coefs)[coef.index[j]]]]), 
    `Minimum coefficient` = apply(coefs[[coef.index[j]]], 2, min),
    `Median coefficient` = apply(coefs[[coef.index[j]]], 2, median),
    `Maximum coefficient` = apply(coefs[[coef.index[j]]], 2, max),
    `Signed -log10(p-value) likelihood ratio test` = -log10(leave.protein.out.COVID.list[[outcomes[outcome.index[j]]]][[names(coefs)[coef.index[j]]]]$pLRTAB)*sign(leave.protein.out.COVID.list[[outcomes[outcome.index[j]]]][[names(coefs)[coef.index[j]]]]$LRTstat),
    pval = leave.protein.out.COVID.list[[outcomes[outcome.index[j]]]][[names(coefs)[coef.index[j]]]]$pLRTAB,
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
  
  scatterplot.COVID.leaveout[[j]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = `Signed -log10(p-value) likelihood ratio test`, y = `Median coefficient`)) +
      geom_pointrange(aes(ymin = `Minimum coefficient`, ymax = `Maximum coefficient`, color = Significance), size = 0.7)  +
      scale_color_manual(values = alpha(c("black", "#ff0000ff", "#00be64ff"), 0.5), breaks = levels(plot.df$Significance)) +
      ggtitle("") + # ggtitle("Chronological age adjusted for sex") +
      ggtitle(names(scatterplot.COVID.leaveout)[j]) +
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
  
  attr(scatterplot.COVID.leaveout[[j]], "data") <- plot.df
  
}

plot(scatterplot.COVID.leaveout[[1]])
plot(scatterplot.COVID.leaveout[[2]])
plot(scatterplot.COVID.leaveout[[3]])
plot(scatterplot.COVID.leaveout[[4]])




