### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
df.test.list <- readRDS(file = paste0(rds.dir, "df_test_list.rds"))

### We are only interested in the associations of these organs with these diseases ###

diseases <- c("dementia", "kidney_failure", "cirrhosis_fibrosis", "COPD")

coef.list <- vector(mode = "list", length = 4)
names(coef.list) <- c("Brain", "Kidney", "Liver", "Lung")

coef.list[["Brain"]] <- read.csv(paste0(dir.gen2.models, "Brain", "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
coef.list[["Kidney"]] <- read.csv(paste0(dir.gen2.models, "Kidney", "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
coef.list[["Liver"]] <- read.csv(paste0(dir.gen2.models, "Liver", "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
coef.list[["Lung"]] <- read.csv(paste0(dir.gen2.models, "Lung", "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)



leave.protein.out.list <- vector(mode = "list", length = length(organ.proteins))
names(leave.protein.out.list) <- names(organ.proteins)

for(k in 1:length(leave.protein.out.list)){
  if(names(organ.proteins)[k] != "Bladder"){
  leave.protein.out.list[[k]] <- vector(mode = "list", length = length(general_hazard_outcomes))
  names(leave.protein.out.list[[k]]) <- names(general_hazard_outcomes)
  
  leave.protein.out.list[[k]] <- lapply(leave.protein.out.list[[k]], function(x){
    x <- data.frame(protein = organ.proteins[[k]],
                    `exp(coef)` = NA,
                    z = NA,
                    `Pr(>|z|)` = NA, 
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

### This takes quite a long time! ###

for(k in 1:length(organ.proteins)){
  # if(names(organ.proteins)[k] != "Bladder"){
    if(names(organ.proteins)[k] %in% c("Brain", "Kidney", "Liver", "Lung")){
    coef.tmp <- read.csv(paste0(dir.gen2.models, names(organ.proteins)[k], "_mortality_coefs_GTEx_4x_FC.csv"), header = TRUE, check.names = FALSE)
    
    for(j in 1:length(coef.tmp)){
      
      tmp <- vector(mode = "list", length = 5)
      for(m in 1:5){
        tmp[[m]] <- rowSums(sweep(as.matrix(df.test.list[[m]][, colnames(coef.tmp)[-j], drop = FALSE]), MARGIN=2, unlist(coef.tmp[m, -j]), `*`))
      }
      
      predicted.ages.oof.minus.j <- (do.call("c", tmp))[as.character(olink_bd_annotation_0$eid)]
      
      for(i in 1:length(names(general_hazard_outcomes))){
        
        # Only the cases we need for now:
        if(
          ((names(organ.proteins)[k] == names(coef.list)[1]) & (names(general_hazard_outcomes)[i] == diseases[1])) |
          ((names(organ.proteins)[k] == names(coef.list)[2]) & (names(general_hazard_outcomes)[i] == diseases[2])) |
          ((names(organ.proteins)[k] == names(coef.list)[3]) & (names(general_hazard_outcomes)[i] == diseases[3])) |
          ((names(organ.proteins)[k] == names(coef.list)[4]) & (names(general_hazard_outcomes)[i] == diseases[4])) 
          ){
          
          df <- data.frame(
            eid = olink_bd_annotation_0$eid,
            status_os = olink_bd_annotation_0[, paste0("status_", names(general_hazard_outcomes)[i])],
            time = olink_bd_annotation_0[, paste0("time_", names(general_hazard_outcomes)[i])],
            chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
            biological_age_full = predicted.ages.oof[["gen2"]][["predicted"]][[k]],
            biological_age_noJ = predicted.ages.oof.minus.j,
            sex = olink_bd_annotation_0$p31
          )
          
          # Both need to be fit with x = TRUE, so that the x matrix is returned in component x!
          cox.model_full <- coxph(Surv(time = time, event = status_os, type='right') ~ chronological_age * sex + biological_age_full, data = df, x = TRUE)
          cox.model_noJ <- coxph(Surv(time = time, event = status_os, type='right') ~ chronological_age * sex + biological_age_noJ, data = df, x = TRUE)
          plr.object <- plrtest(object1 = cox.model_full, object2 = cox.model_noJ, nested = FALSE)
          
          leave.protein.out.list[[k]][[i]][(leave.protein.out.list[[k]][[i]]$protein == colnames(coef.tmp)[j]), c("exp(coef)", "z", "Pr(>|z|)", "ll", "ul", "LRTstat", "pLRTA", "pLRTB", "pLRTAB")] <- c(summary(cox.model_noJ)$coef["biological_age_noJ", c("exp(coef)", "z", "Pr(>|z|)")], exp(confint(cox.model_noJ, parm = "biological_age_noJ")), unlist(plr.object[c("LRTstat", "pLRTA", "pLRTB", "pLRTAB")]))
          
        }
        
      }
      
      
    }
  }
}

plr.object
# Variance test 
# H0: Model 1 and Model 2 are indistinguishable 
# H1: Model 1 and Model 2 are distinguishable 
# Fine: p = 0.5
# Non-nested likelihood ratio test 
# H0: Model fits are equally close to true Model 
# H1A: Model 1 fits better than Model 2 
# z =  NaN,   p = NA
# H1B: Model 2 fits better than Model 1 
# z =  NaN,   p = NA
# H1: Model fits not equally close to true Model 
# z =  NaN,   two-sided p = NA

# pLRTA small and LRTstat positive means that the full model fits better than the model in which the specific protein is removed, as one would expect
# pLRTAB is for the two-sided test
# Hence most appropriate: signed -log10(pLRTAB)

# saveRDS(leave.protein.out.list, file = paste0(rds.dir, "leave_protein_out_list.rds"))
leave.protein.out.list <- readRDS(file = paste0(rds.dir, "leave_protein_out_list.rds"))

leave.protein.out.list[["Brain"]][["dementia"]]
leave.protein.out.list[["Kidney"]][["kidney_failure"]]
leave.protein.out.list[["Liver"]][["cirrhosis_fibrosis"]]
leave.protein.out.list[["Lung"]][["COPD"]]

# Inspect the distributions
hist(-log10(leave.protein.out.list[["Brain"]][["dementia"]]$pLRTAB)*sign(leave.protein.out.list[["Brain"]][["dementia"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.list[["Kidney"]][["kidney_failure"]]$pLRTAB)*sign(leave.protein.out.list[["Kidney"]][["kidney_failure"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.list[["Liver"]][["cirrhosis_fibrosis"]]$pLRTAB)*sign(leave.protein.out.list[["Liver"]][["cirrhosis_fibrosis"]]$LRTstat), breaks = 100)
hist(-log10(leave.protein.out.list[["Lung"]][["COPD"]]$pLRTAB)*sign(leave.protein.out.list[["Lung"]][["COPD"]]$LRTstat), breaks = 100)


scatterplot.leaveout <- vector(mode = "list", length = 4)
names(scatterplot.leaveout) <- c("Brain model - dementia",
                                 "Kidney model - kidney failure", 
                                 "Liver model - liver cirrhosis/fibrosis", 
                                 "Lung model - COPD")

for(j in 1:length(scatterplot.leaveout)){
  
  plot.df <- data.frame(
    protein = organ.proteins[[names(coef.list)[j]]], 
    `Minimum coefficient` = apply(coef.list[[j]], 2, min),
    `Median coefficient` = apply(coef.list[[j]], 2, median),
    `Maximum coefficient` = apply(coef.list[[j]], 2, max),
    `Signed -log10(p-value) partial likelihood ratio test` = -log10(leave.protein.out.list[[names(coef.list)[j]]][[diseases[j]]]$pLRTAB)*sign(leave.protein.out.list[[names(coef.list)[j]]][[diseases[j]]]$LRTstat),
    pval = leave.protein.out.list[[names(coef.list)[j]]][[diseases[j]]]$pLRTAB,
    check.names = FALSE
  )
  
  # Exclude proteins that were not included in the model to begin with
  plot.df <- plot.df[!((plot.df[, "Minimum coefficient"] == 0) & (plot.df[, "Maximum coefficient"] == 0)), ]
  
  plot.df$qval <- p.adjust(plot.df$pval, method = "BH")
  
  plot.df$Significance <- "Not significant"
  plot.df$Significance[((plot.df$`Signed -log10(p-value) partial likelihood ratio test` < 0) & (plot.df$qval < 0.05))] <- "Negative effect & q-value < 5%"
  plot.df$Significance[((plot.df$`Signed -log10(p-value) partial likelihood ratio test` > 0) & (plot.df$qval < 0.05))] <- "Positive effect & q-value < 5%"
  plot.df$Significance <- factor(plot.df$Significance, levels = c("Not significant", "Negative effect & q-value < 5%", "Positive effect & q-value < 5%"))
  
  plot.df <- plot.df %>% arrange(-abs(`Signed -log10(p-value) partial likelihood ratio test`))
  plot.df$label <- ""
  plot.df$label[1:min(c(10, nrow(plot.df)))] <- plot.df$protein[1:min(c(10, nrow(plot.df)))]
  
  scatterplot.leaveout[[j]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = `Signed -log10(p-value) partial likelihood ratio test`, y = `Median coefficient`)) +
      geom_pointrange(aes(ymin = `Minimum coefficient`, ymax = `Maximum coefficient`, color = Significance), size = 0.7)  +
      scale_color_manual(values = alpha(c("black", "#ff0000ff", "#00be64ff"), 0.5), breaks = levels(plot.df$Significance)) +
      ggtitle("") + # ggtitle("Chronological age adjusted for sex") +
      ggtitle(names(scatterplot.leaveout)[j]) +
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
  
  attr(scatterplot.leaveout[[j]], "data") <- plot.df
  
}



plot(scatterplot.leaveout[[1]])
plot(scatterplot.leaveout[[2]])
plot(scatterplot.leaveout[[3]])
plot(scatterplot.leaveout[[4]])
