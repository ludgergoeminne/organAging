
### Fit cox-proportional Hazards models ###

# Note: people who already had the disease are excluded, as well as people for whom the data of the disease was not specified

### This comes from 4_Calculate_Predicted_Resid_Ages.R
training.eids <- readRDS(file = paste0(rds.dir, "training_eids.rds"))
test.eids <- readRDS(file = paste0(rds.dir, "test_eids.rds"))
predicted.ages.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_1.rds"))
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

# See also 1_Prepare_Data.R
general_hazard_outcomes <- get_general_outcomes()
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

HR.list <- vector(mode = "list", length = 2)
names(HR.list) <- c("gen1", "gen2")

HR.list <- lapply(HR.list, function(x){
  x <- vector(mode = "list", length = length(organ.proteins))
  names(x) <- names(organ.proteins)
  
  x <- lapply(x, function(y){
    return(data.frame(outcome = rep(names(general_hazard_outcomes), each = 3),
                      set = factor(rep(c("Training set", "Test set", "Out of fold"), times = length(general_hazard_outcomes)), levels = c("Training set", "Test set", "Out of fold")),
                      `exp(coef)` = NA,
                      coef = NA,
                      z = NA,
                      `Pr(>|z|)` = NA, 
                      ll = NA, 
                      ul = NA, 
                      check.names = FALSE))
  })
  
  return(x)
} )

### This takes a bit of time ###
for(g in 1:length(HR.list)){
  
  for(k in 1:length(HR.list[[g]])){
    
    if(names(HR.list[[g]])[k] != "Bladder"){
      
      for(i in 1:length(names(general_hazard_outcomes))){
        
        df <- data.frame(
          eid = as.character(olink_bd_annotation_0$eid),
          status_os = olink_bd_annotation_0[, paste0("status_", names(general_hazard_outcomes)[i])],
          time = olink_bd_annotation_0[, paste0("time_", names(general_hazard_outcomes)[i])],
          chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
          residual_age_1 = predicted.ages.1[[g]]$residual[[k]],
          residual_age_oof = predicted.ages.oof[[g]]$residual[[k]],
          sex = olink_bd_annotation_0$p31
        )
        
        tmp <- c(df[training.eids, "residual_age_1"]/sd(df[training.eids, "residual_age_1"]), df[test.eids, "residual_age_1"]/sd(df[test.eids, "residual_age_1"]))
        names(tmp) <- c(training.eids, test.eids)
        
        df$norm_residual_age_1 <- tmp[df$eid]
        df$norm_residual_age_oof <- df$residual_age_oof/sd(df$residual_age_oof)
        
        # df$time[which(df$time < 0)] <- 0 # The ones that already had the disease before coming to the center, we put at 0
        # df <- df[!is.na(df$time),]
        # Forward-looking: only include people who have the disease in the future; this also excludes NA values in time
        df <- df[which(df$time >= 0),]
        
        for(j in 1:3){
          
          ### Training set
          if(j == 1){
            cox.model <- coxph(Surv(time = time, event = status_os, type='right') ~ chronological_age * sex + norm_residual_age_1, data = df[(df$eid %in% training.eids), ])
            # summary(cox.model)
            HR.list[[g]][[k]][(HR.list[[g]][[k]]$set == "Training set" & HR.list[[g]][[k]]$outcome == names(general_hazard_outcomes)[i]), c("exp(coef)", "coef", "z", "Pr(>|z|)", "ll", "ul")] <- c(summary(cox.model)$coef["norm_residual_age_1", c("exp(coef)", "coef", "z", "Pr(>|z|)")], exp(confint(cox.model, parm = "norm_residual_age_1")))
          }
          
          ### Test set
          if(j == 2){
            cox.model <- coxph(Surv(time = time, event = status_os, type='right') ~ chronological_age * sex + norm_residual_age_1, data = df[(df$eid %in% test.eids), ])
            # summary(cox.model)
            HR.list[[g]][[k]][(HR.list[[g]][[k]]$set == "Test set" & HR.list[[g]][[k]]$outcome == names(general_hazard_outcomes)[i]), c("exp(coef)", "coef", "z", "Pr(>|z|)", "ll", "ul")] <- c(summary(cox.model)$coef["norm_residual_age_1", c("exp(coef)", "coef", "z", "Pr(>|z|)")], exp(confint(cox.model, parm = "norm_residual_age_1")))
          }
          
          ### Out-of-fold predictions
          if(j == 3){
            cox.model <- coxph(Surv(time = time, event = status_os, type='right') ~ chronological_age * sex + norm_residual_age_oof, data = df)
            # summary(cox.model)
            HR.list[[g]][[k]][(HR.list[[g]][[k]]$set == "Out of fold" & HR.list[[g]][[k]]$outcome == names(general_hazard_outcomes)[i]), c("exp(coef)", "coef", "z", "Pr(>|z|)", "ll", "ul")] <- c(summary(cox.model)$coef["norm_residual_age_oof", c("exp(coef)", "coef", "z", "Pr(>|z|)")], exp(confint(cox.model, parm = "norm_residual_age_oof")))
          }
        }
      }
    }
    
  }
  
}

# Add q-value per group:
for(g in 1:length(HR.list)){
  HR.list[[g]] <- lapply(HR.list[[g]], function(x){
    x <- x %>% group_by(set) %>% mutate(qval = p.adjust(`Pr(>|z|)`, method = "BH"))
    return(x)
  })
}

# saveRDS(HR.list, file = paste0(rds.dir, "HR_list.rds"))
HR.list <- readRDS(file = paste0(rds.dir, "HR_list.rds"))

### Make Hazard ratio plot for gen1 conventional ###

plot.df <- HR.list$gen1$Conventional

plot.df$outcome <- gsub("_", " ", plot.df$outcome)
plot.df$outcome <- paste(toupper(substr(plot.df$outcome, 1, 1)), substr(plot.df$outcome, 2, nchar(plot.df$outcome)), sep="")
plot.df$outcome[plot.df$outcome == "Cirrhosis fibrosis"] <- "Liver cirrhosis/fibrosis"

### Biggest HR in out of fold set first ###
plot.df <- plot.df %>% arrange(`exp(coef)`) %>% arrange(-as.numeric(set))
plot.df$outcome <- factor(plot.df$outcome, levels = unique(plot.df$outcome))

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

lims <- c(0.9, max(plot.df$ul)*1.1)

HR.plot.gen1 <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = `exp(coef)`, y = outcome, label = stars, group = set)) + # aes: , color = set
    geom_pointrange(aes(xmin=ll, xmax=ul), color = "#00BFC4")  +
    facet_wrap(~set) +
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    # scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    guides(color = "none") + # no legend for color
    # labs(fill="NES") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    geom_text(size = 4.5, nudge_y = 0.13, na.rm = TRUE, colour = "black") +
    geom_vline(xintercept = 1, colour = "darkgrey") +
    xlab("Hazard ratio per standard deviation increase") +
    ylab("") +
    xlim(lims) +
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

attr(HR.plot.gen1, "data") <- plot.df

plot(HR.plot.gen1)

plot.df[plot.df$outcome == "Mortality", ]
# Out of fold     HR    1.37

### Make Hazard ratio plot for gen2 conventional ###

# res.df.gen2 <- HR.list$gen2
# for(k in 1:length(res.df.gen2)){
#   res.df.gen2[[k]]$organ <- names(res.df.gen2)[k]
# }
# res.df.gen2 <- do.call("rbind", res.df.gen2)
# 
# plot.df <- res.df.gen2
# 
# plot.df <- res.df.gen2[res.df.gen2$organ == "Conventional",]

plot.df <- HR.list$gen2$Conventional

plot.df$outcome <- gsub("_", " ", plot.df$outcome)
plot.df$outcome <- paste(toupper(substr(plot.df$outcome, 1, 1)), substr(plot.df$outcome, 2, nchar(plot.df$outcome)), sep="")
plot.df$outcome[plot.df$outcome == "Cirrhosis fibrosis"] <- "Liver cirrhosis/fibrosis"

### Biggest HR in out of fold set first ###
plot.df <- plot.df %>% arrange(`exp(coef)`) %>% arrange(-as.numeric(set))
plot.df$outcome <- factor(plot.df$outcome, levels = unique(plot.df$outcome))

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

lims <- c(0.9, max(plot.df$ul)*1.1)

HR.plot.gen2 <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = `exp(coef)`, y = outcome, label = stars, group = set, color = set)) +  # aes: , color = set
    geom_pointrange(aes(xmin=ll, xmax=ul), color = "#F8766D")  +
    facet_wrap(~set) +
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    # scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    guides(color = "none") + # no legend for color
    # labs(fill="NES") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    geom_text(size = 4.5, nudge_y = 0.13, na.rm = TRUE, colour = "black") +
    geom_vline(xintercept = 1, colour = "darkgrey") +
    xlab("Hazard ratio per standard deviation increase") +
    ylab("") +
    xlim(lims) +
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

attr(HR.plot.gen2, "data") <- plot.df

plot(HR.plot.gen2)

plot.df[plot.df$outcome == "Mortality", ]
# Out of fold     HR    2.13

### 3. Make combined Hazard ratio plots ###

HR.plots.organs.combined <- vector(mode = "list", length = length(organ.proteins))
names(HR.plots.organs.combined) <- names(organ.proteins)

organs.to.plot <- c("Liver", "Kidney", "Lung")
HR.plots.organs.combined <- vector(mode = "list", length = length(organs.to.plot))
names(HR.plots.organs.combined) <- organs.to.plot

titles <- organs.to.plot

for(k in 1:length(HR.plots.organs.combined)){
  if(!(names(organ.proteins)[k] %in% c("Bladder"))){
  
  tmp1 <- HR.list[["gen1"]][[names(HR.plots.organs.combined)[k]]]
  tmp1$clock.type <- "1st-generation model"
  tmp2 <- HR.list[["gen2"]][[names(HR.plots.organs.combined)[k]]]
  tmp2$clock.type <- "mortality model"
  
  plot.df <- rbind(tmp1, tmp2)
  plot.df$clock.type <- factor(plot.df$clock.type, levels = rev(c("1st-generation model", "mortality model")))
  
  plot.df$outcome <- gsub("_", " ", plot.df$outcome)
  plot.df$outcome <- paste(toupper(substr(plot.df$outcome, 1, 1)), substr(plot.df$outcome, 2, nchar(plot.df$outcome)), sep="")
  plot.df$outcome[plot.df$outcome == "Cirrhosis fibrosis"] <- "Liver cirrhosis/fibrosis"
  
  ### Biggest HR in out of fold set first ###
  plot.df <- plot.df %>% arrange(-`exp(coef)`) %>% arrange(-as.numeric(set))
  plot.df$outcome <- factor(plot.df$outcome, levels = rev(unique(plot.df$outcome)))
  
  plot.df$stars <- ""
  plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
  plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
  plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
  plot.df$stars[which(plot.df$qval < 0.001)] <- "***"
  
  # lims <- c(0.9, max(plot.df$ul)*1.1)
  lims <- c(0.8, 3.5)
  
  plot.df$y.pos <- as.numeric(plot.df$outcome)
  plot.df[plot.df$clock.type == "1st-generation model", "y.pos"] <- plot.df[plot.df$clock.type == "1st-generation model", "y.pos"]+0.125
  plot.df[plot.df$clock.type == "mortality model", "y.pos"] <- plot.df[plot.df$clock.type == "mortality model", "y.pos"]-0.125
  
  HR.plots.organs.combined[[k]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = `exp(coef)`, y = y.pos, label = stars, group = set, color = clock.type)) + 
      geom_pointrange(aes(xmin=ll, xmax=ul))  +
      facet_grid(~set, scales = "free_y", space="free") +
      scale_y_continuous(breaks = as.numeric(unique(plot.df$outcome)), labels = unique(plot.df$outcome)) +
      guides(color = "none") + # no legend for color
      geom_text(size = 4.5, nudge_y = 0.13, na.rm = TRUE, colour = "black") +
      geom_vline(xintercept = 1, colour = "darkgrey") +
      xlab("Hazard ratio per standard deviation increase") +
      ylab("") +
      xlim(lims) +
      ggtitle(paste0(titles[k], " model")) +
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(), 
            # plot.title = element_text(face = "bold"), # element_blank(),
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
  
  attr(HR.plots.organs.combined[[k]], "data") <- plot.df
  
  }
}

plot(HR.plots.organs.combined[["Liver"]])
plot(HR.plots.organs.combined[["Kidney"]])
plot(HR.plots.organs.combined[["Lung"]])

### Get the exact values to mention in the paper ###

HR.list[["gen1"]][["Liver"]]
HR.list[["gen2"]][["Liver"]]

HR.list[["gen1"]][["Kidney"]]
HR.list[["gen2"]][["Kidney"]]

HR.list[["gen1"]][["Lung"]]
HR.list[["gen2"]][["Lung"]]




