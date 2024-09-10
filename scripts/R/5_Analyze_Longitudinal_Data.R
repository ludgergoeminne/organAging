
### This comes from 1_Prepare_Data.R:
general_hazard_outcomes <- get_general_outcomes()

### This comes from 2_Prepare_Longitudinal_Data.R:
longitudinal.eids <- readRDS(file = paste0(rds.dir, "longitudinal_eids.rds"))

### This comes from 3_Process_GTEx.R:
organ.proteins.longitudinal <- readRDS(paste0(rds.dir, "organ_proteins_longitudinal.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R:
predicted.ages.longitudinal.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_longitudinal_1.rds"))
coefficients.longitudinal <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal.rds"))
predicted.ages.longitudinal.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_longitudinal_oof.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

# Only needed for export to Excel
xlsx.quantiles <- c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99, 0.999, 1)

# Sanity check
all(olink_bd_annotation_list[["instance_0"]]$eid[olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids] == longitudinal.eids)
all(olink_bd_annotation_list[["instance_2"]]$eid[olink_bd_annotation_list[["instance_2"]]$eid %in% longitudinal.eids] == longitudinal.eids)
all(olink_bd_annotation_list[["instance_3"]]$eid[olink_bd_annotation_list[["instance_3"]]$eid %in% longitudinal.eids] == longitudinal.eids)
# Must be all TRUE!

### Prepare slopes.gen2 ###

### Out of fold predictions for the first visit, predictions of the first fold for the two other visits

### Create longitudinal data frames for all clocks ###

longitudinal.df.list <- vector(mode = "list", length = length(organ.proteins.longitudinal))
names(longitudinal.df.list) <- names(organ.proteins.longitudinal)

for(k in 1:length(longitudinal.df.list)){
  
  if(!(names(longitudinal.df.list)[k] %in% c("Bladder", "Thyroid"))){
    # olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids can be used because eids are sorted!!!
    longitudinal.df.list[[k]] <- data.frame(eid = rep(longitudinal.eids, times = 3),
                                            predicted.age = c(predicted.ages.longitudinal.oof$gen2$predicted[[k]][longitudinal.eids], predicted.ages.longitudinal.1$instance_2$gen2$predicted[[k]][longitudinal.eids], predicted.ages.longitudinal.1$instance_3$gen2$predicted[[k]][longitudinal.eids]),
                                            resid.age = c(predicted.ages.longitudinal.oof$gen2$residual[[k]][longitudinal.eids], predicted.ages.longitudinal.1$instance_2$gen2$residual[[k]][longitudinal.eids], predicted.ages.longitudinal.1$instance_3$gen2$residual[[k]][longitudinal.eids]),
                                            age = c(olink_bd_annotation_list[["instance_0"]][olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids, "age_first_visit"], olink_bd_annotation_list[["instance_2"]][olink_bd_annotation_list[["instance_2"]]$eid %in% longitudinal.eids, "age_third_visit"], olink_bd_annotation_list[["instance_3"]][olink_bd_annotation_list[["instance_3"]]$eid %in% longitudinal.eids, "age_fourth_visit"]),
                                            delta_age = c(rep(0, length(longitudinal.eids)), (olink_bd_annotation_list[["instance_2"]][olink_bd_annotation_list[["instance_2"]]$eid %in% longitudinal.eids, "age_third_visit"]-olink_bd_annotation_list[["instance_0"]][olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids, "age_first_visit"]), (olink_bd_annotation_list[["instance_3"]][olink_bd_annotation_list[["instance_3"]]$eid %in% longitudinal.eids, "age_fourth_visit"]-olink_bd_annotation_list[["instance_0"]][olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids, "age_first_visit"])),
                                            sex = rep(factor(ifelse((olink_bd_annotation_list[["instance_0"]][olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men")), times = 3),
                                            visit = factor(rep(c("first", "third", "fourth"), each = length(longitudinal.eids)), levels = c("first", "third", "fourth")))
  }
}


### First the slopes for all mortality-based models ###

slopes.gen2 <- vector(mode = "list", length = length(organ.proteins.longitudinal))
names(slopes.gen2) <- names(organ.proteins.longitudinal)

for(k in 1:length(slopes.gen2)){
  if(!(names(slopes.gen2)[k] %in% c("Bladder", "Thyroid"))){
    
    mixed.model <- lmer(predicted.age ~ delta_age + (1|eid) + (0 + delta_age|eid), data = longitudinal.df.list[[k]])
    # summary(mixed.model)
    # fixef(mixed.model)["delta_age"]
    # ranef(mixed.model)
    slopes.gen2[[k]] <- fixef(mixed.model)["delta_age"]+ranef(mixed.model)$eid[, "delta_age"]
    names(slopes.gen2[[k]]) <- rownames(ranef(mixed.model)$eid)
    
  }
}

### 1. Create the histogram with the longitudinal slopes.gen2 ###

hist.slopes.gen2.longitudinal.gen2 <- vector(mode = "list", length = length(organ.proteins.longitudinal))
names(hist.slopes.gen2.longitudinal.gen2) <- names(organ.proteins.longitudinal)

# To get the range of the slopes for the models of interest
range(do.call("c", slopes.gen2[names(organ.proteins.selected)]))
# -0.04803299  0.18459814
lims <- c(-0.06, 0.21)

for(k in 1:length(slopes.gen2)){
  if(!(names(slopes.gen2)[k] %in% c("Bladder", "Thyroid"))){
    
    data <- ggplot_build(
      ggplot(data = data.frame(slope = slopes.gen2[[k]]), aes(x = slope)) +
        geom_histogram(bins = 30, fill = "#69b3a2", color = "#e9ecef", alpha = 0.9) +
        ggtitle(paste0("Slopes ", gsub("_", "-", paste(tolower(substr(names(slopes.gen2)[k], 1, 1)), substr(names(slopes.gen2)[k], 2, nchar(names(slopes.gen2)[k])), sep="")), " model")) +
        xlab(paste0("DELTA log(mortality hazard) / year"))+
        
        ylab("Counts") +
        theme_bw() +
        scale_x_continuous(limits = lims)+
        
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
    )
    hist.slopes.gen2.longitudinal.gen2[[k]] <- ggplot_gtable(data)
    attr(hist.slopes.gen2.longitudinal.gen2[[k]], "data") <- data$data[[1]]
    
  }
}

### Get numbers ###
sum(slopes.gen2[["Conventional"]] < 0)
# 2

### 2.1. Prepare the heatmap for organ-specific longitudinal aging with sex ###

df.sex.deviations.gen2.longitudinal <- data.frame(
  organ = names(organ.proteins),
  `DELTA log(mortality hazard) / year` = rep(NA, length(organ.proteins)),
  pval = rep(NA, length(organ.proteins)), 
  check.names = FALSE
)

for(k in 1:length(organ.proteins)){
  if(!(names(organ.proteins)[k] %in% c("Bladder", "Thyroid"))){
    
    df <- data.frame(
      eid = olink_bd_annotation_list[["instance_0"]][longitudinal.eids, "eid"],
      chronological_age = olink_bd_annotation_list[["instance_0"]][longitudinal.eids, "age_first_visit"], 
      slope = slopes.gen2[[k]],
      sex = factor(ifelse((olink_bd_annotation_list[["instance_0"]][longitudinal.eids, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men"))
    )
    
    ### Statistics ###
    t.stat <- tryCatch(t.test(slope~sex, data = df, var.equal = TRUE), error = function(e){
      return(structure(
             "Error", 
             class = c("try-error", "character")
         ))
    })
    
    if(class(t.stat)[1] != "try-error"){
      df.sex.deviations.gen2.longitudinal[k, "DELTA log(mortality hazard) / year"] <- t.stat$estimate[2]-t.stat$estimate[1]
      df.sex.deviations.gen2.longitudinal[k, "pval"] <- t.stat$p.value 
    }
    
  }
}

table(olink_bd_annotation_list[["instance_0"]][longitudinal.eids, "p31"])
# Female   Male 
# 531    454 

### 2.1. Make the heatmap for organ-specific longitudinal aging with sex ###

plot.df <- df.sex.deviations.gen2.longitudinal[df.sex.deviations.gen2.longitudinal$organ %in% names(organ.proteins.selected),]
plot.df$qval <- p.adjust(plot.df$pval, method = "hommel")
signif(plot.df$qval, digits = 2)
# Artery       Brain        Immune       Intestine    Kidney       Liver        Lung  Skin Conventional
# 0.41         4e-03        2e-15       0.77           0.16        1e-04        0.77  0.15    1e-17

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

Heatmap_palette <- c("magenta", "white","#4169E1")

range(plot.df$`DELTA log(mortality hazard) / year`)
# -0.0008175698  0.0132242419
lims <- c(-0.015, 0.015)

plot.df$x <- factor(1)
plot.df$organ <- factor(plot.df$organ, levels = rev(names(organ.proteins.selected)))

heatmap.sex.deviations.gen2.longitudinal <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = x, y = organ, label = stars)) + 
    geom_tile(aes(fill = `DELTA log(mortality hazard) / year`)) + 
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + 
    labs(fill="DELTA log(mortality hazard) / year") + 
    geom_text(size = 4.5, nudge_y = -0.13, na.rm = TRUE) +
    ggtitle("") +
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
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.line.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.line.y = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
          axis.ticks.y = element_blank(), # element_line(color="black", linewidth = 0.5),
          strip.text = element_text(face = "bold",size = 12),
          strip.background = element_blank(),
          plot.margin = unit(c(0, 0, 1, 0), # EXCEPTIONAL: increase bottom margin
                             "inches"))
))

attr(heatmap.sex.deviations.gen2.longitudinal, "data") <- plot.df

plot(heatmap.sex.deviations.gen2.longitudinal)

### 3. Plots slopes per category ###

violin.plots.slopes.gen2.category <- vector(mode = "list", length = length(organ.proteins.longitudinal))
names(violin.plots.slopes.gen2.category) <- names(organ.proteins.longitudinal)

for(k in 1:length(slopes.gen2)){
  if(!(names(slopes.gen2)[k] %in% c("Bladder", "Thyroid"))){
    
    plot.df <- data.frame(
      eid = longitudinal.eids,
      sex = factor(ifelse((olink_bd_annotation_list[["instance_0"]][, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men"))[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids)],
      resid.age.first.visit = longitudinal.df.list[[k]][longitudinal.df.list[[k]]$visit == "first", "resid.age"],
      age.category = "others",
      slope = slopes.gen2[[k]]
    )
    
    plot.df$age.category[plot.df$resid.age.first.visit < quantile(plot.df$resid.age.first.visit, 0.05)] <- "5% lowest"
    plot.df$age.category[plot.df$resid.age.first.visit > quantile(plot.df$resid.age.first.visit, 0.95)] <- "5% highest"
    plot.df$age.category <- factor(plot.df$age.category, levels = c("5% lowest", "others", "5% highest"))

    height1 <- (max(plot.df$slope, na.rm = TRUE) - min(plot.df$slope, na.rm = TRUE))*0.8+min(plot.df$slope, na.rm = TRUE)
    height2 <- (max(plot.df$slope, na.rm = TRUE) - min(plot.df$slope, na.rm = TRUE))*0.97+min(plot.df$slope, na.rm = TRUE)
    
    summary <- summary(glht(lm(slope~age.category, data = plot.df), linfct = mcp(age.category =  c("others - `5% lowest` = 0", "`5% highest` - `5% lowest` = 0", "`5% highest` - others = 0"))))
    # summary$test$pvalues

    significance.df <- data.frame(x = c(1.5, 2, 2.5), 
                                  xend = c(1.9, 2.9, 2.9), 
                                  y = c(height1, height2, height1),
                                  pval = summary$test$pvalues,
                                  stars = c("n.s.", "n.s.", "n.s."),
                                  size = c(4.5, 4.5, 4.5),
                                  nudge_y = (max(plot.df$slope, na.rm = TRUE) - min(plot.df$slope, na.rm = TRUE))/13
    )
    
    significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
    significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
    significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
    significance.df$stars[which(significance.df$pval < 0.001)] <- "***"
    
    significance.df[(significance.df$stars != "n.s."), "size"] <- 10
    significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- (max(plot.df$slope, na.rm = TRUE) - min(plot.df$slope, na.rm = TRUE))/30

    violin.plots.slopes.gen2.category[[k]] <- ggplot_gtable(ggplot_build(
      ggplot(plot.df, aes(x = age.category, y = slope)) +
        
        geom_quasirandom(cex = 1.5, aes(x = age.category, color = age.category, alpha = 1)) +
        scale_color_manual(values=c("#70ad00", "#4585b7", "#80521c")) +
        scale_fill_manual(values=c("#70ad00", "#4585b7", "#80521c")) +
        
        geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
        guides(alpha = "none") + # no legend for "alpha"
        
        # geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
        # 
        # geom_segment(aes(x = significance.df$x[1]-0.4, y = significance.df$y[1], xend = significance.df$xend[1], yend = significance.df$y[1])) +
        # geom_segment(aes(x = significance.df$x[2]-0.4, y = significance.df$y[2], xend = significance.df$xend[2], yend = significance.df$y[2])) +
        
        geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black") +
        
        geom_segment(aes(x = significance.df$x[1]-0.4, y = significance.df$y[1], xend = significance.df$xend[1], yend = significance.df$y[1])) +
        geom_segment(aes(x = significance.df$x[2]-0.9, y = significance.df$y[2], xend = significance.df$xend[2], yend = significance.df$y[2])) +
        geom_segment(aes(x = significance.df$x[3]-0.4, y = significance.df$y[3], xend = significance.df$xend[3], yend = significance.df$y[3])) +
        
        ggtitle(paste0(paste(toupper(substr(names(slopes.gen2)[k], 1, 1)), substr(names(slopes.gen2)[k], 2, nchar(names(slopes.gen2)[k])), sep=""), " model")) +
        xlab("Log(mortality hazard) deviation at first visit") +
        ylab("DELTA log(mortality hazard) / year") +
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
    
    tmp <- data.frame(
      `Log(mortality hazard) deviation at first visit` = rep(c("5% lowest", "others", "5% highest"), each = length(xlsx.quantiles)),
      Quantile = rep(paste0(100*xlsx.quantiles, "%"), times = 3),
      check.names = FALSE
    )
    
    tmp[, "DELTA log(mortality hazard) / year"] <- c(
      round(quantile(plot.df[(plot.df$age.category == "5% lowest"), "slope"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$age.category == "others"), "slope"], xlsx.quantiles), 2), 
      round(quantile(plot.df[(plot.df$age.category == "5% highest"), "slope"], xlsx.quantiles), 2)
    )
    
    attr(violin.plots.slopes.gen2.category[[k]], "data") <- tmp
    
  }
}

### Statistics ###
k <- 24
plot.df <- data.frame(
  eid = longitudinal.eids,
  sex = factor(ifelse((olink_bd_annotation_list[["instance_0"]][, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men"))[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.eids)],
  resid.age.first.visit = longitudinal.df.list[[k]][longitudinal.df.list[[k]]$visit == "first", "resid.age"],
  age.category = "others",
  slope = slopes.gen2[[k]]
)

plot.df$age.category[plot.df$resid.age.first.visit < quantile(plot.df$resid.age.first.visit, 0.05)] <- "5% lowest"
plot.df$age.category[plot.df$resid.age.first.visit > quantile(plot.df$resid.age.first.visit, 0.95)] <- "5% highest"
plot.df$age.category <- factor(plot.df$age.category, levels = c("5% lowest", "others", "5% highest"))

summary <- summary(glht(lm(slope~age.category, data = plot.df), linfct = mcp(age.category =  c("others - `5% lowest` = 0", "`5% highest` - `5% lowest` = 0", "`5% highest` - others = 0"))))
summary$test$pvalues
# 1.951810e-05 2.519045e-08 1.802663e-03

### 4. Plots number of ICD-10 events in function of slopes ###

### Plot the correlation between the slopes.gen2 and the number of ICD-10 disease codes ###

all.diseases <- paste0("p", seq(from = 130000, to = 132604, by = 2))
all.diseases <- all.diseases[all.diseases %in% colnames(olink_bd_annotation_list[["instance_0"]])]

tmp <- olink_bd_annotation_list[["instance_0"]]
rownames(tmp) <- olink_bd_annotation_list[["instance_0"]]$eid
n.disease <- rowSums(!is.na(tmp[longitudinal.eids, all.diseases])) # - rowSums(!is.na(tmp[longitudinal.eids, exclude.diseases]))

test <- as.data.frame(do.call("cbind", slopes.gen2))
test$n.disease <- n.disease
plot.df <- pivot_longer(test, cols = 1:(ncol(test)-1), names_to = "organ", values_to = "slope")
plot.df <- plot.df[plot.df$organ == "Conventional", ]
plot.df$sex <- factor(tmp[longitudinal.eids, "p31"], ordered = FALSE)

summary <- summary(lm(n.disease ~ slope, data = plot.df))
resid <- resid(lm(n.disease ~ slope, data = plot.df))
metrics <- paste0("r = ", sprintf("%.2f", round(sign(summary$coefficients["slope", "Estimate"])*sqrt(summary$r.squared), 2)), ", r² = ", sprintf("%.2f", round(summary$r.squared, 2)), ", MAE = ", sprintf("%.2f", round(mean(abs(resid)), 2)))

correlation.slopes.gen2.n.diseases <- ggplot_gtable(ggplot_build(
  ggplot(data = plot.df, aes(x = slope, y = n.disease)) +
    geom_point(size = 1, aes(color = sex, shape = sex, alpha = 0.05)) + # aes(color = Gender, shape = Sequence)
    scale_color_manual(values=c("magenta", "#4169E1")) +
    geom_smooth(method = MASS::rlm , color="turquoise", fill="turquoise", se=TRUE) + # df[df$diet == "CD",]
    ggtitle(paste0(metrics)) +
    xlab("DELTA log(mortality hazard) / year") +
    ylab("Number of ICD-10 events") +
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
      panel.grid.minor = element_blank()
    )
))

plot.df <- as.data.frame(plot.df)

tmp.male <- kMeansEqual(kdat = plot.df[plot.df$sex == "Male", ], 
                        k = floor(nrow(plot.df[plot.df$sex == "Male", ])/10),
                        columns = c("slope", "n.disease"),
                        add.centers = FALSE
)
tmp.female <- kMeansEqual(kdat = plot.df[plot.df$sex == "Female", ], 
                          k = floor(nrow(plot.df[plot.df$sex == "Female", ])/10),
                          columns = c("slope", "n.disease"),
                          add.centers = FALSE
)

tmp.male <- tmp.male %>% group_by(assigned, sex) %>% summarize(
  `N participants` = n(),
  `Avg. DELTA log(mortality hazard) / year` = mean(slope, na.rm = TRUE), 
  `Avg. number of ICD-10 events` = mean(n.disease, na.rm = TRUE)
) 

tmp.female <- tmp.female %>% group_by(assigned, sex) %>% summarize(
  `N participants` = n(),
  `Avg. DELTA log(mortality hazard) / year` = mean(slope, na.rm = TRUE), 
  `Avg. number of ICD-10 events` = mean(n.disease, na.rm = TRUE)
) 

tmp <- rbind(tmp.male, tmp.female) %>% arrange(`Avg. DELTA log(mortality hazard) / year`)
tmp <- tmp[, c("N participants", "sex", "Avg. DELTA log(mortality hazard) / year", "Avg. number of ICD-10 events")]
colnames(tmp)[2] <- "Sex"

attr(correlation.slopes.gen2.n.diseases, "data") <- as.data.frame(tmp)
rm(tmp.male, tmp.female, tmp)

plot(correlation.slopes.gen2.n.diseases)

### Statistics ###
summary
# Call:
#   lm(formula = n.disease ~ slope, data = plot.df)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -12.413  -4.853  -1.302   3.267  37.245 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.6309     0.6701   9.895  < 2e-16 ***
#   slope        56.8030     9.1904   6.181 9.34e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 6.905 on 983 degrees of freedom
# Multiple R-squared:  0.03741,	Adjusted R-squared:  0.03643 
# F-statistic:  38.2 on 1 and 983 DF,  p-value: 9.338e-10

### For Supplementary Fig. S2 ###

### Correlations with chronological age ###

prediction.plots.longitudinal.gen2 <- vector(mode = "list", length = 4)
names(prediction.plots.longitudinal.gen2) <- c("training dataset", "test dataset", "third visit", "fourth visit")

age.visit <- c("first", "first", "third", "fourth")

for(i in 1:4){
  
  if(names(prediction.plots.longitudinal.gen2)[i] == "training dataset"){
    plot.df <- data.frame(
      chronological_age = olink_bd_annotation_list[["instance_0"]]$age_first_visit[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.training.eids)],
      biological_age = predicted.ages.longitudinal.1[["instance_0"]]$gen2$predicted[["Conventional"]][longitudinal.training.eids],
      sex = factor(ifelse((olink_bd_annotation_list[["instance_0"]][, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men"))[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.training.eids)]
      )
  } else if(names(prediction.plots.longitudinal.gen2)[i] == "test dataset"){
    plot.df <- data.frame(
      chronological_age = olink_bd_annotation_list[["instance_0"]]$age_first_visit[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.test.eids)],
      biological_age = predicted.ages.longitudinal.1[["instance_0"]]$gen2$predicted[["Conventional"]][longitudinal.test.eids],
      sex = factor(ifelse((olink_bd_annotation_list[["instance_0"]][, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men"))[(olink_bd_annotation_list[["instance_0"]]$eid %in% longitudinal.test.eids)]
      )
  } else if(names(prediction.plots.longitudinal.gen2)[i] == "third visit"){
    plot.df <- data.frame(
      chronological_age = olink_bd_annotation_list[["instance_2"]]$age_third_visit,
      biological_age = predicted.ages.longitudinal.1[["instance_2"]]$gen2$predicted[["Conventional"]],
      sex = factor(ifelse((olink_bd_annotation_list[["instance_2"]][, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men"))
      )
  } else if(names(prediction.plots.longitudinal.gen2)[i] == "fourth visit"){
    plot.df <- data.frame(
      chronological_age = olink_bd_annotation_list[["instance_3"]]$age_fourth_visit,
      biological_age = predicted.ages.longitudinal.1[["instance_3"]]$gen2$predicted[["Conventional"]],
      sex = factor(ifelse((olink_bd_annotation_list[["instance_3"]][, "p31"]) == "Female", "Women", "Men"), levels = c("Women", "Men"))
      )
  }
  plot.df <- plot.df[!is.na(plot.df$chronological_age),]
  
  summary <- summary(lm(biological_age ~ chronological_age, data = plot.df))
  resid <- resid(lm(biological_age ~ chronological_age, data = plot.df))
  metrics <- paste0("r = ", round(sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared), 2), ", r² = ", round(summary$r.squared, 2), ", MAE = ", round(mean(abs(resid)), 2))
  
  prediction.plots.longitudinal.gen2[[i]] <- ggplot_gtable(ggplot_build(
    ggplot(data = plot.df, aes(x = chronological_age, y = biological_age)) +
      geom_point(size = 1, aes(color = sex, shape = sex, alpha = 0.05)) + 
      scale_color_manual(values=c("magenta", "#4169E1")) +
      geom_smooth(method = MASS::rlm , color="turquoise", fill="turquoise", se=TRUE) + 
      ggtitle(paste0("Mortality-based conventional model: ", names(prediction.plots.longitudinal.gen2)[i], ": ", sum(coefficients.longitudinal[["gen2"]][["Conventional"]] != 0), " / ", length(organ.proteins.longitudinal[["Conventional"]]), " proteins", "\n", metrics)) +
      xlab(paste0("Chronological age at ", age.visit[i], " visit")) +
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
  
  tmp.male <- kMeansEqual(kdat = plot.df[plot.df$sex == "Men", ], 
                          k = floor(nrow(plot.df[plot.df$sex == "Men", ])/10),
                          columns = c("chronological_age", "biological_age"),
                          add.centers = FALSE
  )
  tmp.female <- kMeansEqual(kdat = plot.df[plot.df$sex == "Women", ], 
                            k = floor(nrow(plot.df[plot.df$sex == "Women", ])/10),
                            columns = c("chronological_age", "biological_age"),
                            add.centers = FALSE
  )
  
  tmp.male <- tmp.male %>% group_by(assigned, sex) %>% summarize(
    `N participants` = n(),
    `Avg. chronological age` = mean(chronological_age, na.rm = TRUE), 
    `Avg. predicted rel. log(mortality hazard)` = mean(biological_age, na.rm = TRUE)
  ) 
  
  tmp.female <- tmp.female %>% group_by(assigned, sex) %>% summarize(
    `N participants` = n(),
    `Avg. chronological age` = mean(chronological_age, na.rm = TRUE), 
    `Avg. predicted rel. log(mortality hazard)` = mean(biological_age, na.rm = TRUE)
  ) 
  
  tmp <- rbind(tmp.male, tmp.female) %>% arrange(`Avg. chronological age`)
  tmp <- tmp[, c("N participants", "sex", "Avg. chronological age", "Avg. predicted rel. log(mortality hazard)")]
  colnames(tmp)[2] <- "Sex"
  
  attr(prediction.plots.longitudinal.gen2[[i]], "data") <- as.data.frame(tmp)
  rm(tmp.male, tmp.female, tmp)
  
}

### Get the number of participants ###

dim(olink_bd_annotation_list[["instance_0"]])
# 52341  5692

length(longitudinal.training.eids)
# 41873

length(longitudinal.test.eids)
# 10468

dim(olink_bd_annotation_list[["instance_2"]])
# 1161 5692

dim(olink_bd_annotation_list[["instance_3"]])
# 1113 5692

length(longitudinal.eids)
# 985

### For Supplementary Fig. S4 ###

### Association between Log(mortality hazard) deviation at first visit and rate of aging ###

plots.longitudinal.change.vs.first.visit <- vector(mode = "list", length = length(organ.proteins.longitudinal))
names(plots.longitudinal.change.vs.first.visit) <- names(organ.proteins.longitudinal)

for(k in 1:length(slopes.gen2)){
  if(!(names(slopes.gen2)[k] %in% c("Bladder", "Thyroid"))){
    
    plot.df <- data.frame(
      eid = longitudinal.eids,
      sex = longitudinal.df.list[[k]][longitudinal.df.list[[k]]$visit == "first", "sex"],
      resid.age.first.visit = longitudinal.df.list[[k]][longitudinal.df.list[[k]]$visit == "first", "resid.age"],
      slope = slopes.gen2[[k]]
    )
    
    model <- mgcv::gam(slope~s(resid.age.first.visit), data = plot.df)
    prediction <- predict(model, newdata = plot.df, type = "link", se.fit = TRUE) # , unconditional = TRUE
    
    prediction.df <- data.frame(
      resid.age.first.visit = plot.df$resid.age.first.visit,
      est = prediction$fit,
      upper = prediction$fit + (qnorm(0.975) * prediction$se.fit),
      lower = prediction$fit - (qnorm(0.975) * prediction$se.fit)
    )

    plots.longitudinal.change.vs.first.visit[[k]] <- ggplot_gtable(ggplot_build(
      ggplot(data = plot.df, aes(x = resid.age.first.visit, y = slope)) +
        geom_point(size = 1, aes(color = sex, shape = sex, alpha = 0.05)) + # aes(color = Gender, shape = Sequence)
        scale_color_manual(values=c("magenta", "#4169E1")) +
        geom_ribbon(aes(x = resid.age.first.visit, y = est, ymin = lower, ymax = upper), alpha = 0.3, data = prediction.df, inherit.aes = FALSE, fill = "turquoise") +  
        geom_line(aes(x = resid.age.first.visit, y = est),  data = prediction.df, inherit.aes = FALSE, color = "turquoise") +
        ggtitle(paste0(paste(toupper(substr(names(slopes.gen2)[k], 1, 1)), substr(names(slopes.gen2)[k], 2, nchar(names(slopes.gen2)[k])), sep=""), " model")) +
        xlab("Log(mortality hazard) deviation at first visit") +
        ylab("DELTA log(mortality hazard) / year") +
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
    
    tmp.male <- kMeansEqual(kdat = plot.df[plot.df$sex == "Men", ], 
                            k = floor(nrow(plot.df[plot.df$sex == "Men", ])/10),
                            columns = c("resid.age.first.visit", "slope"),
                            add.centers = FALSE
    )
    tmp.female <- kMeansEqual(kdat = plot.df[plot.df$sex == "Women", ], 
                              k = floor(nrow(plot.df[plot.df$sex == "Women", ])/10),
                              columns = c("resid.age.first.visit", "slope"),
                              add.centers = FALSE
    )
    
    tmp.male <- tmp.male %>% group_by(assigned, sex) %>% summarize(
      `N participants` = n(),
      `Avg. log(mortality hazard) deviation at first visit` = mean(resid.age.first.visit, na.rm = TRUE), 
      `Avg. DELTA log(mortality hazard) / year` = mean(slope, na.rm = TRUE)
    ) 
    
    tmp.female <- tmp.female %>% group_by(assigned, sex) %>% summarize(
      `N participants` = n(),
      `Avg. log(mortality hazard) deviation at first visit` = mean(resid.age.first.visit, na.rm = TRUE), 
      `Avg. DELTA log(mortality hazard) / year` = mean(slope, na.rm = TRUE)
    ) 
    
    tmp <- rbind(tmp.male, tmp.female) %>% arrange(`Avg. log(mortality hazard) deviation at first visit`)
    tmp <- tmp[, c("N participants", "sex", "Avg. log(mortality hazard) deviation at first visit", "Avg. DELTA log(mortality hazard) / year")]
    colnames(tmp)[2] <- "Sex"
    
    attr(plots.longitudinal.change.vs.first.visit[[k]], "data") <- as.data.frame(tmp)
    rm(tmp.male, tmp.female, tmp)
    
  }
}


