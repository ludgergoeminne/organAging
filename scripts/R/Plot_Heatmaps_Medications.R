
### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

medications <- openxlsx::read.xlsx(paste0(input.coding.dir, "uk_biobank_atc_codes.xlsx"))

medications.lvls <- vector(mode = "list", length = 5)
names(medications.lvls) <- paste0("lvl", 1:5)

for(m in 1:length(medications.lvls)){
  medications.lvls[[m]] <- unique(medications[, paste0("level_", m, "_name")])
}

#### To shorten the run time ####
medicine.indices <- vector(mode = "list", length = 5)
names(medicine.indices) <- paste0("lvl", 1:5)

for(m in 1:length(medicine.indices)){
  medicine.indices[[m]] <- vector(mode = "list", length = 4)
  names(medicine.indices[[m]]) <- c("full", "visit2", "visit3", "visit4")
  medicine.indices[[m]] <- lapply(medicine.indices[[m]], function(y){
    y <- vector(mode = "list", length = length(medications.lvls[[m]]))
    names(y) <- medications.lvls[[m]]
    return(y)
  })
}

# This takes a bit of time to run!
# Each iteration takes longer because the list is longer.
for(m in 1:length(medicine.indices)){
  print(paste0("We are now in iteration ", m, " out of ", length(medicine.indices), "."))
  for(i in 1:length(medicine.indices[[m]][[1]])){
    uk_biobank_codings <- medications[medications[, paste0("level_", m, "_name")] == medications.lvls[[m]][i], "uk_biobank_coding"]
    
    for(j in 1:length(medicine.indices[[m]])){
      if(j == 1){
        tmp <- unname(which(apply(olink_bd_annotation_0[, paste0("p20003_i0_a", 1:47)], 1, function(x){return(any(uk_biobank_codings %in% x))})))
        medicine.indices[[m]][[j]][[i]] <- tmp
        names(medicine.indices[[m]][[j]][[i]]) <- olink_bd_annotation_0[tmp, "eid"]
      } else if(j == 2){
        tmp <- unname(which(apply(olink_bd_annotation_0[, paste0("p20003_i1_a", 1:47)], 1, function(x){return(any(uk_biobank_codings %in% x))})))
        medicine.indices[[m]][[j]][[i]] <- tmp
        names(medicine.indices[[m]][[j]][[i]]) <- olink_bd_annotation_0[tmp, "eid"]
      } else if(j == 3){
        tmp <- unname(which(apply(olink_bd_annotation_0[, paste0("p20003_i2_a", 1:47)], 1, function(x){return(any(uk_biobank_codings %in% x))})))
        medicine.indices[[m]][[j]][[i]] <- tmp
        names(medicine.indices[[m]][[j]][[i]]) <- olink_bd_annotation_0[tmp, "eid"]
      } else if(j == 4){
        tmp <- unname(which(apply(olink_bd_annotation_0[, paste0("p20003_i3_a", 1:47)], 1, function(x){return(any(uk_biobank_codings %in% x))})))
        medicine.indices[[m]][[j]][[i]] <- tmp
        names(medicine.indices[[m]][[j]][[i]]) <- olink_bd_annotation_0[tmp, "eid"]
      }
    }
  }
}

# saveRDS(medicine.indices, file = paste0(rds.dir, "medicine_indices.rds"))
medicine.indices <- readRDS(file = paste0(rds.dir, "medicine_indices.rds"))

########

which(unlist(lapply(medicine.indices[[5]]$full, function(x){length(x) > 1}))) %>% length
# 563 medications with more than 1 person

# We only analyze medicines where we have more than one person taking them:
prefiltered.meds <- lapply(medicine.indices, function(x){
  x <- names(which(unlist(lapply(x$full, function(x){length(x) > 1}))))
  return(x)
})

res.list.meds <- vector(mode = "list", length = length(medications.lvls))
names(res.list.meds) <- names(medications.lvls)

for(m in 1:length(res.list.meds)){
  res.list.meds[[m]] <- vector(mode = "list", length = length(organ.proteins))
  names(res.list.meds[[m]]) <- names(organ.proteins)
}

for(m in 1:length(res.list.meds)){
  res.list.meds[[m]] <- lapply(res.list.meds[[m]], function(x){
    return(data.frame(med = medications.lvls[[m]],
                      Estimate = NA,
                      `Std. Error` = NA,
                      `t value` = NA, 
                      df = NA,
                      `Pr(>|t|)` = NA, 
                      ll = NA, 
                      ul = NA, 
                      r = NA, 
                      check.names = FALSE))
  })
}

# This also takes quite some time to run!
for(m in 1:length(res.list.meds)){
  print(paste0("***We are now in iteration ", m, " out of ", length(res.list.meds), ".***"))
  for(k in 1:length(res.list.meds[[m]])){
    print(paste0("We are now in organ ", k, " out of ", length(res.list.meds[[m]]), " (iteration ", m, " out of ", length(res.list.meds), ")."))
    # if(names(organ.proteins)[k] != "Bladder"){
    if(names(organ.proteins)[k] %in% names(organ.proteins.selected)){
      
      df <- data.frame(
        eid = olink_bd_annotation_0$eid,
        med = FALSE, # apply(olink_bd_annotation_0[, paste0("p20003_i0_a", 1:47)], 1, function(x){return(any(uk_biobank_codings %in% x))}),
        chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
        biological_age = predicted.ages.oof$gen2$predicted[[k]],
        sex = factor(olink_bd_annotation_0$p31, levels = c("Female", "Male"), ordered = FALSE), 
        deprivation = olink_bd_annotation_0$p22189, # Townsend deprivation index
        centre = factor(olink_bd_annotation_0$p54_i0),
        IPAQ_activity = factor(olink_bd_annotation_0$p22032_i0, level = c("low", "moderate", "high"), ordered = FALSE),
        smoking_status = factor(olink_bd_annotation_0$p20116_i0, levels = c("Prefer not to answer", "Never", "Previous", "Current"), ordered = FALSE)
      )
      
      for(i in 1:length(medications.lvls[[m]])){
        # We only analyze medicines where we have more than one person taking them:
        if(medications.lvls[[m]][i] %in% prefiltered.meds[[m]]){

          df$med <- FALSE
          df$med[medicine.indices[[m]][["full"]][[i]]] <- TRUE

          model <- lm(biological_age ~ chronological_age * sex + med + deprivation + centre + IPAQ_activity + smoking_status, data = df)
          summary.vals <- tryCatch(summary(model)$coef["medTRUE", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
          res.list.meds[[m]][[k]][(res.list.meds[[m]][[k]]$med == medications.lvls[[m]][i]), c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model, parm = "medTRUE"), sign(model$coef["medTRUE"])*sqrt(summary(model)$r.squared))
          res.list.meds[[m]][[k]][(res.list.meds[[m]][[k]]$med == medications.lvls[[m]][i]), c("df")] <- model$df.residual
        }
      }
      
      gc()
      
      # Add q-value per group:
      res.list.meds[[m]][[k]] <- res.list.meds[[m]][[k]] %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))
      
    }
  }
}

# saveRDS(res.list.meds, file = paste0(rds.dir, "res_list_meds.rds"))
res.list.meds <- readRDS(file = paste0(rds.dir, "res_list_meds.rds"))

res.df.list.meds <- lapply(res.list.meds, function(x){
  # x[["Bladder"]] <- NULL
  x[which(!(names(x) %in% names(organ.proteins.selected)))] <- NULL
  for(k in 1:length(x)){
    x[[k]]$organ <- names(x)[k]
  }
  x <- do.call("rbind", x)
  x$zval <- limma::zscoreT(x = x$`t value`, df = x$df)
  return(x)
})

# res.list.meds[["Bladder"]] <- NULL
# res.df.meds.gen2 <- do.call("rbind", res.list.meds)

# saveRDS(res.df.list.meds, file = paste0(rds.dir, "res_df_list_meds.rds"))
res.df.list.meds <- readRDS(file = paste0(rds.dir, "res_df_list_meds.rds"))

### Overall most significant ones; also seem mostly most significant on the conventional model ###
head((res.df.list.meds[[5]] %>% arrange(-abs(zval))), n = 20)

head((res.df.list.meds[[5]] %>% arrange(zval)), n = 100)
# acarbose and acebutolol seem to make the heart younger
res.df.list.meds[[5]][which((res.df.list.meds[[5]]$qval) < 0.05 & (res.df.list.meds[[5]]$Estimate < 0)),]

### Acebutolol ###

# Beta blockers cause the heart to beat more slowly and with less force. This lowers blood pressure.
# acebutolol: No HDL decrease has been observed. In this regard, it is unlike many other beta-blockers which have this unfavourable property. 
# https://en.wikipedia.org/wiki/Acebutolol

# Q: how do hypertension patients on these drugs compare to hypertension patients on other drugs in terms of sex, age, deprivation index, diet?

### Acarbose ###
# acarbose is significantly more effective in patients eating a relatively high carbohydrate Eastern diet
# Acarbose is cheap and popular in China, but not in the U.S. One physician explains that use in the U.S. is limited because it is not potent enough to justify the side effects of diarrhea and flatulence.
# https://en.wikipedia.org/wiki/Acarbose

# Q: how do type II diabetes patients on these drugs compare to type II diabetes patients on other drugs in terms of sex, age, deprivation index, diet?

### Plot level 5: the only 2 that are significantly down, as well as the 10 most significantly up in conventional
res.df.list.meds[[5]][which((res.df.list.meds[[5]]$qval) < 0.05 & (res.df.list.meds[[5]]$Estimate < 0)),]

head((res.df.list.meds[[5]][(res.df.list.meds[[5]]$organ == "Conventional"), ] %>% arrange(-abs(zval))), n = 20)
head((res.df.list.meds[[5]][(res.df.list.meds[[5]]$organ == "Conventional"), ] %>% arrange(zval)), n = 10)
head((res.df.list.meds[[5]][(res.df.list.meds[[5]]$organ == "Conventional"), ] %>% arrange(-zval)), n = 10)

length(medicine.indices[["lvl5"]][["full"]][["glucosamine"]])
length(medicine.indices[["lvl5"]][["full"]][["xylometazoline, combinations"]])

top.medications <- c((res.df.list.meds[[5]][(res.df.list.meds[[5]]$organ == "Conventional"), ] %>% arrange(-zval))$med[1:10],
                     rev((res.df.list.meds[[5]][(res.df.list.meds[[5]]$organ == "Conventional"), ] %>% arrange(zval))$med[1:10]))

# saveRDS(top.medications, file = paste0(rds.dir, "top_medications.rds"))
top.medications <- readRDS(file = paste0(rds.dir, "top_medications.rds"))

### Heatmap top medications lvl 5 ###

plot.df <- res.df.list.meds[[5]]
plot.df <- plot.df[plot.df$organ %in% names(organ.proteins.selected),]

# Add q-value per group:
plot.df <- plot.df %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))

# This should come after calculating the q-values!
plot.df <- plot.df[(plot.df$med %in% top.medications),]

# Arrange data frame by top medications (to turn into a factor later)
plot.df <- plot.df %>% arrange(match(med, top.medications))
all(unique(plot.df$med) == top.medications)
# TRUE
plot.df$med <- paste(toupper(substr(plot.df$med, 1, 1)), substr(plot.df$med, 2, nchar(plot.df$med)), sep="")
plot.df$med <- gsub(" and ", " & ", plot.df$med)

plot.df$stars <- ""
plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
plot.df$stars[which(plot.df$qval < 0.001)] <- "***"

lims <- c(-1.5, 1.5)

range(plot.df$Estimate, na.rm = TRUE)
# -1.103937  1.256279

c(quantile(plot.df$Estimate, 0.05, na.rm = TRUE), quantile(plot.df$Estimate, 0.95, na.rm = TRUE))
# -0.5873923  0.8118928

# Remove Bladder
plot.df <- plot.df[plot.df$organ != "Bladder",]

# Only here change to factor!

plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))

plot.df$med = factor(plot.df$med, levels = rev(unique(plot.df$med)))

Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")

heatmap.top.medications.gen2 <- ggplot_gtable(ggplot_build(
  ggplot(plot.df, aes(x = organ, y = med, label = stars)) + 
    geom_tile(aes(fill = Estimate)) + # size = minuslog10p , color = p.adjust
    # scale_shape_manual(values=c(22, 24, 25)) +
    # facet_grid(pathway~diet, scales = "free_y", space="free") +
    scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
    # theme_bw(base_size = 14) +
    # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
    # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
    guides(color = "none") + # no legend for color
    labs(fill="Effect size") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
    geom_text(size = 4.5, nudge_y = -0.13, na.rm = TRUE) +
    xlab("") +
    ylab("") +
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

attr(heatmap.top.medications.gen2, "data") <- plot.df

plot(heatmap.top.medications.gen2)

unlist(lapply(medicine.indices$lvl5[["full"]][top.medications], length))


