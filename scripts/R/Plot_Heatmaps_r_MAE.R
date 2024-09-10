
### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from 4_calculate_predicted_resid_ages.R ###
coefficients <- readRDS(paste0(rds.dir, "coefficients.rds"))
predicted.ages.1 <- readRDS(file = paste0(rds.dir, "predicted_ages_1.rds"))
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))
training.eids <- readRDS(file = paste0(rds.dir, "training_eids.rds"))
test.eids <- readRDS(file = paste0(rds.dir, "test_eids.rds"))

### This comes from 1_Prepare_Data.R ###
sds <- readRDS(file = paste0(rds.dir, "standard_deviations.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))


### Import the MESA dataset (Bild et al., 2002) ###
meta <- read.table(paste0(input.MESA.dir, "phs001416.v3.pht010511.v1.p1.c1.TOPMed_MESA_Proteomics_Sample_Attributes.HMB.txt"), 
                   sep = "\t", skip = 10, header = TRUE)

prot <- read.table(paste0(input.MESA.dir, "MESA_ProteomicsDataMergedRunlistKey_ProteinInfo_All_20220426.txt"), 
                   sep = "\t", header = TRUE, quote = "")

MESA.df <- read.table(paste0(input.MESA.dir, "MESA_ProteomicsDataMergedRunlistKey_DS_HMB_20220426.txt"), 
                 sep = "\t", header = TRUE, quote = "")

MESA.df[1:10, 1:10]
all(colnames(MESA.df)[-1] == prot$SomaId)
# TRUE
colnames(MESA.df)[-1] <- prot$EntrezGeneSymbol
MESA.df[1:10, 1:10]

# Rescale the data because MESA is SomaLogic data!
MESA.df[, -1] <- log2(MESA.df[, -1])
tmp2 <- apply(MESA.df[, -1], 2, function(x){return((x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))})
MESA.df <- cbind(data.frame(TOP_ID = MESA.df[, 1]), sweep(as.matrix(tmp2[, names(sds)[names(sds) %in% colnames(tmp2)]]), MARGIN=2, sds[names(sds) %in% colnames(tmp2)], `*`))

# Arrange meta by TOP_ID from df
meta <- meta[match(MESA.df$TOP_ID, meta$SAMPLE_ID),]
all(meta$SAMPLE_ID == MESA.df$TOP_ID)

sum(colnames(MESA.df)[-1] %in% names(coefficients$gen1$Conventional[coefficients$gen1$Conventional != 0]))/sum(coefficients$gen1$Conventional != 0)
# 27% identified
sum(colnames(MESA.df)[-1] %in% names(coefficients$gen2$Conventional[coefficients$gen2$Conventional != 0]))/sum(coefficients$gen2$Conventional != 0)
# 39% identified

### Put the data together in a data frame ###

sets <- c("Training set", "Test set", "Out of fold", "MESA visit 1", "MESA visit 2")

summary.list <- vector(mode = "list", length = 2)
names(summary.list) <- c("1st-generation models", "Mortality-based models")

summary.list <- lapply(summary.list, function(x){
  x <- vector(mode = "list", length = 5)
  names(x) <- c("r", "r²", "MAE", "MSE", "pval")
  
  x <- lapply(x, function(x){
    y <- data.frame(
      organ = rep(names(organ.proteins), each = length(sets)),
      set = rep(sets, times = length(organ.proteins)),
      value = NA
    )
    return(y)
  })
  
  return(x)
})

for(g in 1:length(summary.list)){
    
    for(k in 1:length(organ.proteins)){
      if(names(organ.proteins)[k] != "Bladder"){
      
        for(i in 1:length(sets)){

          if(sets[i] == "Training set"){
            
            df <- data.frame(
              chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% training.eids)],
              biological_age = predicted.ages.1[[g]]$predicted[[k]][training.eids])
            
          } else if(sets[i] == "Test set"){
            
            df <- data.frame(
              chronological_age = olink_bd_annotation_0$age_first_visit[(olink_bd_annotation_0$eid %in% test.eids)],
              biological_age = predicted.ages.1[[g]]$predicted[[k]][test.eids])
            
          } else if(sets[i] == "Out of fold"){
            
            df <- data.frame(
              chronological_age = olink_bd_annotation_0$age_first_visit,
              biological_age = predicted.ages.oof[[g]]$predicted[[k]])
            
          } else if(sets[i] %in% c("MESA visit 1", "MESA visit 2")){
            
            if(sets[i] == "MESA visit 1"){
              df.tmp <- MESA.df[meta$COLLECTION_VISIT == 1,]
              meta.tmp <- meta[meta$COLLECTION_VISIT == 1,]
            } else if(sets[i] == "MESA visit 2"){
              df.tmp <- MESA.df[meta$COLLECTION_VISIT == 5,]
              meta.tmp <- meta[meta$COLLECTION_VISIT == 5,]
            }
            if(!all(meta.tmp$SAMPLE_ID == df.tmp$TOP_ID)){stop("MESA dataset is not properly ordered!")}
            
            df <- data.frame(
              chronological_age = meta.tmp$AGE_AT_COLLECTION
              )
            if(g == 1){
              coefs <- coefficients[[g]][[k]][coefficients[[g]][[k]] != 0]
              coefs <- coefs[names(coefs) %in% colnames(df.tmp)]
              intercept <- coefficients[[g]][[k]]["Intercept"]
              
              tmp <- df.tmp[, colnames(df.tmp) %in% names(coefs), drop = FALSE]
              tmp <- tmp[, names(coefs), drop = FALSE]
              
              df$biological_age = intercept + rowSums(sweep(as.matrix(tmp), MARGIN=2, coefs, `*`))
              
            } else if(g == 2){
              coefs <- coefficients[[g]][[k]][coefficients[[g]][[k]] != 0]
              coefs <- coefs[names(coefs) %in% colnames(df.tmp)]
              
              tmp <- df.tmp[, colnames(df.tmp) %in% names(coefs), drop = FALSE]
              tmp <- tmp[, names(coefs), drop = FALSE]
              
              df$biological_age = rowSums(sweep(as.matrix(tmp), MARGIN=2, coefs, `*`))
              
            }
            
          }
          
          model <- lm(biological_age ~ chronological_age, data = df)
          summary <- summary(model)
          
          for(j in 1:length(summary.list[[g]])){
            
            if(names(summary.list[[g]])[j] == "r"){
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- sign(summary$coefficients["chronological_age", "Estimate"])*sqrt(summary$r.squared)
            } else if(names(summary.list[[g]])[j] == "r²"){
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- summary$r.squared
            } else if(names(summary.list[[g]])[j] == "MAE"){
              resid <- resid(model)
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- mean(abs(resid))
            } else if(names(summary.list[[g]])[j] == "MSE"){
              resid <- resid(model)
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- mean(resid^2)
            } else if(names(summary.list[[g]])[j] == "pval"){
              summary.list[[g]][[j]][(summary.list[[g]][[j]][, "organ"] == names(organ.proteins)[k]) & (summary.list[[g]][[j]][, "set"] == sets[i]), "value"] <- summary$coefficients["chronological_age", "Pr(>|t|)"]
            }
            
          }
          
        }
      
      }
    }

}

### Get number of participants ###
length(training.eids)
# 35962
length(test.eids)
# 8990
dim(MESA.df[meta$COLLECTION_VISIT == 1,])
# 921 756
dim(MESA.df[meta$COLLECTION_VISIT == 5,])
# 921 756

summary.df.list <- lapply(summary.list, function(x){
  x <- do.call("cbind", x)
  x <- x[, c(1, 2, seq(3, ncol(x), by = 3))]
  colnames(x)[c(1, 2)] <- gsub("^.+?\\.", "", colnames(x)[c(1, 2)])
  colnames(x)[-c(1, 2)] <- gsub("\\..+$", "", colnames(x)[-c(1, 2)])
  return(x)
})

# Needed for Export_Supp_Tables.R
# saveRDS(summary.df.list, file = paste0(rds.dir, "summary_df_list.rds"))
summary.df.list <- readRDS(file = paste0(rds.dir, "summary_df_list.rds"))

### Plot the correlation heatmaps for the non-longitudinal models ###

heatmaps.r <- vector(mode = "list", length = length(summary.df.list))
names(heatmaps.r) <- names(summary.df.list)

for(g in 1:length(heatmaps.r)){
  
  plot.df <- summary.df.list[[g]]
  plot.df <- plot.df[plot.df$organ %in% names(organ.proteins.selected),]
  
  # Add q-value per group:
  plot.df <- plot.df %>% mutate(qval = p.adjust(`pval`, method = "BH"))
  
  plot.df$stars <- ""
  plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
  plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
  plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
  plot.df$stars[which(plot.df$qval < 0.001)] <- "***"
  
  lims <- c(-1, 1)
  
  # Only here change to factor!
  plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))
  plot.df$set = factor(plot.df$set, levels = rev(unique(plot.df$set)))
  
  Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  
  heatmaps.r[[g]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = organ, y = set, label = stars)) + 
      geom_tile(aes(fill = r)) + # size = minuslog10p , color = p.adjust
      scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
      guides(color = "none") + # no legend for color
      labs(fill = "r") + 
      geom_text(size = 4.5, nudge_y = -0.13, na.rm = TRUE) +
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
  attr(heatmaps.r[[g]], "data") <- plot.df
}

plot(heatmaps.r[[1]])
plot(heatmaps.r[[2]])

### Plot the MAE heatmaps for the non-longitudinal models ###

heatmaps.MAE <- vector(mode = "list", length = length(summary.df.list))
names(heatmaps.MAE) <- names(summary.df.list)

lims.max <- c(3.5, 0.6)

for(g in 1:length(heatmaps.MAE)){
  
  plot.df <- summary.df.list[[g]]
  plot.df <- plot.df[plot.df$organ %in% names(organ.proteins.selected),]
  
  lims <- c(0, lims.max[g])
  
  # Only here change to factor!
  plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ))
  plot.df$set = factor(plot.df$set, levels = rev(unique(plot.df$set)))
  
  Heatmap_palette <- c("white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  
  heatmaps.MAE[[g]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = organ, y = set)) + 
      geom_tile(aes(fill = MAE)) + # size = minuslog10p , color = p.adjust
      scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
      guides(color = "none") + # no legend for color
      labs(fill = "MAE") + 
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
  attr(heatmaps.MAE[[g]], "data") <- plot.df
}

plot(heatmaps.MAE[[1]])
plot(heatmaps.MAE[[2]])

