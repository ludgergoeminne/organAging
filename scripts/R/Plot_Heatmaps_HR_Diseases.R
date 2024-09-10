
### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

# See also 1_Prepare_Data.R
general_hazard_outcomes <- get_general_outcomes()

# 2. calculate HR as in res.organs.list from organ_specific_HR_plots.R
# Now we calculate the effects of organ-specific age minus conventional age

res.difforgan.list <- vector(mode = "list", length = 2)
names(res.difforgan.list) <- c("gen1", "gen2")

res.difforgan.list <- lapply(res.difforgan.list, function(x){
  x <- vector(mode = "list", length = length(organ.proteins))
  names(x) <- names(organ.proteins)
  
  x <- lapply(x, function(y){
    return(data.frame(outcome = names(general_hazard_outcomes),
                      `exp(coef)` = NA,
                      z = NA,
                      `Pr(>|z|)` = NA, 
                      ll = NA, 
                      ul = NA, 
                      check.names = FALSE))
  })
  
  return(x)
} )




for(g in 1:length(res.difforgan.list)){
  
  for(k in 1:length(res.difforgan.list[[g]])){
    
    # Conventional cannot be included now because we aim to regress it out.
    if(!(names(res.difforgan.list[[g]])[k] %in% c("Bladder", "Conventional"))){
      
      for(i in 1:length(names(general_hazard_outcomes))){
        
        df <- data.frame(
          eid = olink_bd_annotation_0$eid,
          status_os = olink_bd_annotation_0[, paste0("status_", names(general_hazard_outcomes)[i])],
          time = olink_bd_annotation_0[, paste0("time_", names(general_hazard_outcomes)[i])],
          chronological_age = olink_bd_annotation_0$age_first_visit, # == df.test$age_first_visit,
          biological_age = predicted.ages.oof[[g]]$predicted[[k]],
          biological_age_conv = predicted.ages.oof[[g]]$predicted[["Conventional"]],
          sex = olink_bd_annotation_0$p31
        )
        
            cox.model <- coxph(Surv(time = time, event = status_os, type='right') ~ chronological_age * sex + biological_age + biological_age_conv, data = df)
            # summary(cox.model)
            if(is.na(summary(cox.model)$coef["biological_age_conv", "coef"])){stop("This model did NOT correct for conventional age!")}
            res.difforgan.list[[g]][[k]][(res.difforgan.list[[g]][[k]]$outcome == names(general_hazard_outcomes)[i]), c("exp(coef)", "z", "Pr(>|z|)", "ll", "ul")] <- c(summary(cox.model)$coef["biological_age", c("exp(coef)", "z", "Pr(>|z|)")], exp(confint(cox.model, parm = "biological_age")))
      }
    }
    
  }
  
}

# Add q-value per group:
for(g in 1:length(res.difforgan.list)){
  res.difforgan.list[[g]] <- lapply(res.difforgan.list[[g]], function(x){
    x <- x %>% mutate(qval = p.adjust(`Pr(>|z|)`, method = "BH"))
    return(x)
  })
}

res.difforgan.list[["gen2"]][["Heart"]]
res.difforgan.list[["gen2"]][["Brain"]]
res.difforgan.list[["gen2"]][["Lung"]]
res.difforgan.list[["gen2"]][["Intestine"]]

# saveRDS(res.difforgan.list, file = paste0(rds.dir, "res_difforgan_list.rds"))
res.difforgan.list <- readRDS(file = paste0(rds.dir, "res_difforgan_list.rds"))


### Heatmap 2: statistics ###
### A very interesting plot !!! ###

res.difforgan.list$gen1$Conventional <- NULL
res.difforgan.list$gen2$Conventional <- NULL

res.difforgan.list <- lapply(res.difforgan.list, function(x){
  for(k in 1:length(x)){
    x[[k]]$organ <- names(x)[k]
  }
  return(x)
})



res.difforgan.df <- lapply(res.difforgan.list, function(x){
  x <- do.call("rbind", x)
  x <- x[x$organ %in% names(organ.proteins.selected),]
  return(x)
})



heatmaps.HR.disease <- vector(mode = "list", length = 2)
names(heatmaps.HR.disease) <- c("gen1", "gen2")

for(g in 1:length(heatmaps.HR.disease)){
  
  plot.df <- res.difforgan.df[[g]]
  plot.df$coef <- log(plot.df$`exp(coef)`)
  
  plot.df$outcome <- gsub("_", " ", plot.df$outcome)
  plot.df$outcome <- paste(toupper(substr(plot.df$outcome, 1, 1)), substr(plot.df$outcome, 2, nchar(plot.df$outcome)), sep="")
  plot.df$outcome[plot.df$outcome == "Cirrhosis fibrosis"] <- "Liver cirrhosis/fibrosis"
  
  ### Clustering ###
  data <- scale(as.matrix((plot.df[,c("outcome", "organ", "coef")] %>% pivot_wider(names_from = "outcome", values_from = "coef"))[, -1]))
  rownames(data) <- unique(plot.df$organ)
  clust.organ <- hclust(dist(data, method = "euclidean"), method = "ward.D")
  ord <- clust.organ$order
  ord
  plot.df$organ = factor(plot.df$organ, levels = unique(plot.df$organ)[ord])
  
  data <- scale(as.matrix((plot.df[,c("organ", "outcome", "coef")] %>% pivot_wider(names_from = "organ", values_from = "coef"))[, -1]))
  rownames(data) <- unique(plot.df$outcome)
  clust.outcome <- hclust(dist(data, method = "euclidean"), method = "ward.D")
  ord <- clust.outcome$order
  ord
  plot.df$outcome = factor(plot.df$outcome, levels = rev(unique(plot.df$outcome)[ord]))
  clust.outcome$order <- rev(clust.outcome$order)
  
  plot.df$stars <- ""
  plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
  plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
  plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
  plot.df$stars[which(plot.df$qval < 0.001)] <- "***"
  
  lims <- c(-1,1)
  
  range(plot.df$coef, na.rm = TRUE)
  # -0.4140808  0.8994797
  
  c(quantile(plot.df$coef, 0.05, na.rm = TRUE), quantile(plot.df$coef, 0.95, na.rm = TRUE))
  # -0.2528320  0.7908659
  
  Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  
  tmp <- plot.df # Needed for Excel file, where we don't want the numeric values
  rownames(tmp) <- NULL
  # Needed for geom_dendro!!!
  plot.df$organ <- as.numeric(plot.df$organ)
  plot.df$outcome <- as.numeric(plot.df$outcome)
  
  heatmaps.HR.disease[[g]] <- ggplot_gtable(ggplot_build(
    ggplot(plot.df, aes(x = organ, y = outcome, label = stars)) + 
      geom_tile(aes(fill = coef)) + 
      scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
      guides(color = "none") + # no legend for color
      geom_dendro(clust.organ, ylim=c(0, -1)) + #upper dendrogram
      geom_dendro(clust.outcome, xlim=c(0, -1), pointing="side") + #side dendrogram
      labs(fill = "Log(hazard ratio)") + 
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
            plot.title = element_blank(), 
            axis.title.x = element_text(size = 12, color = "black"),
            axis.title.y = element_text(size = 12, color = "black"),
            axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.line.x = element_blank(), 
            axis.line.y = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.ticks.y = element_blank(), 
            strip.text = element_text(face = "bold",size = 12),
            strip.background = element_blank())
  ))
  
  attr(heatmaps.HR.disease[[g]], "data") <- tmp # I.e. plot.df, but without the numeric values
  
}

plot(heatmaps.HR.disease$gen1)
plot(heatmaps.HR.disease$gen2)

