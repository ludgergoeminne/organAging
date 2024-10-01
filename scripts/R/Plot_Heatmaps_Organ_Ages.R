
### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
predicted.ages.oof <- readRDS(file = paste0(rds.dir, "predicted_ages_oof.rds"))

### This comes from 1_Prepare_Data.R ###
olink_bd_annotation_0 <- readRDS(file = paste0(rds.dir, "olink_bd_annotation_0.rds"))

### This comes from 3_Process_GTEx.R ###
organ.proteins <- readRDS(file = paste0(rds.dir, "organ_proteins.rds"))

### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

# Import PhenoAge and Locomotor Age
tmp <- read.csv(file = paste0(external.models.dir, phenoage.file), header = TRUE, quote = "")
PhenoAge <- tmp$PhenoAge
names(PhenoAge) <- as.character(tmp$UID)
resid.PhenoAge <- resid(lm(PhenoAge[as.character(olink_bd_annotation_0[, "eid"])]~olink_bd_annotation_0[, "age_first_visit"]))

tmp <- read.csv(file = paste0(external.models.dir, locoage.file), header = TRUE, quote = "")
LocoAge <- tmp$LocoBioage
names(LocoAge) <- as.character(tmp$UID)
resid.LocoAge <- resid(lm(LocoAge[as.character(olink_bd_annotation_0[, "eid"])]~olink_bd_annotation_0[, "age_first_visit"]))

### 1. Create correlation data frames ###

res.df.organ.list <- vector(mode = "list", length = 2)
names(res.df.organ.list) <- c("gen1", "gen2")

res.df.organ.list <- lapply(res.df.organ.list, function(x){
  x <- vector(mode = "list", length = 2)
  names(x) <- c("predicted", "residual")
  return(x)
} )


for(g in 1:length(res.df.organ.list)){
  
  for(j in 1:length(res.df.organ.list[[g]])){
    
    bio.ages <- as.data.frame(do.call("cbind", predicted.ages.oof[[g]][[j]]))
    # cbind here chronological age (not for residuals), phenoAge,...!!!!
    if(names(res.df.organ.list[[g]])[j] != "residual"){
      bio.ages$`chronological age` <- olink_bd_annotation_0[, "age_first_visit"]
      # cbind here PhenoAge and LocoAge
      bio.ages$`phenotypic age` <- PhenoAge[as.character(olink_bd_annotation_0[, "eid"])]
      bio.ages$LocoAge <- LocoAge[as.character(olink_bd_annotation_0[, "eid"])]
    } else{
      # cbind here residual PhenoAge and LocoAge
      bio.ages$`phenotypic age` <- resid.PhenoAge[as.character(olink_bd_annotation_0[, "eid"])]
      bio.ages$LocoAge <- resid.LocoAge[as.character(olink_bd_annotation_0[, "eid"])]
    }
    
    colnames(bio.ages) <- paste(toupper(substr(colnames(bio.ages), 1, 1)), substr(colnames(bio.ages), 2, nchar(colnames(bio.ages))), sep="")
    colnames(bio.ages) <- gsub("_", " ", colnames(bio.ages))
    correlation.traits <- colnames(bio.ages)[!(colnames(bio.ages) %in% c("Bladder"))]
    
    res.df.organ.list[[g]][[j]] <- data.frame(trait.x = rep(correlation.traits, each = length(correlation.traits)),
                                              trait.y = rep(correlation.traits, times = length(correlation.traits)),
                                              Estimate = NA,
                                              `Std. Error` = NA,
                                              `t value` = NA, 
                                              df = NA, 
                                              `Pr(>|t|)` = NA, 
                                              r = NA, 
                                              check.names = FALSE)
    
    for(k in 1:nrow(res.df.organ.list[[g]][[j]])){
      model <- lm(formula(get(res.df.organ.list[[g]][[j]]$trait.x[k])~get(res.df.organ.list[[g]][[j]]$trait.y[k])), data = bio.ages)
      summary <- summary(model)
      res.df.organ.list[[g]][[j]][k, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")] <- summary$coef["get(res.df.organ.list[[g]][[j]]$trait.y[k])", ]
      res.df.organ.list[[g]][[j]][k, c("df")] <- model$df.residual
      res.df.organ.list[[g]][[j]][k, c("r")] <- sign(summary$coef["get(res.df.organ.list[[g]][[j]]$trait.y[k])", "Estimate"])*sqrt(summary$r.squared)
    }
    
    res.df.organ.list[[g]][[j]]$zval <- limma::zscoreT(res.df.organ.list[[g]][[j]]$`t value`, df = res.df.organ.list[[g]][[j]]$df, approx=FALSE, method = "bailey") # qnorm(res.df.organ.list[[g]][[j]]$`Pr(>|t|)`, lower.tail = FALSE) * sign(res.df.organ.list[[g]][[j]]$Estimate)
    res.df.organ.list[[g]][[j]] <- res.df.organ.list[[g]][[j]] %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))
    
  }
}

# saveRDS(res.df.organ.list, file = paste0(rds.dir, "res_df_organ_list.rds"))
res.df.organ.list <- readRDS(file = paste0(rds.dir, "res_df_organ_list.rds"))

### Number of observations ###
sum(!is.na(PhenoAge[as.character(olink_bd_annotation_0[, "eid"])]))
# 38025
sum(!is.na(LocoAge[as.character(olink_bd_annotation_0[, "eid"])]))
# 8660

### 2. Make the correlation plots ###

heatmap.corr.organs <- vector(mode = "list", length = length(res.df.organ.list))
names(heatmap.corr.organs) <- names(res.df.organ.list)

heatmap.corr.organs <- lapply(heatmap.corr.organs, function(x){
  x <- vector(mode = "list", length = length(res.df.organ.list)[[1]])
  names(x) <- names(res.df.organ.list[[1]])
  return(x)
} )


### Change "LocoAge" to "Locomotor age"
res.df.organ.list <- lapply(res.df.organ.list, function(x){
  x <- lapply(x, function(y){
    y[y$trait.x == "LocoAge", "trait.x"] <- "Locomotor age"
    y[y$trait.y == "LocoAge", "trait.y"] <- "Locomotor age"
    return(y)
  })
  return(x)
})

for(g in 1:length(res.df.organ.list)){
  for(j in 1:length(res.df.organ.list[[g]])){
    
    plot.df <- res.df.organ.list[[g]][[j]]
    
    ### Select subset ###
    to.exclude <- names(organ.proteins)[!(names(organ.proteins) %in% names(organ.proteins.selected))]
    to.exclude <- paste(toupper(substr(to.exclude, 1, 1)), substr(to.exclude, 2, nchar(to.exclude)), sep="")
    to.exclude <- gsub("_", " ", to.exclude)
    plot.df <- plot.df[!(plot.df$trait.x %in% to.exclude), ]
    plot.df <- plot.df[!(plot.df$trait.y %in% to.exclude), ]
    ########################
    
    plot.df$trait.x <- gsub("_", " ", plot.df$trait.x)
    plot.df$trait.x <- paste(toupper(substr(plot.df$trait.x, 1, 1)), substr(plot.df$trait.x, 2, nchar(plot.df$trait.x)), sep="")
    plot.df$trait.y <- gsub("_", " ", plot.df$trait.y)
    plot.df$trait.y <- paste(toupper(substr(plot.df$trait.y, 1, 1)), substr(plot.df$trait.y, 2, nchar(plot.df$trait.y)), sep="")
    
    ### Clustering ###
    data <- scale(as.matrix((plot.df[,c("trait.x", "trait.y", "r")] %>% pivot_wider(names_from = "trait.x", values_from = "r"))[, -1]))
    rownames(data) <- colnames(data)
    clust.trait.x <- hclust(dist(data, method = "euclidean"), method = "ward.D")
    ord <- clust.trait.x$order
    ord
    plot.df$trait.x = factor(plot.df$trait.x, levels = unique(plot.df$trait.x)[ord])
    
    data <- scale(as.matrix((plot.df[,c("trait.y", "trait.x", "r")] %>% pivot_wider(names_from = "trait.y", values_from = "r"))[, -1]))
    rownames(data) <- colnames(data)
    clust.trait.y <- hclust(dist(data, method = "euclidean"), method = "ward.D")
    ord <- clust.trait.y$order
    ord
    plot.df$trait.y = factor(plot.df$trait.y, levels = rev(unique(plot.df$trait.y)[ord]))
    clust.trait.y$order <- rev(clust.trait.y$order)
    
    ### Plot only the lower triangle ###
    plot.df <- plot.df[order(plot.df$trait.x),]
    plot.df <- plot.df[order(plot.df$trait.y),]
    pasted.traits <- apply(plot.df[, c("trait.x", "trait.y")], 1, function(x){return(sort(unlist(x)))}) %>% apply(2, function(x){return(paste0(x, collapse = ""))})
    plot.df <- plot.df[which(!duplicated(pasted.traits)),]
    
    # plot.df <- plot.df %>% arrange(`exp(coef)`) %>% arrange(set)
    # plot.df$outcome <- factor(plot.df$outcome, levels = unique(plot.df$outcome))
    
    plot.df$stars <- ""
    plot.df$stars[which(plot.df$qval < 0.1)] <- "'"
    plot.df$stars[which(plot.df$qval < 0.05)] <- "*"
    plot.df$stars[which(plot.df$qval < 0.01)] <- "**"
    plot.df$stars[which(plot.df$qval < 0.001)] <- "***"
    
    lims <- c(-1,1)
    
    range(plot.df$r, na.rm = TRUE)
    # -0.02973752  1.00000000
    
    c(quantile(plot.df$r, 0.05, na.rm = TRUE), quantile(plot.df$r, 0.95, na.rm = TRUE))
    # 0.004233909 0.926340640 
    
    # Only here change to factor!
    # plot.df$trait.x = factor(plot.df$trait.x, levels = unique(plot.df$trait.x))
    # plot.df$trait.y = factor(plot.df$trait.y, levels = rev(unique(plot.df$trait.y)))
    
    # Turn the plot into a lower triangle
    # plot.df <- plot.df[nrow(plot.df):1, ]
    plot.df <- plot.df[!duplicated(apply(apply(plot.df[, c("trait.x", "trait.y")], 1, function(x){sort(x)}), 2, function(x){return(paste0(x, collapse = ""))})),]
    
    Heatmap_palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "white","#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
    
    tmp <- plot.df # Needed for Excel file, where we don't want the numeric values
    rownames(tmp) <- NULL
    # Needed for geom_dendro!!!
    plot.df$trait.x <- as.numeric(plot.df$trait.x)
    plot.df$trait.y <- as.numeric(plot.df$trait.y)
    
    heatmap.corr.organs[[g]][[j]] <- ggplot_gtable(ggplot_build(
      ggplot(plot.df, aes(x = trait.x, y = trait.y, label = stars)) + 
        geom_tile(aes(fill = r)) + # size = minuslog10p , color = p.adjust
        # scale_shape_manual(values=c(22, 24, 25)) +
        # facet_grid(pathway~diet, scales = "free_y", space="free") +
        scale_fill_gradientn(colours = Heatmap_palette, oob = scales::squish, na.value = 'darkgrey', limits = lims) + # , limits = lims
        # theme_bw(base_size = 14) +
        # scale_fill_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + # , breaks = breaks, labels = 10^(-breaks)
        # scale_color_viridis(option = "plasma", breaks = breaks, labels = 10^(-breaks), limits = c(0, 12)) + 
        guides(color = "none") + # no legend for color
        geom_dendro(clust.trait.x, ylim=c(0, -1)) +                             #lower dendrogram
        geom_dendro(clust.trait.y, xlim=c(0, -1), pointing="side") +             #side dendrogram
        labs(fill="r") + # size="-log10(q-value)", # fgsea uses "BH" to correct for multiple testing
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
              axis.line.x = element_blank(), # element_line(color="black", linewidth = 0.5),
              axis.line.y = element_blank(), # element_line(color="black", linewidth = 0.5),
              axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
              axis.ticks.y = element_blank(), # element_line(color="black", linewidth = 0.5),
              strip.text = element_text(face = "bold",size = 12),
              strip.background = element_blank())
    ))
    attr(heatmap.corr.organs[[g]][[j]], "data") <- tmp # I.e. plot.df, but without the numeric values
  }
}

plot(heatmap.corr.organs$gen1$residual)
plot(heatmap.corr.organs$gen2$residual)
plot(heatmap.corr.organs$gen1$predicted)
plot(heatmap.corr.organs$gen2$predicted)
