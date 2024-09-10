

### Barplot smoking: current vs never ###

res.smoking.current.never.gen2 <- vector(mode = "list", length = length(organ.proteins))
names(res.smoking.current.never.gen2) <- names(organ.proteins)

res.smoking.current.never.gen2 <- lapply(res.smoking.current.never.gen2, function(x){
  return(data.frame(Estimate = NA,
                    `Std. Error` = NA,
                    `t value` = NA, 
                    df = NA,
                    `Pr(>|t|)` = NA, 
                    ll = NA, 
                    ul = NA, 
                    r = NA, 
                    check.names = FALSE))
})

df <- data.frame(
  eid = olink_bd_annotation_0$eid,
  smoking = NA, 
  chronological_age = olink_bd_annotation_0$age_first_visit, 
  biological_age = NA,
  sex = factor(olink_bd_annotation_0$p31, levels = c("Female", "Male"), ordered = FALSE),
  deprivation = olink_bd_annotation_0$p22189, # Townsend deprivation index
  centre = factor(olink_bd_annotation_0$p54_i0),
  IPAQ_activity = factor(olink_bd_annotation_0$p22032_i0, level = c("low", "moderate", "high"), ordered = FALSE)
)

df$smoking[olink_bd_annotation_0$p20116_i0 == "Current"] <- 1
df$smoking[olink_bd_annotation_0$p20116_i0 == "Never"] <- 0

for(k in 1:length(res.smoking.current.never.gen2)){
  
  if(names(organ.proteins)[k] != "Bladder"){
    
    df$biological_age <- predicted.ages.oof$gen2$predicted[[k]]

        model <- lm(biological_age ~ chronological_age * sex + smoking + deprivation + centre + IPAQ_activity, data = df)
        summary.vals <- tryCatch(summary(model)$coef["smoking", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
        res.smoking.current.never.gen2[[k]][, c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model, parm = "smoking"), sign(model$coef["smoking"])*sqrt(summary(model)$r.squared))
        res.smoking.current.never.gen2[[k]][, "df"] <- model$df.residual
    
    # Add q-value per group:
    res.smoking.current.never.gen2[[k]] <- res.smoking.current.never.gen2[[k]] %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))
    
  }
  
}

for(k in 1:length(res.smoking.current.never.gen2)){
  res.smoking.current.never.gen2[[k]]$organ <- names(res.smoking.current.never.gen2)[k]
}

res.smoking.current.never.gen2[["Bladder"]] <- NULL
res.df.smoking.current.never.gen2 <- do.call("rbind", res.smoking.current.never.gen2)
res.df.smoking.current.never.gen2$zval <- limma::zscoreT(res.df.smoking.current.never.gen2$`t value`, df = res.df.smoking.current.never.gen2$df, approx=FALSE, method = "bailey")

# saveRDS(res.df.smoking.current.never.gen2, file = paste0(rds.dir, "res_df_smoking_current_never_gen2.rds"))
res.df.smoking.current.never.gen2 <- readRDS(file = paste0(rds.dir, "res_df_smoking_current_never_gen2.rds"))

### Now make the barplot with whiskers ###

plot.df <- res.df.smoking.current.never.gen2[(res.df.smoking.current.never.gen2$organ %in% names(organ.proteins.selected)), ]
plot.df$qval <- p.adjust(plot.df$`Pr(>|t|)`, method = "BH")

plot.df$stars <- ""
plot.df$stars[plot.df$qval < 0.1] <- "'"
plot.df$stars[plot.df$qval < 0.05] <- "*"
plot.df$stars[plot.df$qval < 0.01] <- "**"
plot.df$stars[plot.df$qval < 0.001] <- "***"

plot.df <- plot.df %>% arrange(-abs(zval))

rescale_colors <- function(x){
  x <- c(0,x,1)
  out <- scales::rescale(log(x+0.0001))
  return(out[-c(1,length(out))])
}
labels_scale <- c(0, 0.001, 0.01, 0.05, 0.1, 1)

# To make sure that a zero after the comma is also printed out
plot.df$numbers <- unlist(rapply(list(plot.df$Estimate), sprintf, fmt = "%0.2f", how = "replace"))

x_pos_stars <- rep(NA, nrow(plot.df))
x_pos_stars[plot.df$Estimate > 0] <- plot.df$Estimate[plot.df$Estimate > 0]+0.1
x_pos_stars[plot.df$Estimate < 0] <- plot.df$Estimate[plot.df$Estimate < 0]-0.1

x_pos_numbers <- rep(NA, nrow(plot.df))
x_pos_numbers <- apply(cbind(plot.df$Estimate, 0), 1, mean)

main = paste0("")
xlim = c(-0.9, 0.9)

plot.df$organ <- factor(plot.df$organ, levels = rev(unique(plot.df$organ)))

barplot.smoking.current.never.gen2 <- ggplot_gtable(ggplot_build(
  ggplot(data = plot.df, aes(x = Estimate, y = organ)) +
    geom_bar(stat="identity", aes(fill = rescale_colors(qval))) + # fill = colors
    scale_fill_viridis(option = "plasma", name = "q-value", limits = c(0, 1), 
                       values = rescale_colors(labels_scale), 
                       breaks = rescale_colors(labels_scale), 
                       labels = labels_scale,
                       direction = -1,
                       guide = guide_colorbar(reverse = TRUE)) +
    geom_text(aes(label = numbers), x = x_pos_numbers, color="black", size=3, fontface = "bold") + # Not useful, all q-values of 0.01
    geom_text(aes(label = stars), x = x_pos_stars, color="black", size=5) + # , fontface = "bold"
    xlim(xlim[1], xlim[2]) +
    geom_vline(xintercept = 0, colour = "darkgrey") +
    geom_errorbarh(aes(xmax = ul, xmin = ll), height = 0.3) +
    ggtitle("Smoking") +
    xlab("Effect of smoking on predicted log(mortality hazard)") +
    ylab("") +
    theme_bw() +
    theme(panel.border = element_blank(),
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

attr(barplot.smoking.current.never.gen2, "data") <- plot.df

plot(barplot.smoking.current.never.gen2)

# log10(5) == log(5)/log(10)
# Get -log10(p-values)
-pt(abs(plot.df$`t value`), df = plot.df$df, lower.tail = FALSE, log.p = TRUE)/log(10)
# "Lung" "Conventional" "Immune" "Brain"  "Liver" "Kidney" "Intestine" "Skin" "Artery"  
# 1096.11903  755.56091  390.42154  363.09938  202.19484  127.88795  109.85434   85.65967   37.04795

### Drinking ###

res.drinking.gen2 <- vector(mode = "list", length = length(organ.proteins))
names(res.drinking.gen2) <- names(organ.proteins)

res.drinking.gen2 <- lapply(res.drinking.gen2, function(x){
  return(data.frame(Estimate = NA,
                    `Std. Error` = NA,
                    `t value` = NA, 
                    df = NA,
                    `Pr(>|t|)` = NA, 
                    ll = NA, 
                    ul = NA, 
                    r = NA, 
                    check.names = FALSE))
})


for(k in 1:length(res.drinking.gen2)){
  
  if(names(organ.proteins)[k] != "Bladder"){
    
    df <- data.frame(
      eid = olink_bd_annotation_0$eid,
      drinking = olink_bd_annotation_0$p26030_i0, 
      chronological_age = olink_bd_annotation_0$age_first_visit, 
      biological_age = predicted.ages.oof$gen2$predicted[[k]],
      sex = olink_bd_annotation_0$p31,
      deprivation = olink_bd_annotation_0$p22189, # Townsend deprivation index
      centre = factor(olink_bd_annotation_0$p54_i0),
      IPAQ_activity = factor(olink_bd_annotation_0$p22032_i0, level = c("low", "moderate", "high"), ordered = FALSE),
      smoking_status = factor(olink_bd_annotation_0$p20116_i0, levels = c("Prefer not to answer", "Never", "Previous", "Current"), ordered = FALSE)
    )
        model <- lm(biological_age ~ chronological_age * sex + drinking + deprivation + centre + IPAQ_activity + smoking_status, data = df)
        # summary(model)
        summary.vals <- tryCatch(summary(model)$coef["drinking", c("Estimate", "Std. Error", "t value", "Pr(>|t|)")], error = function(e){return(rep(NA, 4))})
        res.drinking.gen2[[k]][, c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "ll", "ul", "r")] <- c(summary.vals, confint(model, parm = "drinking"), sign(model$coef["drinking"])*sqrt(summary(model)$r.squared))
        res.drinking.gen2[[k]][, "df"] <- model$df.residual
    
    # Add q-value per group:
    res.drinking.gen2[[k]] <- res.drinking.gen2[[k]] %>% mutate(qval = p.adjust(`Pr(>|t|)`, method = "BH"))
    
  }
  
}

for(k in 1:length(res.drinking.gen2)){
  res.drinking.gen2[[k]]$organ <- names(res.drinking.gen2)[k]
}

res.drinking.gen2[["Bladder"]] <- NULL
res.df.drinking.gen2 <- do.call("rbind", res.drinking.gen2)

# saveRDS(res.df.drinking.gen2, file = paste0(rds.dir, "res_df_drinking_gen2.rds"))
res.df.drinking.gen2 <- readRDS(file = paste0(rds.dir, "res_df_drinking_gen2.rds"))

### Now make the barplot with whiskers ###

plot.df <- res.df.drinking.gen2[(res.df.drinking.gen2$organ %in% names(organ.proteins.selected)), ]
plot.df$qval <- p.adjust(plot.df$`Pr(>|t|)`, method = "BH")

plot.df$organ <- paste(toupper(substr(plot.df$organ, 1, 1)), substr(plot.df$organ, 2, nchar(plot.df$organ)), sep="")

plot.df$stars <- ""
plot.df$stars[plot.df$qval < 0.1] <- "'"
plot.df$stars[plot.df$qval < 0.05] <- "*"
plot.df$stars[plot.df$qval < 0.01] <- "**"
plot.df$stars[plot.df$qval < 0.001] <- "***"

plot.df <- plot.df %>% arrange(`Pr(>|t|)`)

rescale_colors <- function(x){
  x <- c(0,x,1)
  out <- scales::rescale(log(x+0.0001))
  return(out[-c(1,length(out))])
}
labels_scale <- c(0, 0.001, 0.01, 0.05, 0.1, 1)

# To make sure that a zero after the comma is also printed out
plot.df$numbers <- unlist(rapply(list(plot.df$Estimate), sprintf, fmt = "%0.4f", how = "replace"))

x_pos_stars <- rep(NA, nrow(plot.df))
x_pos_stars[plot.df$Estimate > 0] <- plot.df$Estimate[plot.df$Estimate > 0]+0.0002
x_pos_stars[plot.df$Estimate < 0] <- plot.df$Estimate[plot.df$Estimate < 0]-0.0002

x_pos_numbers <- rep(NA, nrow(plot.df))
x_pos_numbers[plot.df$Estimate > 0] <- plot.df$ul[plot.df$Estimate > 0]+0.0003
x_pos_numbers[plot.df$Estimate < 0] <- plot.df$ll[plot.df$Estimate < 0]-0.00035

main = paste0("")
xlim = c(-0.002, 0.0025)

plot.df$organ <- factor(plot.df$organ, levels = rev(unique(plot.df$organ)))

barplot.drinking.gen2 <- ggplot_gtable(ggplot_build(
  ggplot(data = plot.df, aes(x = Estimate, y = organ)) +
    geom_bar(stat="identity", aes(fill = rescale_colors(qval))) + # fill = colors
    scale_fill_viridis(option = "plasma", name = "q-value", limits = c(0, 1), 
                       values = rescale_colors(labels_scale), 
                       breaks = rescale_colors(labels_scale), 
                       labels = labels_scale,
                       direction = -1,
                       guide = guide_colorbar(reverse = TRUE)) +
    geom_text(aes(label = numbers), x = x_pos_numbers, color="black", size=3, fontface = "bold") + 
    geom_text(aes(label = stars), x = x_pos_stars, color="black", size=5) + 
    xlim(xlim[1], xlim[2]) +
    geom_vline(xintercept = 0, colour = "darkgrey") +
    geom_errorbarh(aes(xmax = ul, xmin = ll), height = 0.3) +
    ggtitle("Drinking") +
    xlab("Effect per g of alcohol on predicted log(mortality hazard)") +
    ylab("") +
    theme_bw() +
    theme(panel.border = element_blank(),
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

attr(barplot.drinking.gen2, "data") <- plot.df

plot(barplot.drinking.gen2)





