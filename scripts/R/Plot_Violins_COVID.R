
### The Filbin et al. (2021) COVID dataset ###

### This comes from 4_Calculate_Predicted_Resid_Ages.R ###
coefficients.longitudinal <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal.rds"))

### This comes from 2_Prepare_Longitudinal_data.R ###
longitudinal.proteins <- readRDS(file = paste0(rds.dir, "longitudinal_proteins.rds"))
# sds.longitudinal <- readRDS(file = "standard_deviations_longitudinal.rds")

# Origin of the data:
# https://data.mendeley.com/datasets/nf853r8xsj/2
# https://www.sciencedirect.com/science/article/pii/S2666379121001154?via%3Dihub
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8091031/

# Clinical_Metadata.xlsx tab "Annotations" gives the info
# Clinical_Metadata.xlsx tab "Subject-level metadata" gives the clincal metadata
# Olink_Proteomics.xlsx tab "Olink Proteomics" contains the Olink normalized values
Clinical_Metadata <- openxlsx::read.xlsx(paste0(input.COVID.dir, "Clinical_Metadata.xlsx"), sheet = "Subject-level metadata")
Olink_Proteomics <- openxlsx::read.xlsx(paste0(input.COVID.dir, "Olink_Proteomics.xlsx"), sheet = "Olink Proteomics")

conversion.COVID <- openxlsx::read.xlsx(paste0(input.COVID.dir, "Supplemental-Table-3-Olink-Models-All.xlsx"), sheet = "Global F-test")
length(unique(conversion.COVID$OlinkID))
# 1429
conversion.COVID <- conversion.COVID[, c("Assay", "OlinkID")]
conversion.COVID <- conversion.COVID[!duplicated(conversion.COVID$OlinkID),]
length(conversion.COVID$OlinkID)
# 1429
all(colnames(Olink_Proteomics)[-c(1,2)] %in% conversion.COVID$OlinkID)
# TRUE
all(conversion.COVID$OlinkID %in% colnames(Olink_Proteomics)[-c(1,2)])
# TRUE

conversion.COVID <- conversion.COVID[match(colnames(Olink_Proteomics)[-c(1,2)], conversion.COVID$OlinkID),]
dim(conversion.COVID)
# TRUE
all(conversion.COVID$OlinkID == colnames(Olink_Proteomics)[-c(1,2)])
# TRUE

colnames(Olink_Proteomics)[-c(1,2)] <- conversion.COVID$Assay
Olink_Proteomics[1:10,1:10]

# Prepare metadata
tmp <- cbind(data.frame(Unique.ID = Olink_Proteomics$Public.ID), Olink_Proteomics[, c(1,2)])
tmp$Public.ID <- as.numeric(gsub("_D.*", "", tmp$Public.ID)) # To be able to merge with Clinical_Metadata
# Clinical_Metadata$Public.ID[!(Clinical_Metadata$Public.ID %in% tmp$Public.ID)]
# "180": one missing ID
metadata.COVID <- dplyr::left_join(tmp, Clinical_Metadata, c("Public.ID" = "Public.ID"))
dim(metadata.COVID)
# 784  46
all(Olink_Proteomics$Public.ID == metadata.COVID$Unique.ID)
# TRUE

Olink_Proteomics <- Olink_Proteomics[, -c(1,2)]
rownames(Olink_Proteomics) <- metadata.COVID$Unique.ID
Olink_Proteomics[1:10,1:10]

# No need to rescale with the standard deviations, this is Olink data
# Olink_Proteomics <- apply(Olink_Proteomics, 2, function(x){return((x-mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))})
# Olink_Proteomics <- sweep(as.matrix(Olink_Proteomics[, names(sds.longitudinal)[names(sds.longitudinal) %in% colnames(Olink_Proteomics)]]), MARGIN=2, sds.longitudinal[names(sds.longitudinal) %in% colnames(Olink_Proteomics)], `*`)

# Only filter out bad particpants
# When predicting, little data is better than no data
# Only filter out bad participants
Olink_Proteomics.COVID <- filter_olink_data(Olink_Proteomics, participant_cutoff = 0.49948682860075, protein_cutoff = 1)
Olink_Proteomics.COVID <- Olink_Proteomics.COVID[, colnames(Olink_Proteomics.COVID) %in% longitudinal.proteins]
Olink_Proteomics.COVID <- Olink_Proteomics.COVID[, sort(colnames(Olink_Proteomics.COVID))]
Olink_Proteomics.COVID <- impute::impute.knn(as.matrix(Olink_Proteomics.COVID), k = 10, rowmax = 0.5, colmax = 0.5, maxp = 1500, rng.seed = 1714933057)$data
Olink_Proteomics.COVID[1:10, 1:10]
dim(Olink_Proteomics.COVID)
# 784 1400 # No participants are filtered out

# saveRDS(Olink_Proteomics.COVID, file = paste0(rds.dir, "Olink_Proteomics_COVID.rds"))
Olink_Proteomics.COVID <- readRDS(file = paste0(rds.dir, "Olink_Proteomics_COVID.rds"))

df.COVID.list <- vector(mode = "list", length = length(coefficients.longitudinal$gen2))
names(df.COVID.list) <- names(coefficients.longitudinal$gen2)

coefficients.longitudinal.gen2.COVID <- vector(mode = "list", length = length(coefficients.longitudinal$gen2))
names(coefficients.longitudinal.gen2.COVID) <- names(coefficients.longitudinal$gen2)

for(k in 1:length(coefficients.longitudinal$gen2)){
  
  if(!(names(coefficients.longitudinal$gen2)[k] %in% c("Bladder", "Thyroid"))){
    
    coef.tmp <- coefficients.longitudinal$gen2[[k]]
    coefficients.longitudinal.gen2.COVID[[k]] <- coef.tmp[names(coef.tmp) %in% colnames(Olink_Proteomics.COVID), drop = FALSE]
    coef.tmp <- coefficients.longitudinal.gen2.COVID[[k]][coefficients.longitudinal.gen2.COVID[[k]] != 0]
    tmp2 <- Olink_Proteomics.COVID[, names(coef.tmp), drop = FALSE]
    
    df.predictions <- data.frame(
      Unique.ID = rownames(tmp2),
      predicted.ages = rowSums(sweep(tmp2, 2, coef.tmp, FUN="*"))
    )
    
    df.predictions <- dplyr::left_join(metadata.COVID, df.predictions, by = c("Unique.ID" = "Unique.ID"))
    
    df.predictions$age.numeric <- NA
    df.predictions$age.numeric[df.predictions$Age.cat == 1] <- median(c(20, 34)) # 27
    df.predictions$age.numeric[df.predictions$Age.cat == 2] <- median(c(36, 49)) # 42.5
    df.predictions$age.numeric[df.predictions$Age.cat == 3] <- median(c(50, 64)) # 57
    df.predictions$age.numeric[df.predictions$Age.cat == 4] <- median(c(65, 79)) # 72
    df.predictions$age.numeric[df.predictions$Age.cat == 5] <- 85
    
    # Add the residual ages, based on age categories
    # Not perfect, there will still be age information in here!!!
    df.predictions$resid.ages <- resid(lm(predicted.ages ~ age.numeric, data = df.predictions))
    
    df.COVID.list[[k]] <- df.predictions
    
  }
  
}

# saveRDS(coefficients.longitudinal.gen2.COVID, file = paste0(rds.dir, "coefficients_longitudinal_gen2_COVID.rds"))
coefficients.longitudinal.gen2.COVID <- readRDS(file = paste0(rds.dir, "coefficients_longitudinal_gen2_COVID.rds"))

# saveRDS(df.COVID.list, file = paste0(rds.dir, "df_COVID_list.rds"))
df.COVID.list <- readRDS(file = paste0(rds.dir, "df_COVID_list.rds"))


### Plots ###

violin.plots.COVID <- vector(mode = "list", length = 3)
names(violin.plots.COVID) <- c("Age.cat", "COVID", "Acuity.max")

for(m in 1:length(violin.plots.COVID)){
  violin.plots.COVID[[m]] <- vector(mode = "list", length = length(coefficients.longitudinal.gen2.COVID))
  names(violin.plots.COVID[[m]]) <- names(coefficients.longitudinal.gen2.COVID)
}

for(m in 1:length(violin.plots.COVID)){
  
  predictor <- names(violin.plots.COVID)[m]
  
  if(predictor == "Age.cat"){
    colors <- c("#70ad00", "#5Baa5C", "#4585b7", "#aa521c", "#80521c")
    conversion <- 1:5
    names(conversion) <- c("20-34", "36-49", "50-64", "65-79", "80+")
    outcome <- "predicted.ages"
    xlab <- "Age category"
    ylab <- "Predicted rel. log(mortality hazard)"
  } else if(predictor == "COVID"){
    colors <- c("#4189ee", "#aa410e")
    conversion <- c(1,2)
    names(conversion) <- c("Neg.", "Pos.")
    outcome <- "resid.ages"
    xlab <- "COVID status"
    ylab <- "Log(mortality hazard) deviation"
  } else if(predictor == "Acuity.max"){
    colors <- c("#133bff", "#0086a8", "#f6c200", "#d04e00", "#a00e00")
    conversion <- 1:5
    names(conversion) <- c("Discharged", "Hospitalized - O2", "Hospitalized + O2", "Intubated/ventilated", "Death")
    outcome <- "resid.ages"
    xlab <- "Acuity"
    ylab <- "Log(mortality hazard) deviation"
  }
  
  for(k in 1:length(violin.plots.COVID[[m]])){
    if(!(names(violin.plots.COVID[[m]])[k] %in% c("Bladder", "Thyroid"))){
      
      df.0 <- df.COVID.list[[k]][df.COVID.list[[k]]$Day == 0,]
      df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
      df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
      df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
      df.0$outcome <- df.0[, outcome]
      
      if(predictor == "COVID"){
        height <- (max(df.0$resid.ages, na.rm = TRUE) - min(df.0$resid.ages, na.rm = TRUE))*0.9+min(df.0$resid.ages, na.rm = TRUE)
        significance.df <- data.frame(x =1.5, 
                                      y = height, 2, 
                                      pval = summary(lm(predicted.ages ~ COVID + age.numeric, data = df.0))$coefficients["COVID", "Pr(>|t|)"],
                                      stars = "n.s.",
                                      size = 4.5,
                                      nudge_y = (max(df.0$resid.ages, na.rm = TRUE) - min(df.0$resid.ages, na.rm = TRUE))/13
        )
        significance.df$stars[which(significance.df$pval < 0.1)] <- "'"
        significance.df$stars[which(significance.df$pval < 0.05)] <- "*"
        significance.df$stars[which(significance.df$pval < 0.01)] <- "**"
        significance.df$stars[which(significance.df$pval < 0.001)] <- "***"
        
        significance.df[(significance.df$stars != "n.s."), "size"] <- 10
        significance.df[(significance.df$stars != "n.s."), "nudge_y"] <- (max(df.0$resid.ages, na.rm = TRUE) - min(df.0$resid.ages, na.rm = TRUE))/30
      }
      
      violin.plots.COVID[[m]][[k]] <- ggplot_gtable(ggplot_build(
        ggplot(df.0, aes(x = predictor, y = outcome)) +
          
          geom_quasirandom(cex = 1.5, aes(x = predictor, color = predictor, alpha = 1)) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          
          geom_boxplot(width = 0.5, alpha = 0, outlier.shape = NA) + # fill = Diet, alpha = 1
          guides(alpha = "none") + # no legend for "alpha"
          
          {if(predictor == "COVID") geom_text(aes(x, y, label = stars), size = significance.df$size, nudge_y = significance.df$nudge_y, data = significance.df, na.rm = TRUE, colour = "black")} +
          {if(predictor == "COVID") geom_segment(aes(x = 1.1, y = height, xend = 1.9, yend = height))} +
          
          ggtitle(paste0(gsub("_", "-", paste(toupper(substr(names(violin.plots.COVID[[m]])[k], 1, 1)), substr(names(violin.plots.COVID[[m]])[k], 2, nchar(names(violin.plots.COVID[[m]])[k])), sep="")), " model")) +
          xlab(xlab) + # No label on the horizontal axis
          ylab(ylab) +
          theme_bw() +
          theme(legend.title = element_text(size = 12), # legend title size
                legend.text = element_text(size = 12), # legend text size
                legend.position = "none", # No legend
                panel.border = element_blank(),
                panel.background = element_blank(),
                plot.title = element_text(size = 12),
                strip.text.x = element_text(size = 12), # the size of the facet labels
                axis.title.x = element_text(size = 12, color = "black"), # grey color: "#808080"
                axis.title.y = element_text(size = 12, color = "black"),
                axis.text.x = element_text(size = 12, color = "black"),
                axis.text.y = element_text(size = 12, color = "black"),
                axis.line.x = element_line(color="black", linewidth = 0.5),
                axis.line.y = element_line(color="black", linewidth = 0.5),
                axis.ticks.x = element_blank(), # element_line(color="black", linewidth = 0.5),
                axis.ticks.y = element_line(color="black", linewidth = 0.5),
                strip.background = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())
      ))
      tmp <- df.0[, c("predictor", "outcome")]
      colnames(tmp) <- c(xlab, ylab)
      attr(violin.plots.COVID[[m]][[k]], "data") <- tmp
      
    }
  }
}

### Get dimensions ###
dim(df.COVID.list[[k]][df.COVID.list[[k]]$Day == 0,])
# 383  48

table(df.COVID.list[[k]][df.COVID.list[[k]]$Day == 0, "COVID"])
# 0   1 
# 78 305

### Statistics COVID status ###

predictor <- "COVID"
conversion <- c(1,2)
names(conversion) <- c("Neg.", "Pos.")
outcome <- "resid.ages"
df.0 <- df.COVID.list[["Conventional"]][df.COVID.list[["Conventional"]]$Day == 0,]
df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
df.0$outcome <- df.0[, outcome]

summary(lm(predicted.ages ~ COVID + age.numeric, data = df.0))
summary(lm(resid.ages ~ COVID, data = df.0))
t.test(resid.ages ~ COVID, data = df.0, var.equal = TRUE)
# p = 0.004 **

df.0 <- df.COVID.list[["Lung"]][df.COVID.list[["Lung"]]$Day == 0,]
df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
df.0$outcome <- df.0[, outcome]

summary(lm(predicted.ages ~ COVID + age.numeric, data = df.0))
summary(lm(resid.ages ~ COVID, data = df.0))
t.test(resid.ages ~ COVID, data = df.0, var.equal = TRUE)
# p = 0.04 *

### Statistics association with age ###

predictor <- "Age.cat"
conversion <- 1:5
names(conversion) <- c("20-34", "36-49", "50-64", "65-79", "80+")
outcome <- "predicted.ages"
df.0 <- df.COVID.list[["Conventional"]][df.COVID.list[["Conventional"]]$Day == 0,]
df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
df.0$outcome <- df.0[, outcome]

summary(lm(predicted.ages ~ age.numeric, data = df.0))
# Call:
#   lm(formula = predicted.ages ~ age.numeric, data = df.0)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.1595 -1.1987 -0.2751  1.0237  4.8121 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.190693   0.296850  -4.011 7.28e-05 ***
#   age.numeric  0.091626   0.004746  19.306  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.626 on 381 degrees of freedom
# Multiple R-squared:  0.4945,	Adjusted R-squared:  0.4932 
# F-statistic: 372.7 on 1 and 381 DF,  p-value: < 2.2e-16

df.0 <- df.COVID.list[["Lung"]][df.COVID.list[["Lung"]]$Day == 0,]
df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
df.0$outcome <- df.0[, outcome]

summary(lm(predicted.ages ~ age.numeric, data = df.0))
# Call:
#   lm(formula = predicted.ages ~ age.numeric, data = df.0)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.12416 -0.43780 -0.04191  0.43989  2.47248 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 3.879677   0.124990   31.04  < 2e-16 ***
#   age.numeric 0.011370   0.001998    5.69 2.54e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6848 on 381 degrees of freedom
# Multiple R-squared:  0.07832,	Adjusted R-squared:  0.0759 
# F-statistic: 32.37 on 1 and 381 DF,  p-value: 2.539e-08

### Statistics association with acuity ###

predictor <- "Acuity.max"
conversion <- 1:5
names(conversion) <- c("Discharged", "Hospitalized - O2", "Hospitalized + O2", "Intubated/ventilated", "Death")
outcome <- "resid.ages"
pvals <- rep(NA, length(organ.proteins.selected))
names(pvals) <- names(organ.proteins.selected)

for(k in 1:length(pvals)){
  
  df.0 <- df.COVID.list[[names(organ.proteins.selected)[k]]][df.COVID.list[[names(organ.proteins.selected)[k]]]$Day == 0,]
  df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
  df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
  df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
  df.0$outcome <- df.0[, outcome]
  
  pvals[k] <- summary(lm(resid.ages ~ Acuity.max, data = df.0))$coefficients["Acuity.max", "Pr(>|t|)"]
  
}

signif(p.adjust(pvals, method = "hommel"), digits = 1)
# Conventional  Brain       Artery        Liver    Intestine       Immune       Kidney         Skin         Lung 
# 3e-13         3e-11        3e-10        4e-09        1e-04        4e-12        8e-06        1e-03        2e-01 

### r-value for Lung ###
k <- 9
df.0 <- df.COVID.list[[names(organ.proteins.selected)[k]]][df.COVID.list[[names(organ.proteins.selected)[k]]]$Day == 0,]
df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
df.0$outcome <- df.0[, outcome]
summary(lm(resid.ages ~ Acuity.max, data = df.0))
# 0.003

### r-value for Conventional ###
k <- 1
df.0 <- df.COVID.list[[names(organ.proteins.selected)[k]]][df.COVID.list[[names(organ.proteins.selected)[k]]]$Day == 0,]
df.0$COVID <- df.0$COVID+1 # Needs to be in 1,2 format, not 0,1
df.0$Acuity.max <- 6-df.0$Acuity.max # Needs to go from 1 to 5, not 5 to 1
df.0$predictor <- factor(names(conversion[df.0[, predictor]]), levels = names(conversion))
df.0$outcome <- df.0[, outcome]
summary(lm(resid.ages ~ Acuity.max, data = df.0))
# 0.14
