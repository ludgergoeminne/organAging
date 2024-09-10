
### This comes from Plot_Barplots_Correlations.R ###
organ.proteins.selected <- readRDS(file = paste0(rds.dir, "organ_proteins_selected.rds"))

kfolds_test <- lapply(1:5, function(x) read.csv(paste0(dir.Python.input, "test_imputed_",x, ".csv"), row.names = 1))

organs <- names(organ.proteins.selected)
organs <- organs[!organs %in% c("Conventional")]

aft_predictions <- vector(mode = "list", length = length(organs))
names(aft_predictions) <- organs

for (organ in organs){
  aft_predictions[[organ]] <- c()
  for (kf in 1:5){
    aft_coefs <- read.csv(paste0(dir.AFT.models, organ, "_aft_", kf, "_coefs_GTEx_4x_FC.csv"))
    aft_coefs <- aft_coefs[, colnames(aft_coefs)!= "Intercept", drop=FALSE]
    data_prediction_subset <- kfolds_test[[kf]][, colnames(aft_coefs)]
    predicted <- (as.matrix(data_prediction_subset) %*% t(as.matrix(aft_coefs)))*(-1)
    aft_cor <- cor(predicted, kfolds_test[[kf]][["Age"]])
    aft_predictions[[organ]] <- append(aft_predictions[[organ]], aft_cor)
  }
}

aft_correlations_median <- lapply(aft_predictions, median)

elastic_predictions <- vector(mode = "list", length = length(organs))
names(elastic_predictions) <- organs

for (organ in organs){
  elastic_predictions[[organ]] <- c()
  coefs <- read.csv(paste0(dir.gen2.models, organ, "_mortality_coefs_GTEx_4x_FC.csv"))
  for (kf in 1:5){
    data_prediction_subset <- kfolds_test[[kf]][, colnames(coefs)]
    predicted <- (as.matrix(data_prediction_subset) %*% t(as.matrix(coefs[kf,])))
    elastic_cor <- cor(predicted, kfolds_test[[kf]][["Age"]])
    elastic_predictions[[organ]] <- append(elastic_predictions[[organ]], elastic_cor)
  }
}

elastic_correlations_median <- lapply(elastic_predictions, median)

plot.df <- data.frame("Organ" = c(names(aft_correlations_median), names(elastic_correlations_median)), 
                   "Model" = c(rep("Accelerated Failure Time", length(aft_correlations_median)), rep("Cox Elastic Net", length(aft_correlations_median))), 
                   "Correlation" = c(unlist(aft_correlations_median), unlist(elastic_correlations_median)))

plot.df[["Correlation"]] <- as.numeric(plot.df[["Correlation"]])
plot.df <- as.data.frame(plot.df)

barplot.other.models <- ggplot(plot.df, aes(fill = as.factor(Model), x = Correlation, y = reorder(Organ, Correlation))) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("burlywood1", "salmon")) +
  xlab("Correlation with chronological age") +
  ylab("") +
  labs(fill = "Model")+
  theme_bw() +
  theme(legend.position = c(0.85, 0.5), # The coordinates for legend.position are x- and y- offsets from the bottom-left of the plot, ranging from 0 - 1
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        # plot.title = element_text(face = "bold"), # element_blank(),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"), # , angle = 45, hjust = 1, vjust = 1
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        axis.ticks.x = element_line(color="black", linewidth = 0.5),
        axis.ticks.y = element_line(color="black", linewidth = 0.5),
        strip.text = element_text(face = "bold",size = 12),
        strip.background = element_blank())

attr(barplot.other.models, "data") <- plot.df

plot(barplot.other.models)




