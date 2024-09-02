# Plasma protein-based organ-specific aging and mortality models unveil diseases as accelerated aging of organismal systems

This repository accompanies our Cell Metabolism manuscript "Plasma protein-based organ-specific aging and mortality models unveil diseases as accelerated aging of organismal systems". The repository contains all required code to reproduce the analyses for our pulication which is under revision at *Cell Metabolism*, but for which a preprint is available on *medRxiv* (see https://doi.org/10.1101/2024.04.08.24305469).

## Overview of the repository

This repository contains the following main folders:

- `scripts`: this folder contains a subfolder `R` and a subfolder `Python`. The `R` folder contains all `R` scripts. These scripts were used for data preprocessing, data analysis and plotting. The `Python` folder contains all `Python` scripts. These scripts were using to train the models.

- `data`: this folder contains all intermediary data files, as long as they can be made publicly available (i.e. they do not contain individual-level UK Biobank data).

## How to reproduce our results

All analyses (except for model training, which was done in Python) were done in R, running in RStudio. If you want to reproduce our analysis, you will need an installation of R (https://www.r-project.org) and possibly RStudio (https://www.rstudio.com/products/rstudio/download/). Next, download this project to your computer.

If you want to reproduce the full analysis based on individual-level data, you will also need to apply for access to the following datasets:
- The UK Biobank (Sudlow *et al.* (2015), apply from https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access)
- Bild *et al.* (2002) (apply from https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001416.v3.p1)
- Dammer *et al.* (2022) (apply from https://www.synapse.org/#!Synapse:syn30549757/files/)

Additionally, the individual-level data from the following datasets can be downloaded freely:
- Filbin *et al.* (2021) (download from https://doi.org/10.17632/nf853r8xsj.2)
- Damsky *et al.* (2022) (download from https://www.ncbi.xyz/geo/query/acc.cgi?acc=GSE169148)

Note that many temporary results (as long as they do not contain individual-level data from the restricted datasets) have been saved in the folder `/data/rds` and are loaded in the scripts by default to speed up the workflow. However, you can always check the code that was used to generate them. Thanks to the saved files, each file is created in such a way that it can be run independently of the other files.

The most important files are numbered so that if you don't want to make use of the saved intermediary results, you can run them in ascending order just like we did to generate the results.

### Reproducing the figures

To reproducing our figures, please look into the file `/scripts/Arrange_Figures.R` (if you want to reproduce a panel from a main figure) and/or `/scripts/Arrange_Supp_Figures.R` (if you want to reproduce a panel from a supplementary figure). For each figure object that is being assembled, the names of the scripts where the figure panels come from are given. For example, if one would like to reproduced Figure 4A, you will find the following comment in the file `Arrange_Figures.R`:

```
### Arrange, plot, and save Fig. 4 ###
# Plot_Hazard_Ratios.R
# Plot_Heatmaps_HR_Diseases.R
# Plot_Scatter_Leave_Proteins_Out.R
```
This example helps you find that the code to reproduce this hazard ratio plot can be found in the file `Plot_Hazard_Ratios.R`.

The Excel files that from the Data S1 file were assembled with the scripts in file `/scripts/Export_Figures_Excel.R`. To find the scripts where the data was generated, we again suggest to look into the files `/scripts/Arrange_Figures.R` and `/scripts/Arrange_Supp_Figures.R`.

## How to apply the organ-specific models to your own dataset

The best approach depends on the type of data you are using.

### 1. Olink Explore 3072 data

When using data obtained with the Olink Explore 3072 platform, we suggest to use our full models.
The coefficients for these models can be found in Supp. Table 1A (for the 1st-generation models, which were trained to predict chronological age) and in Supp. Table 1C (for the mortality-based models, which were trained to predict time to death).
In principle, the coefficients from any of the five folds can be used, but we suggest to always simply use the first fold (unless you work with the same UK Biobank data, in which case you will have some theoretical potential for overfitting, although in practice, the out-of-fold predictions are very similar to the in-fold predictions).
Then simply multiply the coefficients with the corresponding Olink NPX-normalized intensities and sum them up.
For the 1st-generation models, add the intercept to this value.
The result will be the predicted biological age or the predicted relative log(mortality hazard).

### 2. Olink Explore 1536 data

When using data obtained with the less extensive Olink Explore 1536 platform, we suggest to use our feature-reduced models, as these models were specifically trained on the data from this platform.
The coefficients for these models can be found in Supp. Table 6A (for the 1st-generation models) and in Supp. Table 6C (for the mortality-based models).
Then apply the same procedure as described above.

### 3. Other data types

If you are using a different Olink platform, you should assess which of both aforementioned platforms most closely resembles your data.

If you are using protein-level non-Olink data (e.g. SomaLogic data, Alamar data, mass spectrometry-based proteomics data, ...), we suggest to transform the data so that its distribution is more or less symmetrical (e.g. by log-transforming raw intensity values). Then, for each protein, subtract the mean and divide by the standard deviation. Then multiply these values with the standard deviations from the Olink NPX normalized values for either the full models or the feature-reduced models, as given in Supp. Table 3.

### What about missing values?

In our manuscript, we ignored missing values when applying our models to external datasets. 
However, note that the more proteins are missing (especially those with strong non-zero coefficients), the poorer the performance of our models will be.
One could consider imputing missing values, especially when the dataset is rather large, and there are not too many missing values for each proteins.
This may somewhat improve the performance of the models, but would need to be evaluated on a case-by-case basis.
When in doubt, it is probably better not to impute.

## Help

If anything is unclear or does not work, please do not hesitate to contact lgoeminne@bwh.harvard.edu or raise an issue on GitHub.

## References

1. Sudlow, C., Gallacher, J., Allen, N., Beral, V., Burton, P., Danesh, J., Downey, P., Elliott, P., Green, J., Landray, M., et al. (2015). UK Biobank: An Open Access Resource for Identifying the Causes of a Wide Range of Complex Diseases of Middle and Old Age. PLoS Med 12, e1001779. https://doi.org/10.1371/JOURNAL.PMED.1001779
2. Bild, D.E., Bluemke, D.A., Burke, G.L., Detrano, R., Diez Roux, A. V., Folsom, A.R., Greenland, P., Jacobs, D.R., Kronmal, R., Liu, K., et al. (2002). Multi-Ethnic Study of Atherosclerosis: Objectives and Design. Am J Epidemiol 156, 871–881. https://doi.org/10.1093/AJE/KWF113.
3. Dammer, E.B., Ping, L., Duong, D.M., Modeste, E.S., Seyfried, N.T., Lah, J.J., Levey, A.I., and Johnson, E.C.B. (2022). Multi-platform proteomic analysis of Alzheimer’s disease cerebrospinal fluid and plasma reveals network biomarkers associated with proteostasis and the matrisome. Alzheimers Res Ther 14, 1–32. https://doi.org/10.1186/S13195-022-01113-5/FIGURES/10.
4. Filbin, M.R., Mehta, A., Schneider, A.M., Kays, K.R., Guess, J.R., Gentili, M., Fenyves, B.G., Charland, N.C., Gonye, A.L.K., Gushterova, I., et al. (2021). Longitudinal proteomic analysis of severe COVID-19 reveals survival-associated signatures, tissue-specific cell death, and cell-cell interactions. Cell Rep Med 2, 100287. https://doi.org/10.1016/J.XCRM.2021.100287.
5. Damsky, W., Wang, A., Kim, D.J., Young, B.D., Singh, K., Murphy, M.J., Daccache, J., Clark, A., Ayasun, R., Ryu, C., et al. (2022). Inhibition of type 1 immunity with tofacitinib is associated with marked improvement in longstanding sarcoidosis. Nature Communications 2022 13:1 13, 1–17. https://doi.org/10.1038/s41467-022-30615-x.

## Citation

If you make use of the functions in this repository, please refer to: 

Ludger J.E. Goeminne, Alec Eames, Alexander Tyshkovskiy, M. Austin Argentieri, Kejun Ying, Mahdi Moqri, Vadim N. Gladyshev (2024). **Plasma-based organ-specific aging and mortality models unveil diseases as accelerated aging of organismal systems**. *medRxiv*. doi: https://doi.org/10.1101/2024.04.08.24305469.
