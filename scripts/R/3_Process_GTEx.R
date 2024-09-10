
### This comes from 1_Prepare_Data.R ###
olink_data_0 <- readRDS(file = paste0(rds.dir, "olink_data_0.rds"))
coding143 <- readRDS(file = paste0(rds.dir, "coding143.rds"))

### This comes from 2_Prepare_Longitudinal_Data.R ###
longitudinal.proteins <- readRDS(file = paste0(rds.dir, "longitudinal_proteins.rds"))

gene_reads <- as.data.frame(data.table::fread(paste0(input.GTEx.dir, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")))
gene_reads[1:10, 1:10]

SampleAttributes <- as.data.frame(data.table::fread(paste0(input.GTEx.dir, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")))
head(SampleAttributes)
unique(SampleAttributes$SMTSD)
rownames(SampleAttributes) <- SampleAttributes$SAMPID

# A data frame with the average gene expression counts per tissue
tissues <- unique(SampleAttributes$SMTSD)
tissue.list <- vector(mode = "list", length = length(tissues))
names(tissue.list) <- tissues

for(j in 1:length(tissue.list)){
  column.names <- SampleAttributes$SAMPID[SampleAttributes$SMTSD == tissues[j]]
  tissue.list[[j]] <- which(colnames(gene_reads) %in% column.names)
}

tissue.list <- tissue.list[which(unlist(lapply(tissue.list, function(x){length(x) != 0})))]
tissues <- tissues[which(unlist(lapply(tissue.list, function(x){length(x) != 0})))]

avg.counts <- cbind(gene_reads[, c(1, 2)], as.data.frame(matrix(NA, nrow = nrow(gene_reads), ncol = length(tissues))))
colnames(avg.counts) <- c(colnames(gene_reads)[c(1, 2)], tissues)

for(j in 1:length(tissue.list)){
  avg.counts[, (j+2)] <- rowMeans(gene_reads[, tissue.list[[j]]])
}

### Now group per tissue ###
GTEx.grouping <- openxlsx::read.xlsx(paste0(input.GTEx.dir, "GTEx_grouping.xlsx"))

organs <- na.omit(unique(GTEx.grouping$organ))
organ.counts <- cbind(gene_reads[, c(2)], as.data.frame(matrix(NA, nrow = nrow(gene_reads), ncol = length(organs))))
colnames(organ.counts) <- c(colnames(gene_reads)[c(2)], organs)

organ.list <- vector(mode = "list", length = length(organs))
names(organ.list) <- organs

for(j in 1:length(organ.list)){
  column.names <- GTEx.grouping[which(GTEx.grouping$organ == organs[j]), "tissue"]
  organ.list[[j]] <- which(make.names(colnames(avg.counts)) %in% column.names)
  if(!all(column.names %in% make.names(colnames(avg.counts)))){stop("Some GTEx tissues are missing!")}
}

for(j in 1:length(organ.list)){
  organ.counts[, (j+1)] <- apply(avg.counts[, organ.list[[j]], drop = FALSE], 1, max)
}
dim(organ.counts)
# 56200    22
# Let's aggregate duplicate genes with different ENSG numbers:
organ.counts <- organ.counts %>% group_by(Description) %>% summarize_all(.funs = sum)
organ.counts <- as.data.frame(organ.counts)
dim(organ.counts)
# 54592    22

maxima <- apply(organ.counts[, -c(1)], 1, which.max)
sort(unique(maxima))

n <- ncol(organ.counts[, -c(1)])
# This is the second largest value
second.largest <- apply(organ.counts[, -c(1)], 1, function(x){
  x <- unlist(x)
  # return(which(x == sort(x, partial = n-1)[n-1])[1])
  return(sort(x, partial=n-1)[n-1][1])
})

# Then determine organ enrichment and make plots for number of proteins per organ

enriched.proteins <- vector(mode = "list", length = length(organs))
names(enriched.proteins) <- organs

organ.counts[41,]

for(j in 1:length(enriched.proteins)){
  
  first <- organ.counts[which(maxima == j),(j+1)]
  second <- second.largest[which(maxima == j)]
  
  enriched.proteins[[j]] <- organ.counts[which(maxima == j),"Description"][which((first >= 4*second) & (first != 0))]
  
  if(!all(first >= second)){stop("Error: there are maxima that are smaller than the second largest value: something went wrong before.")}
}

### Import the Uniprot data to identify proteins that are missing based on their Ensembl IDs: ###
Uniprot <- data.table::fread(paste0(input.uniprot.dir, uniprot.file))
head(Uniprot)

Uniprot$transcripts <- gsub("\\..*", "", Uniprot$Ensembl)

indices <- which(Uniprot$transcripts != "")

# ### Check ###
# test <- grepl("^ENST", Uniprot$transcripts[indices])
# all(test)
# # TRUE

# define biomart object
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org/")
# query biomart
results <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"),
                 filters = "ensembl_transcript_id", values = Uniprot$transcripts[indices],
                 mart = mart)
head(results)

Uniprot <- dplyr::left_join(Uniprot, results, by = c("transcripts" = "ensembl_transcript_id"))

# saveRDS(Uniprot, file = paste0(rds.dir, "Uniprot.rds"))
Uniprot <- readRDS(file = paste0(rds.dir, "Uniprot.rds"))

gene.names <- gene_reads[, c(1, 2)]
gene.names$ensembl_gene_id <- gsub("\\..*", "", gene.names$Name)
gene.names <- gene.names[, c("Description", "ensembl_gene_id")]
gene.names <- gene.names[!duplicated(gene.names), ]
gene.names$ensembl_gene_id[which(duplicated(gene.names$ensembl_gene_id))]

Uniprot <- Uniprot[!is.na(Uniprot$ensembl_gene_id),]
Uniprot <- Uniprot[, c("Gene Names", "ensembl_gene_id")]
Uniprot <- Uniprot[!duplicated(Uniprot), ]
Uniprot$ensembl_gene_id[which(duplicated(Uniprot$ensembl_gene_id))]
Uniprot$`Gene Names` <- unlist(lapply(strsplit(Uniprot$`Gene Names`, " "), function(x){return(x[1])})) # The Uniprot gene name
Uniprot <- Uniprot[!duplicated(Uniprot), ]
Uniprot$ensembl_gene_id[which(duplicated(Uniprot$ensembl_gene_id))] # No duplicates anymore

gene.names <- dplyr::left_join(gene.names, Uniprot)

head(gene.names)
gene.names$Olink.name <- gene.names$`Gene Names`
gene.names$Olink.name[is.na(gene.names$Olink.name)] <- gene.names$Description[is.na(gene.names$Olink.name)]

# Need to fix this, otherwise 2x HSFY1:
# gene.names[gene.names$Olink.name == "HSFY1", ]
# Uniprot[Uniprot$ensembl_gene_id %in% c("ENSG00000172468", "ENSG00000169953"),]
gene.names[gene.names$Description == "HSFY2", "Olink.name"] <- "HSFY2"

sum(coding143$protein %in% gene.names$Description)
sum(coding143$protein %in% gene.names$Olink.name)

coding143$protein[!(coding143$protein %in% gene.names$Olink.name)]

# The missing Olink names, need to do them manually:
Olink.names <- coding143$protein[!(coding143$protein %in% gene.names$Olink.name)]
# Cannot find KIR2DL2, nor its synonyms, nor its ENSG number (gene.names[gene.names$ensembl_gene_id == "ENSG00000273661", ])
# Same for LILRA3 (ENSG00000275841)
# NTproBNP is the N-terminus of the B-type natriuretic peptide (BNP)
GTEx.names <- c("AMY1A", "ARNTL", "BOLA2", "COL4A3BP", "CGB3",
                "CKMT1A", "CTAG1A", "DDX58", "DEFA1", "DEFB103A",
                "DEFB104A", "DEFB4A", "EBI3", "FAM172A", "FUT3", "GBA", 
                "GP1BB", "C10orf99", "ICAM4", "IL12A", NA, 
                "LGALS7", NA, "MICB", "DPCR1", "MYLPF", 
                NA, "SKIV2L", "SLC9A3R1", "SLC9A3R2", "SPACA5", 
                "TDGF1", "TEX33", "WARS")

# "GP1BB" and "ICAM4" are in GTEx, but have no bulk tissue expression available:
GTEx.names[!(GTEx.names %in% gene_reads$Description)]
GTEx.names[GTEx.names %in% c("GP1BB", "ICAM4")] <- NA

indices <- which(gene.names$Description %in% GTEx.names[!is.na(GTEx.names)])
df <- gene.names[indices,]
dim(df)
# 29 4
# To sort back again:
ensembl_gene_ids <- df$ensembl_gene_id
# Sort by GTEx.names
df <- df[match(GTEx.names[!is.na(GTEx.names)], df$Description),]
# Then put the Olink names for those GTEx names to the Olink.names object
# Description is the GTEx identifier
# Olink.name is the Olink identifier
df$Olink.name <- Olink.names[!is.na(GTEx.names)]

# Sort to the original order and put those back
df <- df[match(ensembl_gene_ids, df$ensembl_gene_id),]
gene.names[indices,] <- df

coding143$protein[!(coding143$protein %in% gene.names$Olink.name)]
# "GP1BB"    "ICAM4"    "KIR2DL2"  "LILRA3"   "NTproBNP" # as discussed above

enriched.proteins[[1]] %in% gene.names$Description
unlist(lapply(enriched.proteins, function(x){return(all(x %in% gene.names$Description))}))

### Now convert enriched.proteins to Olink IDs ###

subsets.Olink <- vector(mode = "list", length = length(enriched.proteins))
names(subsets.Olink) <- names(enriched.proteins)

for(i in 1:length(subsets.Olink)){
  subsets.Olink[[i]] <- unique(gene.names[gene.names$Description %in% enriched.proteins[[i]],"Olink.name"])
}

# saveRDS(subsets.Olink, file = paste0(rds.dir, "subsets_Olink.rds"))
subsets.Olink <- readRDS(file = paste0(rds.dir, "subsets_Olink.rds"))

sum(subsets.Olink[[9]] %in% coding143$protein)

# Tony Wyss-Corrat has 856 proteins organ-enriched, 12 organs in total; ~5000 proteins in total; means ~70 proteins per organ if equally distributed
# Some seem only fitted on 3 (!) proteins, see ST16!
n.proteins <- unlist(lapply(subsets.Olink, function(x){return(sum(x %in% coding143$protein))}))
n.proteins

### 1. Get the set of organ-enriched proteins per organ for the non-longitudinal models ###

proteins <- coding143[match(colnames(olink_data_0), coding143$protein), "protein"]
organ.proteins <- lapply(subsets.Olink, function(x){
  x <- sort(x[x %in% coding143$protein]) # Sort alphabetically; otherwise, we have the random order from GTEx
  x <- x[x %in% proteins]
  return(x)
})
organ.proteins$Organismal <- proteins[!(proteins %in% unlist(subsets.Olink)[unlist(subsets.Olink) %in% coding143$protein])]
organ.proteins$`Multi-organ` <- proteins[proteins %in% unlist(subsets.Olink)[unlist(subsets.Olink) %in% coding143$protein]]
organ.proteins$Conventional <- proteins

# saveRDS(organ.proteins, paste0(rds.dir, "organ_proteins.rds"))
organ.proteins <- readRDS(paste0(rds.dir, "organ_proteins.rds"))

ListJSON <- toJSON(organ.proteins, pretty = TRUE, auto_unbox = TRUE)
write(ListJSON, paste0(dir.4x.FC.genes, "GTEx_4x_FC_genes.json"))

# After this, run the Python scripts
# Where you will find the output:
# dir.gen1.models
# dir.gen2.models

### 2. Get the set of organ-enriched proteins per organ for the longitudinal models ###

proteins <- coding143[match(colnames(olink_data_list[["instance_0"]]), coding143$protein), "protein"]
organ.proteins.longitudinal <- lapply(subsets.Olink, function(x){
  x <- sort(x[x %in% coding143$protein]) # Sort alphabetically; otherwise, we have the random order from GTEx
  x <- x[x %in% proteins]
  return(x)
})
organ.proteins.longitudinal$Organismal <- proteins[!(proteins %in% unlist(subsets.Olink)[unlist(subsets.Olink) %in% coding143$protein])]
organ.proteins.longitudinal$`Multi-organ` <- proteins[proteins %in% unlist(subsets.Olink)[unlist(subsets.Olink) %in% coding143$protein]]
organ.proteins.longitudinal$Conventional <- proteins

# saveRDS(organ.proteins.longitudinal, paste0(rds.dir, "organ_proteins_longitudinal.rds"))
organ.proteins.longitudinal <- readRDS(paste0(rds.dir, "organ_proteins_longitudinal.rds"))

ListJSON.longitudinal <- toJSON(organ.proteins.longitudinal, pretty = TRUE, auto_unbox = TRUE)
write(ListJSON.longitudinal, paste0(dir.4x.FC.genes, "GTEx_4x_FC_genes_feature_reduced.json"))

# After this, run the Python scripts
# Where you will find the output:
# dir.gen1.models.longitudinal
# dir.gen2.models.longitudinal

