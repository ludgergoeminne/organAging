# This function performs functional GSEA with KEGG, REACTOME and GO BP ontology simultaneously for
# multiple comparisons (using logFC and p-value matrix as an input)

Run_functional_GSEA_from_matrix <- function(logFC_table,pvalue_table=NULL,
                                            species="Mus musculus",convert_to_symbol=T,
                                            n_permutations=1000,
                                            min_size=15,max_size=500,
                                            with_MF=F,with_hallmarks=F,with_ChEA=F,
                                            gseaParam=1,multilevel=T){

  if (!is.null(pvalue_table)){
    #A) Align logFC and pvalue tables
    common_names <- intersect(colnames(logFC_table),colnames(pvalue_table))
    common_genes <- intersect(rownames(logFC_table),rownames(pvalue_table))
    logFC_table <- logFC_table[common_genes,common_names]
    pvalue_table <- pvalue_table[common_genes,common_names]
  }

  #B) Create list of ranked sets
  ranked_list <- list()
  for (i in colnames(logFC_table)){
    temp_name <- i
    print(temp_name)
    if (!is.null(pvalue_table)){
      temp_signature <- cbind(logFC_table[,i,drop=F],pvalue_table[,i,drop=F])
      rownames(temp_signature) <- rownames(logFC_table)
      temp_signature <- temp_signature[complete.cases(temp_signature),]
      colnames(temp_signature) <- c("logFC","PValue")
      temp_signature <- data.frame(temp_signature)
      temp_signature <- temp_signature[temp_signature$logFC!="NA",]
      
      temp_ranked_list <- Ranked_list_obtain(pvalue_table = as.numeric(as.character(temp_signature$PValue)),
                                             logFC_table = as.numeric(as.character(temp_signature$logFC)),
                                             tolower(rownames(temp_signature)),
                                             convert_to_symbol = convert_to_symbol,species=species)
    }else{
      temp_signature <- data.frame(logFC=logFC_table[,i])
      rownames(temp_signature) <- rownames(logFC_table)
      temp_signature <- temp_signature[complete.cases(temp_signature),,F]
      colnames(temp_signature) <- c("logFC")
      temp_signature <- temp_signature[temp_signature$logFC!="NA",,F]

      temp_ranked_list <- Ranked_list_obtain(pvalue_table = NULL,
                                             logFC_table = as.numeric(as.character(temp_signature$logFC)),
                                             tolower(rownames(temp_signature)),
                                             convert_to_symbol = convert_to_symbol,species=species)
      
    }
    ranked_list[[temp_name]] <- temp_ranked_list
    print(dim(temp_ranked_list))
    print(head(temp_ranked_list))
  }
  list_ranks <- ranked_list
  print("Ranked lists constructed")

  
  #C) Prepare list of functions
  m_df_bp = msigdbr(species = species, category = "C5",subcategory = "BP")
  m_df_kegg = msigdbr(species = species, category = "C2",subcategory = "CP:KEGG")
  m_df_reactome = msigdbr(species = species, category = "C2",subcategory = "CP:REACTOME")
  m_df <- rbind(rbind(m_df_bp,m_df_kegg),m_df_reactome)
  if (with_MF==T){
    m_df_mf <- msigdbr(species = species, category = "C5",subcategory = "MF")
    m_df <- rbind(m_df,m_df_mf)
  }
  if (with_hallmarks==T){
    m_df_hall <- msigdbr(species = species, category = "H")
    m_df <- rbind(m_df,m_df_hall)
  }
  m_df = m_df[,c("gs_cat","gs_subcat","gs_name","gene_symbol")]

  if (with_ChEA==T){
    m_df_chea = enrichr_gsets("ChEA_2022", db="Enrichr")
    if (species=="Mus musculus"){
      m_df_chea = m_df_chea$genesets[grepl("Mouse",names(m_df_chea$genesets))]
    }else if(species=="Homo sapiens"){
      m_df_chea = m_df_chea$genesets[grepl("Human",names(m_df_chea$genesets))]
    }else if(species=="Rattus norvegicus"){
      m_df_chea = m_df_chea$genesets[grepl("Rat",names(m_df_chea$genesets))]
    }
    for (temp_factor in names(m_df_chea)){
      temp_name = paste0(strsplit(temp_factor," ")[[1]][1],"_",
                         strsplit(temp_factor," ")[[1]][4],"_",
                         gsub("-","",strsplit(temp_factor," ")[[1]][3]))
      temp_name = gsub("-","",temp_name)
      temp_list = m_df_chea[[temp_factor]]
      temp_dataset = data.frame(gs_cat="ChEA",
                                gs_subcat=paste0("ChEA_2022.",strsplit(temp_factor," ")[[1]][3]),
                                gs_name=paste0("ChEA_",temp_name),gene_symbol=temp_list)
      m_df = rbind(m_df,temp_dataset)
    }
  }

  functions_list <- levels(factor(m_df$gs_name))
  functions_annotation <- list()
  for (i in functions_list){
    temp_data <- m_df[m_df$gs_name==i,]
    temp_data <- toupper(as.character(temp_data$gene_symbol))
    functions_annotation[[i]] <- temp_data
  }
  print("Functional lists constructed")
  
  #D) Run GSEA
  print("GSEA is running:")
  time_0 <- Sys.time()
  gsea_output <- list()
  for (i in names(list_ranks)){
    print(i)
    
    ranked_list <- list_ranks[[i]]
    ranked_list[,1] <- toupper(as.character(ranked_list[,1]))
    ranked_vector <- as.numeric(as.character(ranked_list[,2]))
    names(ranked_vector) <- as.character(ranked_list[,1])

    if (multilevel==T){
      temp <- fgseaMultilevel(pathways = functions_annotation, 
                              stats    = ranked_vector,
                              minSize  = min_size,
                              maxSize  = max_size,
                              gseaParam=gseaParam)
    }else{
      temp <- fgsea(pathways = functions_annotation, 
                        stats    = ranked_vector,
                        nperm = n_permutations,
                        minSize  = min_size,
                        maxSize  = max_size,
                        gseaParam=gseaParam)
    }

    gsea_output[[i]] <- temp
  }
  print("GSEA is completed")
  print(Sys.time()-time_0)

    
  #E) Make single table
  functions_all <- c()
  for (i in names(gsea_output)){
    functions_all <- union(functions_all,as.character(gsea_output[[i]]$pathway))
  }
  function_table <- matrix(nrow=length(functions_all),ncol=length(gsea_output))
  rownames(function_table) <- functions_all
  colnames(function_table) <- names(gsea_output)
  function_table_NES <- function_table
  function_table_pvalue <- function_table
  for (i in names(gsea_output)){
    temp_data <- gsea_output[[i]]
    ontologies <- unique(sapply(as.character(temp_data$pathway),function(x){strsplit(x,"_")[[1]][1]}))
    temp_data$P.Adjusted <- 1
    for (temp_ontology in ontologies){
      temp_data[grepl(paste0("^",temp_ontology,"_"),temp_data$pathway),]$P.Adjusted <- as.numeric(p.adjust(temp_data[grepl(paste0("^",temp_ontology,"_"),temp_data$pathway),]$pval,"BH"))
    }
    temp <- -log10(temp_data$P.Adjusted)*sign(temp_data$NES)
    function_table[as.character(temp_data$pathway),i] <- temp
    function_table_NES[as.character(temp_data$pathway),i] <- temp_data$NES
    function_table_pvalue[as.character(temp_data$pathway),i] <- temp_data$pval
  }
  GSEA_list <- list(NES=function_table_NES,Zscore=function_table,
                    P.Value=function_table_pvalue,
                    Complete_output=gsea_output)

  return(GSEA_list)
}



