# Ranked_list_obtain is a function for creation of ranked lists (input for GSEA)
# based on logFC, FDR and names of corresponding genes

Ranked_list_obtain <- function(logFC_table,pvalue_table=NULL,gene_names,
                               convert_to_symbol=T,species="Mus musculus"){
        
        if (convert_to_symbol){
          if (species=="Mus musculus" | species=="Mus" | species=="mouse" | species=="mus" | species=="Mouse"){
            gene_symbols <- getSYMBOL(as.character(gene_names),
                                      data='org.Mm.eg')[gene_names]
          }else{
            gene_symbols <- getSYMBOL(as.character(gene_names),
                                      data='org.Hs.eg')[gene_names]
          }
        }else{
                gene_symbols <- gene_names
        }
        
        if (!is.null(pvalue_table)){
                ranked_list <- data.frame(Gene_symbol=gene_symbols,
                                          Value=-sign(logFC_table)*
                                                  log10(pvalue_table))
        }else{
                ranked_list <- data.frame(Gene_symbol=gene_symbols,
                                          Value=logFC_table)
        }
        #rownames(ranked_list) <- ranked_list$Gene_symbol
        ranked_list <- ranked_list[complete.cases(ranked_list),]
        ranked_list$Gene_symbol <- toupper(ranked_list$Gene_symbol)
        ranked_list <- ranked_list[order(ranked_list$Value,decreasing = T),]
        ranked_list <- ranked_list[!duplicated(ranked_list$Gene_symbol),]
        return(ranked_list)
}