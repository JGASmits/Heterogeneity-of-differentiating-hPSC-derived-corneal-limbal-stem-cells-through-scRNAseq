GOenrichment_plot_scrna_JS <- function(GO_object, output_file, topn_GOterms, topn_genes, res, seu, plot_width=20,plot_height=10){
#   #generates a dotplot fo a GO_object and visualize the counts 
#   # give a seurat object, a results table of markergenes (of which the most sig genes are picked) 
#   # generating topn_genes contributing to each of the topn_GOTerms
#   #finally, add some metadata to group
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(patchwork)

  GO_df <-  as.data.frame(GO_object)
  str_split(GO_df$BgRatio, '/')
  GO_df <- GO_df %>%
    separate(GeneRatio, c("genes_found", "GO_size"), "/")
  GO_df$genes_found <- as.numeric(GO_df$genes_found)
  #GO_df <- GO_df[order(GO_df$genes_found, decreasing = T),]
  rownames(GO_df) <- NULL

  pdf(output_file ,width=plot_width,height=plot_height,paper='special')
  print(dotplot(GO_object,showCategory = topn_GOterms, color = 'qvalue'))#print overview dotplot
  for (cluster in unique(GO_df$cluster)){
    print(cluster)
    GO_df_subset <- GO_df[(GO_df$cluster == cluster),]
    GO_df_subset$qvalue
    if(dim(GO_df_subset)[1] > topn_GOterms){max_examples <- topn_GOterms}else{max_examples <- dim(GO_df_subset)[1]}
    for (i in 1:max_examples){
     genes_enriched_res <- strsplit(GO_df_subset[i,]$geneID, "/")[[1]]#take the top enrichment gene names
     #sort the enrichment genes based on marker gene expression
     res_subset <- res[res$gene %in% genes_enriched_res,]
     res_subset <- res_subset[res_subset$cluster == cluster,]
     if(dim(res_subset)[1] > topn_genes){max_genes <- topn_genes}else{max_genes <- dim(res_subset)[1]}
     res_subset <- res_subset[order(res_subset$avg_log2FC, decreasing = T),]
     res_subset <- res_subset[1:max_genes,]
     p <- FeaturePlot(seu, features = list(res_subset$gene)[[1]], label = TRUE)
     print(p +   plot_annotation(
       title = paste0(paste0(paste0('enriched go term: ', GO_df_subset[i,]$Description), ' from cluster: '),cluster),
       caption = paste0('cluster: ', cluster),
       theme = theme(plot.title = element_text(size = 16))
     ))
    }}
    dev.off()
}
