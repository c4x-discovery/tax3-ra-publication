# Plot Disease Association Heatmap and corresponding .txt file
# if no genes or gene_list is provided, SNP-to-gene mapping for all Tax3 analyses for the first disease will be used
# gene_list must correspond to a GeneList.listName in Neo4J
# disease_names must correspond to disease ids and first disease name will be used for tissueWeights for Tax3 SNP-to-gene mapping
# include zeros will return all genes in the plot, otherwise genes that score 0 for all associations will not be included
# hide_genes will ensure the gene names are not included on the plot (for anonymity)

library(neo4jshell)
library(heatmap3)
library(WGCNA)

source('neo4j_connection.R')

da_heatmap <- function(disease_ids, disease_names, output_prefix, credentials_neo4j_txt,
                       gene_list=NULL, genes=NULL, include_zeros=TRUE, hide_genes=FALSE,
                       png_width=1000, png_height=1000) {
  con_neo4j <- get_neo4j_connection(credentials_neo4j_txt)
  if (length(disease_ids) != length(disease_names)) {
    return('disease ids and names must match')
  }
  # Get Tax3 Genes
  if (is.null(genes)) {
    if (is.null(gene_list)) {
      tax3_genes_raw <- neo4j_query(
        qry=paste("MATCH (d:Disease {diseaseId: \"", disease_ids[1], "\"})<--(s:Study)<--
                  (t:Tax3Result {threshold: \"sidak\"})-->(v:Variant)-[s2g]->(g:Gene)
                  WHERE s2g.confidence > 0 AND g.geneType = 'protein_coding' AND s2g.tissueWeights = \"", disease_names[1], "\"
                  RETURN DISTINCT g.geneName,g.ensId,v.rsId,v.chromosome,g.mhcRegion,s2g.confidence;", sep=''),
        con=con_neo4j
                  )
      tax3_genes_raw$s2g.confidence <- as.numeric(tax3_genes_raw$s2g.confidence)
      # filter for highest topsis SNP-gene pairs
      index <- c()
      for (gene in unique(tax3_genes_raw$g.ensId)) {
        if (startsWith(gene, 'ENSG')){
          gene_rows <- tax3_genes_raw[tax3_genes_raw$g.ensId == gene,]
          gene_rows <- gene_rows[order(gene_rows$s2g.confidence, decreasing=TRUE),]
          index <- c(index, as.numeric(rownames(gene_rows)[1]))
        }
      }
      tax3_genes <- tax3_genes_raw[index,]
    } else {
      # Retreive using GeneList
      tax3_genes_raw <- neo4j_query(
        qry=paste("MATCH (g:Gene)-->(gl:GeneList {listName: \"", gene_list, "\"})
                    OPTIONAL MATCH (v:Variant)-[s2g]->(g:Gene)
                    WHERE s2g.tissueWeights = \"", disease_names[1], "\"
                  RETURN DISTINCT g.geneName,g.ensId,v.rsId,v.chromosome,g.mhcRegion,s2g.confidence;", sep=''),
        con=con_neo4j
      )
      tax3_genes_raw$s2g.confidence[is.null(tax3_genes_raw$s2g.confidence)] <- 0
      tax3_genes_raw$s2g.confidence <- as.numeric(tax3_genes_raw$s2g.confidence)
      # filter for highest topsis SNP-gene pairs
      index <- c()
      for (gene in unique(tax3_genes_raw$g.ensId)) {
        if (startsWith(gene, 'ENSG')){
          gene_rows <- tax3_genes_raw[tax3_genes_raw$g.ensId == gene,]
          gene_rows <- gene_rows[order(gene_rows$s2g.confidence, decreasing=TRUE),]
          index <- c(index, as.numeric(rownames(gene_rows)[1]))
        }
      }
      tax3_genes <- tax3_genes_raw[index,]
    }
  } else {
    # Retrieve from list of ensIds
    tax3_genes_raw <- neo4j_query(
      qry=paste("MATCH (g:Gene)
                WHERE g.ensId IN [\"", paste(genes, collapse="\", \""), "\"]
                OPTIONAL MATCH (v:Variant)-[s2g]->(g:Gene)
                WHERE s2g.tissueWeights = \"", disease_names[1], "\"
                RETURN DISTINCT g.geneName,g.ensId,v.rsId,v.chromosome,g.mhcRegion,s2g.confidence;", sep=''),
      con=con_neo4j
    )
    tax3_genes_raw$s2g.confidence[is.null(tax3_genes_raw$s2g.confidence)] <- 0
    tax3_genes_raw$s2g.confidence <- as.numeric(tax3_genes_raw$s2g.confidence)
    # filter for highest topsis SNP-gene pairs
    index <- c()
    for (gene in unique(tax3_genes_raw$g.ensId)) {
      if (startsWith(gene, 'ENSG')){
        gene_rows <- tax3_genes_raw[tax3_genes_raw$g.ensId == gene,]
        gene_rows <- gene_rows[order(gene_rows$s2g.confidence, decreasing=TRUE),]
        index <- c(index, as.numeric(rownames(gene_rows)[1]))
      }
    }
    tax3_genes <- tax3_genes_raw[index,]
  }
  print(paste(nrow(tax3_genes), 'genes retrieved'))
  for (i in 1:length(disease_ids)) {
    print(paste('Collecting data  for', disease_names[i]))
    da_query <- paste("MATCH (d:Disease {diseaseId: \"", disease_ids[i],
                      "\"})-[ad:ASSOCIATED_TARGET]->(g:Gene {geneType: \"protein_coding\"})
                       RETURN DISTINCT g.ensId,ad.openTargetsScore,ad.openTargetsGeneticassociation,
                              ad.openTargetsKnownDrug,ad.openTargetsAffectedPathway,ad.openTargetsRnaExpression,
                              ad.openTargetsLiterature,ad.openTargetsAnimalModel;", sep='')
    da_df <- neo4j_query(
      qry=da_query,
      con=con_neo4j
    )
    if (is.null(da_df)) {
      da_query <- paste("MATCH (d:Disease {diseaseId: \"", disease_ids[i],
                        "\"})<-[ad:ASSOCIATED_DISEASE]-(g:Gene {geneType: \"protein_coding\"})
                       RETURN DISTINCT g.ensId,ad.openTargetsScore,ad.openTargetsGeneticassociation,
                              ad.openTargetsKnownDrug,ad.openTargetsAffectedPathway,ad.openTargetsRnaExpression,
                              ad.openTargetsLiterature,ad.openTargetsAnimalModel;", sep='')
      da_df <- neo4j_query(
        qry=da_query,
        con=con_neo4j
      )
    }
    if (is.null(da_df) == FALSE) {
      names(da_df) <- c('g.ensId', paste(
        disease_names[i], c('_OT_Score', '_OT_Genetic', '_OT_Drug', '_OT_Pathway', '_OT_RNA', '_OT_Literature', '_OT_Animal'), sep=''))
      tax3_genes <- merge(tax3_genes, da_df, all.x=TRUE)
    }
  }
  tax3_genes[is.na(tax3_genes)] <- 0
  for (i in 5:ncol(tax3_genes)) {
    tax3_genes[, i] <- as.numeric(tax3_genes[, i])
  }
  if (include_zeros == FALSE) {
    index <- c()
    for (i in 1:nrow(tax3_genes)) {
      if (sum(as.numeric(tax3_genes[i,7:ncol(tax3_genes)])) > 0) {
        index <- c(index, i)
      }
    }
    tax3_genes <- tax3_genes[index, ]
  }
  print(paste(nrow(tax3_genes), 'genes displayed in output'))
  hm_mat <- as.matrix(tax3_genes[,c(5:ncol(tax3_genes))])
  if (hide_genes) {
    rownames(hm_mat) <- as.character(seq(length(tax3_genes$g.geneName)))
  } else {
    rownames(hm_mat) <- tax3_genes$g.geneName
  }
  
  gene_clust <- hclust(dist(hm_mat[,2:ncol(hm_mat)]), method='ward.D2')
  write.table(tax3_genes[rev(gene_clust$order),],
              paste(output_prefix, '.txt', sep=''),
              sep='\t', quote=FALSE, row.names=FALSE)
  
  png(width=png_width, height=png_height, filename=paste(output_prefix, '.png', sep=''))
  par(mar=c(1,1,1,1))
  heatmap3(x=hm_mat[gene_clust$order,],
           col=colorRampPalette(c("#FBF8FE", "#2A0D43"))(n=100),
           Colv=NA, Rowv=NA, balanceColor=FALSE, scale="none",
           cexRow=0.7, cexCol=0.7)
  dev.off()
}
