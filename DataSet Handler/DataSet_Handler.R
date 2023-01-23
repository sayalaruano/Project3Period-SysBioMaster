require(rstudioapi)
require(biomaRt)
require(clusterProfiler)
require(Biobase)
require(GEOquery)
require(org.Hs.eg.db)
require(caroline)

setwd(paste0(dirname(
  rstudioapi::getSourceEditorContext()$path
), "/Data"))
gxData <-
  read.delim("MAGNET_GeneExpressionData_CPM_19112020.txt", row.names = 1)

#### DATASET EXTRACTION #######################################################
# Extract databases and transform them to a common object type (done manually)
extract_dataset = function(dataset_name) {
  setwd(paste0(
    dirname(rstudioapi::getSourceEditorContext()$path),
    "/DataSets"
  ))
  # Extract dataset
  dataset <- read.delim(dataset_name)
  return(dataset)
}

#### ROSETTA ID CREATION #######################################################
create_Rosetta = function(dataset_ID_type) {
  # Extract Entrez ID for translation between our gx data and the datasets
  if(dataset_ID_type != "ENSEMBL"){
  geneID_gxData = bitr(
    rownames(gxData),
    fromType = "ENSEMBL",
    toType = dataset_ID_type,
    OrgDb = org.Hs.eg.db
  )
  
  # Remove duplicates
  geneID_gxData = geneID_gxData[!duplicated(geneID_gxData[, 1]), ]
  geneID_gxData = geneID_gxData[!duplicated(geneID_gxData[, 2]), ]
  rownames(geneID_gxData) = geneID_gxData[,1]
  # Standardize colnames
  colnames(geneID_gxData) = c("ENSEMBL_ID", paste0(dataset_ID_type, "_ID"))
  }else{
    geneID_gxData = data.frame(ENSEMBL_ID = rownames(gxData),
                               ENSEMBL2_ID = rownames(gxData),
                               row.names = rownames(gxData))
  }
  return(geneID_gxData)
}

#### FIT TO ONLY gxData GENES ##################################################
fit_to_gxData = function(dataset, rosetta, genes.in.rows = TRUE) {
  # Get the list of genes
  gxData_genelist = rosetta[, 2]
  if (genes.in.rows) {
    expData_genelist = rownames(dataset)
  } else{
    expData_genelist = dataset[, 1]
  }
  
  # Get the genes that appear in gxData and the dataset
  common_genes = intersect(gxData_genelist, expData_genelist)
  # Get the genes that appear in gxData but not in the dataset
  missing_genes = setdiff(gxData_genelist, expData_genelist)
  # Create an empty dataframe with the missing genes
  missing_expression = as.data.frame(matrix(
    nrow = length(missing_genes),
    ncol = ncol(dataset),
    dimnames = list(rosetta[missing_genes, 1], colnames(dataset))
  ))
  
  # Get the genes that have no translation
  genes_no_translation = setdiff(rownames(gxData), rosetta[,1])
  # Create an empty dataframe with the untranslated genes
  no_translation_expression = data.frame(matrix(
    nrow = length(genes_no_translation),
    ncol = ncol(dataset),
    dimnames = list(genes_no_translation, colnames(dataset))
  ))
  # Combine the genes in common with the genes missing
  if (genes.in.rows) {
    fit_dataset = rbind(dataset[common_genes, ], missing_expression)
  } else{
    missing_expression[,1] = rownames(missing_expression)
    fit_dataset = rbind(dataset , missing_expression)
    rownames(fit_dataset) = fit_dataset[, 1]
    fit_dataset = fit_dataset[,-1]
  }
  
  # Order them according to the genes in gxData
  fit_dataset = fit_dataset[rosetta[,2],]
  rownames(fit_dataset) = rosetta[,1]
  fit_dataset = rbind(fit_dataset, no_translation_expression)
  return(fit_dataset)
}
#### EXPORT MANIPULATED DATA ###################################################
export_Data = function(dataset){
  setwd(paste0(
    dirname(rstudioapi::getSourceEditorContext()$path),
    "/Translated_DataSets"
  ))
  
  write.delim(dataset, paste0(deparse(substitute(dataset)), ".txt"), row.names = TRUE)
}


#### FEATURE NORMALIZATION #####################################################
rescale_features = function(dataset){
  max_gxData = rowMax(as.matrix(gxData))
  min_gxData = rowMin(as.matrix(gxData))
  
  dataset[is.na(dataset)] = -Inf
  max_dataset = rowMax(as.matrix(dataset))
  dataset[dataset == -Inf] = Inf
  min_dataset = rowMin(as.matrix(dataset))
  dataset[dataset == Inf] = NA
  
  rescaled_dataset = ((dataset - min_dataset)/(max_dataset - min_dataset)*(max_gxData - min_gxData))+min_gxData
  rescaled_dataset[is.nan(as.matrix(rescaled_dataset))] = NA
  
  return(rescaled_dataset)
}
