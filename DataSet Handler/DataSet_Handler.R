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
  
  # Transform to the same class
  if (class(dataset) == "list") {
    dataset = dataset[[1]]
  } else if (class(dataset) == "GDS") {
    dataset <- GDS2eSet(dataset)
  }
  return(dataset)
}

#### ROSETTA ID CREATION #######################################################
create_Rosetta = function(dataset,
                          dataset_ID_type,
                          genes.in.rows = FALSE) {
  # Extract Entrez ID for translation between our gx data and the datasets
  geneID_gxData = bitr(
    rownames(gxData),
    fromType = "ENSEMBL",
    toType = dataset_ID_type,
    OrgDb = org.Hs.eg.db
  )
  
  # Remove duplicates
  geneID_gxData = geneID_gxData[!duplicated(geneID_gxData[, 1]), ]
  geneID_gxData = geneID_gxData[!duplicated(geneID_gxData[, 2]), ]
  # Standardize colnames
  colnames(geneID_gxData) = c("ENSEMBL_ID", paste0(dataset_ID_type, "_ID"))
  rownames(geneID_gxData) = geneID_gxData[, 2]
  
  # Create data frame with the ID of our datasets
  if (genes.in.rows) {
    geneID_dataset = rownames(dataset)
  } else{
    geneID_dataset = dataset[, 1]
  }
  # Remove duplicates
  geneID_dataset = data.frame(geneID_dataset[!duplicated(geneID_dataset)])
  rownames(geneID_dataset) = geneID_dataset[, 1]
  colnames(geneID_dataset) = paste0(dataset_ID_type, "_ID")
  
  # Select genes in both datasets
  common_genes = intersect(geneID_dataset[, 1], geneID_gxData[, 2])
  
  # Create a new object with all the common genes and their three IDs
  geneID_Rosetta = geneID_gxData[common_genes, ]
}

#### FIT TO ONLY gxData GENES ##################################################
fit_to_gxData = function(dataset, rosetta, genes.in.rows = FALSE) {
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
    dimnames = list(missing_genes, colnames(dataset))
  ))
  
  # Combine the genes in common with the genes missing
  if (genes.in.rows) {
    fit_dataset = rbind(dataset[common_genes, ], missing_expression)
    gene_list = data.frame(Gene = rownames(fit_dataset), row.names = rownames(fit_dataset))
    fit_dataset = cbind(gene_list, fit_dataset)
  } else{
    rownames(dataset) = dataset[, 1]
    fit_dataset = rbind(dataset[common_genes, ], missing_expression)
  }
  
  # Order them according to the genes in gxData
  fit_dataset = fit_dataset[rownames(gxData), ]
}

#### PUTTING IT ALL TOGETHER ###################################################
translate_dataset = function(dataset_name,
                             dataset_ID_type,
                             genes.in.rows = FALSE) {
  dataset = extract_dataset(dataset_name)
  
  rosetta = create_Rosetta(dataset, dataset_ID_type, genes.in.rows)
  
  fit_dataset = fit_to_gxData(dataset, rosetta, genes.in.rows)
  
  setwd(paste0(
    dirname(rstudioapi::getSourceEditorContext()$path),
    "/Translated_DataSets"
  ))
  write.delim(fit_dataset,
              paste0("expData_",".txt"),
              row.names = F)
  
  return(fit_dataset)
}
