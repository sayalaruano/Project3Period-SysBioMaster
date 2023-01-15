require(rstudioapi)
require(biomaRt)
require(clusterProfiler)
require(Biobase)
require(GEOquery)
require(org.Hs.eg.db)

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
  dataset <- getGEO(dataset_name, destdir = ".")
  
  # Transform to the same class
  if (class(dataset) == "list") {
    dataset = dataset[[1]]
  } else if (class(dataset) == "GDS") {
    dataset <- GDS2eSet(dataset)
  }
  return(dataset)
}

#gse1869 = extract_dataset("GSE1869")
#gds4772 = extract_dataset("GDS4772")
#gse5406 = extract_dataset("GSE5406")

# Extract Entrez ID for translation between our gx data and the datasets
geneID_gxData = bitr(
  rownames(gxData),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
colnames(geneID_gxData) = c("ENSEMBL_ID", "ENTREZ_ID")
# Remove duplicates
geneID_gxData = geneID_gxData[!duplicated(geneID_gxData$ENTREZ_ID),]

#### ROSETTA ID CREATION #######################################################
create_Rosetta = function(dataset,
                          DATASETID_colname = "ID",
                          ENTREZID_colname = "ENTREZ_GENE_ID") {
  # Create Data frame with the Entrez ID of our datasets
  geneID_dataset = data.frame(DATASET_ID = as.character(fData(dataset)[, DATASETID_colname]),
                              ENTREZ_ID = as.character(fData(dataset)[, ENTREZID_colname]))
  # Remove duplicates
  geneID_dataset = geneID_dataset[!duplicated(geneID_dataset$ENTREZ_ID),]
  
  rownames(geneID_gxData) = geneID_gxData$ENTREZ_ID
  rownames(geneID_dataset) = geneID_dataset$ENTREZ_ID
  
  # Select genes in both datasets
  common_genes = intersect(rownames(geneID_gxData), rownames(geneID_dataset))
  geneID_gxData = geneID_gxData[common_genes,]
  geneID_dataset = geneID_dataset[common_genes,]
  
  
  # Create a new object with all the common genes and their three IDs
  geneID_Rosetta = merge(geneID_gxData, geneID_dataset, by = "ENTREZ_ID")
}

#rosetta4772 = create_Rosetta(gds4772, "ID", "Gene ID")
#rosetta1869 = create_Rosetta(gse1869, "ID", "ENTREZ_GENE_ID")
#rosetta5406 = create_Rosetta(gse5406, "ID", "ENTREZ_GENE_ID")

#### EXPRESSION DATA EXTRACTION ################################################
extract_expData = function(dataset, rosetta) {
  # Extract and translate to Ensembl the genes the dataset
  expression_data = assayData(dataset)$exprs[rosetta$DATASET_ID,]
  rownames(expression_data) = rosetta$ENSEMBL_ID
  
  # Add the genes not in the dataset
  missing_genes = setdiff(rownames(gxData), rownames(expression_data))
  missing_expression = as.data.frame(matrix(
    nrow = length(missing_genes),
    ncol = ncol(expression_data),
    dimnames = list(missing_genes, colnames(expression_data))
  ))
  
  # Putting it together and sorting according to MAGNet gene order
  expression_data = rbind(expression_data, missing_expression)
  expression_data = expression_data[rownames(gxData),]
}

#expData4772 = extract_expData(gds4772, rosetta4772)
#expData1869 = extract_expData(gse1869, rosetta1869)
#expData5406 = extract_expData(gse5406, rosetta5406)

#### PUTTING IT ALL TOGETHER ###################################################
translate_dataset = function(dataset_name) {
  dataset = extract_dataset(dataset_name)
  
  rosetta = create_Rosetta(dataset)
  
  expData = extract_expData(dataset, rosetta)
}