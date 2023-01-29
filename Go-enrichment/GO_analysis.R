# Imports
library(biomaRt)
library(readr)
library(DOSE)

# Load data
most_imp_genes = read_csv("Results_important_genes/DCM/Consensus/Consensus3methods_top700important_genes_DCM.csv")
list_genes = read_csv("Results_important_genes/DCM/Important_genes_DCM_RandomForest_sorted.csv")

# Select the ENSEMBL ids
most_imp_genes = most_imp_genes['Feature']
list_genes = as.data.frame(list_genes['Feature'])

# Choose the Ensembl database and set the host
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# Search the specific attibute names of gene names, gene symbols, and chromosome names
entrezgeneid_att = searchAttributes(mart = ensembl, pattern = "entrezgene")

# Retrieve gene symbols, gene names, and chromosome names for the specified Ensembl gene identifiers
ensembl_out_impgenes = getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                         filters = "ensembl_gene_id",
                         values = most_imp_genes,
                         mart = ensembl)

write_csv(ensembl_out_impgenes['hgnc_symbol'], 'hgnc_symbols_impgenes_700_DCM.csv')

impgenes_entrezids = na.omit(ensembl_out_impgenes['entrezgene_id'])
impgenes_entrezids = lapply(impgenes_entrezids, as.character)
impgenes_entrezids = impgenes_entrezids$entrezgene_id

ensembl_out_allgenes = getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                              filters = "ensembl_gene_id",
                              values = list_genes,
                              mart = ensembl)

allgenes_entrezids = na.omit(ensembl_out_allgenes['entrezgene_id'])
allgenes_entrezids = lapply(allgenes_entrezids, as.character)
allgenes_entrezids = allgenes_entrezids$entrezgene_id

# Disease over-representation analysis
x <- enrichDO(gene          = impgenes_entrezids,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = allgenes_entrezids,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
x

dgn <- enrichDGN(impgenes_entrezids) 

dgn

y <- gseDO(impgenes_entrezids,
           minGSSize     = 120,
           pvalueCutoff  = 0.2,
           pAdjustMethod = "BH",
           verbose       = FALSE)
