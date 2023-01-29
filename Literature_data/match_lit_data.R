# Imports
library(readr)

# Load data
lit_genes = read_csv("Literature_data/DCM_lit_genes.csv")
most_imp_genes = read_csv("Results_important_genes/DCM/Consensus/Consensus3methods_top300important_genes_DCM.csv")

# Select the ENSEMBL ids
most_imp_genes = most_imp_genes['Feature']

# Get the consensus
df = merge(x = lit_genes, y = most_imp_genes, by = "Feature")
