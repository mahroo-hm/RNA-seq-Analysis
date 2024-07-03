

cat('\14')
library("readxl")
#library("edgeR")
library("rtracklayer")
library("plyr")


setwd('/Users/macbook/Desktop/Bioinf/Part-c')

Dat1 = read.table('GSE104836_gene_exp.txt',sep='\t', header = TRUE)
Dat2 = read.table("merged_data.txt",sep='\t', header = TRUE)

main_study_matrix=Dat1[!duplicated(Dat1$gene), ] # remove repeated gene symbols

col1=Dat2[,1]
subset_matrix=Dat2[,-1]
row.names(subset_matrix)=col1

########################

gencode_file = 'gencode.v41.annotation.gtf.gz'

gtf = import.gff(gencode_file, format = 'gtf', genome = 'GRCh38.p13', feature.type = 'exon')

grl = reduce(split(gtf, elementMetadata(gtf)$gene_id))
gene_lengths = ldply(grl, function(x) {
  #sum up the length of individual exons
  return(c('gene_length' = sum(width(x))))
}, .id = 'ensembl_gene_id')


genetype = unique(elementMetadata(gtf)[, c('gene_id', 'gene_type','gene_name')])
colnames(genetype)[1] = 'ensembl_gene_id'
gene_lengths = merge(genetype, gene_lengths)

gene_lengths$ensembl_gene_id = gsub('\\.[0-9]*', '', gene_lengths$ensembl_gene_id)
gene_lengths


# Assuming Dat1 has gene symbols as row names
gene_symbols = row.names(Dat1)

# Find indices where gene symbols from your data match the gene symbols in the genetype data frame
indx = which(genetype$gene_name %in% gene_symbols)

# Extract the corresponding Ensembl gene IDs
ensembl_gene_ids = genetype$ensembl_gene_id[indx]
ensembl_gene_ids <- as.data.frame(ensembl_gene_ids)
row.names(ensembl_gene_ids) <- genetype$gene_name[indx] # Set gene symbols as row names

# Get the corresponding Gene Types based on matched indices
Gene_Type = genetype$gene_type[indx]
Gene_Type <- as.data.frame(Gene_Type)
row.names(Gene_Type) <- genetype$gene_name[indx] # Set gene symbols as row names

# Combine your data with the Ensembl IDs and Gene Types
Dat1_with_ensembl <- cbind(ensembl_gene_ids, Gene_Type, Dat1)

# Print the data frame to see the structure
print(head(Dat1_with_ensembl))

###################################################
#Extract Expression and add Gene Type and Gene_Symbol  
row1=row.names(subset_matrix)
indx=which(gene_lengths@listData$ensembl_gene_id %in% row1 )
ens=gene_lengths@listData$ensembl_gene_id[indx]

indx2=which(row1 %in% ens )

Gene_Symbol  <- gene_lengths@listData$gene_name[indx]
Gene_Symbol <-as.data.frame(Gene_Symbol)
Gene_Type  <- gene_lengths@listData$gene_type[indx]
Gene_Type <-as.data.frame(Gene_Type)

full=cbind(Gene_Type,Gene_Symbol,subset_matrix[indx2,])
subset_matrix <- full

subset_matrix=subset_matrix[!duplicated(subset_matrix$Gene_Symbol), ] # remove repeated gene symbols


# Make sure both matrices have ensemble gene ids as row names
rownames(main_study_matrix) <- main_study_matrix$gene
rownames(subset_matrix) <- subset_matrix$Gene_Symbol

# Subset the main study matrix to get the corresponding columns
# Replace 'X101C_COUNT' and 'X101N_COUNT' with the actual columns names you want to compare from the main study
comparison_matrix <- main_study_matrix[, c('gene', 'X48C_COUNT', 'X48N_COUNT')]

# Ensure both submatrices only contain rows that are in both data sets
common_genes <- intersect(rownames(subset_matrix), rownames(comparison_matrix))
write.table(common_genes, "common_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

subset_matrix <- subset_matrix[common_genes,]
comparison_matrix <- comparison_matrix[common_genes,]

# Compare raw counts (assuming you are comparing raw counts)
# This could be a simple difference, a ratio, etc. Here is a basic example of obtaining a ratio
comparison_ratio <- subset_matrix[,'X48C_COUNT'] / comparison_matrix[,'X48C_COUNT']
# Analyze the differences
summary(comparison_ratio)

comparison_ratio <- subset_matrix[,'X48N_COUNT'] / comparison_matrix[,'X48N_COUNT']
summary(comparison_ratio)

