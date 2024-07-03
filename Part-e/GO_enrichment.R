cat('\14')

library(goseq)
library(org.Hs.eg.db)
library(GO.db)
library(ggplot2)
library(dplyr)
library(calibrate)
library(AnnotationDbi)
library(pathview)
library(clusterProfiler)


setwd('/Users/macbook/Desktop/Bioinf/Part-e')

DEG = read.csv("DEG_FULL_et.csv");


#-------------------- Volcano-Plot --------------------#

jpeg("Volcano_Plot.jpg", width = 9, height = 5, units = 'in', res = 500)
# Create the plot framework
with(DEG, plot(logFC, -log10(FDR), pch=20, main="Volcano plot", xlim=c(-5,7), cex=.8))
# Green points for genes with FDR < 0.1
with(subset(DEG, FDR < 0.1), points(logFC, -log10(FDR), pch=20, col="green", cex=0.8))
# Red points for genes with FDR < 0.1 and logFC > 1.5
with(subset(DEG, FDR < 0.1 & logFC > 1.5), points(logFC, -log10(FDR), pch=20, col="red", cex=0.8))
# Blue points for genes with FDR < 0.1 and logFC < -1.5
with(subset(DEG, FDR < 0.1 & logFC < -1.5), points(logFC, -log10(FDR), pch=20, col="blue", cex=0.8))
# Adding a legend to the plot
legend("topright", legend=c("logFC > 1.5", "logFC < -1.5", "FDR < 0.1 not DE"),
       col=c("red", "blue", "green"), pch=20, cex=0.8, text.width = 2, text.font = 2, x.intersp=0.8, y.intersp=0.8)

dev.off()

# Filter the DEG dataset
DEG_filtered <- subset(DEG, FDR < 0.1 & abs(logFC) > 1.5)

write.csv(DEG_filtered, "DEG_Filtered.csv", row.names = FALSE)
write.table(DEG_filtered, "DEG_Filtered.txt", row.names = FALSE, sep = "\t", quote = FALSE)



#-------------------- ENTREZID Conversion --------------------#

# Map gene symbols to Entrez IDs
entrez_ids <- unlist(mget(DEG_filtered$X, org.Hs.egSYMBOL2EG, ifnotfound = NA))

# Get the Entrez ID for CXCL8
cxcl8_entrez_id <- unlist(mget("CXCL8", org.Hs.egSYMBOL2EG, ifnotfound = NA))
# Replace the IL8 value with CXCL8's Entrez ID
entrez_ids["IL8"] <- cxcl8_entrez_id


#-------------------- GOseq --------------------#

# Getting the length for pwd function
gene_info <- read.table("GSE104836_gene_exp.txt", header = TRUE, sep = "\t")
match_indices <- match(DEG_filtered$X, gene_info$gene)
gene_lengths <- ifelse(is.na(gene_info$Length[match_indices]), 1000, gene_info$Length[match_indices])
names(gene_lengths) <- entrez_ids

# Array indicating whether each gene is significant of not
de_genes <- rep(1, length(entrez_ids))
names(de_genes) <- entrez_ids

deg_list <- names(de_genes)
deg_list

#-------------------- KEGG --------------------#

# perform the enrichment
kegg_enriched <- enrichKEGG(gene = deg_list, organism = 'hsa')
# Extract the results
kegg_results <- as.data.frame(kegg_enriched)
# Order results by the number of genes in descending order for better visualization
kegg_results_ordered <- kegg_results[order(-kegg_results$Count),]

# First, we assign the plot to a variable
kegg_plot <- ggplot(kegg_results_ordered, aes(x = reorder(Description, Count), y = Count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  coord_flip() + # Flip coordinates for horizontal bars
  xlab("KEGG Pathways") +
  ylab("Count of DEGs") +
  ggtitle("DEGs Count Mapped to KEGG Pathways") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot to a JPEG file
ggsave("kegg_plot.jpg", plot = kegg_plot, width = 12, height = 9, units = "in", dpi = 300)


#-------------------- GO enrichment --------------------#

# Function to perform GO enrichment and plot the results
plot_go_enrichment <- function(degs, OrgDb, ont, title_suffix) {
  go_enriched <- enrichGO(gene = degs,
                          OrgDb = OrgDb,
                          keyType = "ENTREZID",
                          ont = ont,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)
  
  go_results <- as.data.frame(go_enriched)
  go_results_ordered <- go_results[order(-go_results$Count),]
  
  ggplot(go_results_ordered, aes(x=reorder(Description, Count), y=Count)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    coord_flip() +
    xlab(paste0(ont, " Term")) +
    ylab("Count of DEGs") +
    ggtitle(paste0("DEGs Count Mapped to GO:", title_suffix)) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Biological Process (BP)
bp_plot <- plot_go_enrichment(deg_list, org.Hs.eg.db, "BP", " Biological Process")
ggsave("GO_BP_plot.jpg", plot = bp_plot, width = 12, height = 9, units = "in", dpi = 300)

# Molecular Function (MF)
mf_plot <- plot_go_enrichment(deg_list, org.Hs.eg.db, "MF", " Molecular Function")
ggsave("GO_MF_plot.jpg", plot = mf_plot, width = 12, height = 9, units = "in", dpi = 300)

# Cellular Component (CC)
cc_plot <- plot_go_enrichment(deg_list, org.Hs.eg.db, "CC", " Cellular Component")
ggsave("GO_CC_plot.jpg", plot = cc_plot, width = 12, height = 9, units = "in", dpi = 300)




# Prepare the pwf object with gene lengths
pwf <- nullp(de_genes, bias.data=gene_lengths, genome='hg19', id='knownGene')
GO.wall = goseq(pwf,"hg19","knownGene")
GO.wall <- GO.wall[GO.wall$numInCat >= 3, ] #only show annotations where more than 2 DEGs has been mapped to it

# Subset the dataframe for each of the categories (BP, MF, CC)
CC_results <- subset(GO.wall, ontology == "CC")

# Define a function to create plots
create_goseq_plot <- function(data, title) {
 ggplot(data, aes(x=reorder(term, -numDEInCat), y=numDEInCat)) +
   geom_bar(stat="identity") +
   coord_flip() + # Flips the coordinates to make the terms readable on the y-axis
   labs(title = title, x = "Number of differentially expressed genes", y = "GO/KEGG Term") +
   theme(plot.title = element_text(hjust = 0.5))
}

# Create plots
cc_plot <- create_goseq_plot(CC_results, "Cellular Component")

# Optionally, save the plots as image files
ggsave("cc_plot.jpg", cc_plot, width = 10, height = 20)

