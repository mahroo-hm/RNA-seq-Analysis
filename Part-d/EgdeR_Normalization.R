

################## RNASeq Normalization (edgeR) ################## 


cat('\14')
library("readxl")
library("edgeR")


setwd('/Users/macbook/Desktop/Bioinf/Part-d')
Dat=read.table('GSE104836_gene_exp.txt',sep='\t', header = TRUE)

Dat1=Dat[!duplicated(Dat$gene), ] # remove repeated gene symbols

# set the first column as gene symbols
col1<-Dat1$gene
Dat2<-Dat1[,-1]
Dat2 = as.data.frame(Dat2)
row.names(Dat2) <- col1

# For the raw count data
count_columns <- grep("_COUNT$", names(Dat2), value = TRUE)
DatRaw <- Dat2[, c(count_columns)]

# For the normalized TPM data
tpm_columns <- grep("_TPM$", names(Dat2), value = TRUE)
DatNorm <- Dat2[, c(tpm_columns)]

Expression_Raw = as.matrix(DatRaw)
Expression_Norm = as.matrix(DatNorm)

means <- rowMeans(Expression_Raw)
filter <- means >= 1
table(filter)
keepSamples = (filter==TRUE)
geneCountHigh_M_Raw <- Expression_Raw[keepSamples,]
dim(geneCountHigh_M_Raw)

means <- rowMeans(Expression_Norm)
filter <- means >= 1
table(filter)
keepSamples = (filter==TRUE)
geneCountHigh_M_Norm <- Expression_Norm[keepSamples,]
dim(geneCountHigh_M_Norm)


#-------------------- TRAIT --------------------#

header_names <- names(DatRaw)
# Use regular expressions to create a table with the header and the associated tag
labels <- sub(".*([CN])_.*", "\\1", header_names)
Group_Raw <- data.frame(Header = header_names, Label = labels)

header_names <- names(DatNorm)
labels <- sub(".*([CN])_.*", "\\1", header_names)
Group_Norm <- data.frame(Header = header_names, Label = labels)


#-------------------- TMM-Normalization --------------------#
dgellist <- DGEList(counts=geneCountHigh_M_Norm, group=factor(Group_Norm$Label))
dgellist <- calcNormFactors(dgellist,method = "TMM") #method ="upperquartile"
dgellist <- estimateCommonDisp(dgellist)
dgellist <-estimateTagwiseDisp(dgellist,trend="movingave")
Normalexpr <-dgellist$pseudo.counts

lNormexpr <- log(Normalexpr + 1)


#-------------------- QC & Plots --------------------#

pdf(file='Heatmap.pdf',width = 10,height=10)
heatmap(cor(lNormexpr))
dev.off()

# Boxplot Before Normalization
pdf(file='B_boxplot.pdf',width = 10,height=10)
par(mar=c(14,5,1,1))
boxplot(Expression_Raw,las = 2,ylim = c(0, 20),labels=FALSE,col="green")
mtext("Normal samples \n Before normalization",side=2,line = 2)
dev.off()

# Boxplot After Normalization
pdf(file='A_boxplot.pdf',width = 30,height=10)
par(mar=c(14,5,1,1))
boxplot(lNormexpr,las = 2,col="green",ylim = c(0, 20))
mtext("After normalization",side=2,line = 2)
dev.off()


pdf(file='MDS.pdf',width = 30,height=10)
plotMDS(dgellist)
dev.off()
 
#-------------------- Save 2 File --------------------#

write.table(lNormexpr, "Full_Normalized.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(lNormexpr, "Full_Normalized.txt", sep="\t",row.names=T, col.names=T, quote=F)  



#-------------------- DEG Analysis (exact test) --------------------#

dge <- DGEList(counts=lNormexpr,group=factor(Group_Norm$Label))
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

# Perform the exact test
et <- exactTest(dge, pair=c(1, 2))

# Get the top DEGs across all genes
DEG <- topTags(et, n=Inf)$table

significant_DEGs <- subset(DEG, FDR <= 0.05)
n_significant <- nrow(significant_DEGs)
print(n_significant)

# Filter DEGs for those with |log2FC| > 1.5
DEG_1_5 <- DEG[abs(DEG$logFC) > 1.5, ]

# Calculate the percentage
percentage_DEG_1_5 <- nrow(DEG_1_5) / nrow(DEG) * 100
print(percentage_DEG_1_5)

# Save the full DEG table and the filtered DEG table to files
write.table(DEG, "DEG_FULL_et.csv", row.names=TRUE, col.names=NA, sep=",")
write.table(DEG, "DEG_FULL_et.txt", row.names=TRUE, col.names=NA, sep="\t", quote=FALSE)



