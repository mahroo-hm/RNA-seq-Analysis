

cat('\14')
library("readxl")
#library("edgeR")
library("rtracklayer")
library("plyr")

##########################
setwd('/Users/macbook/Desktop/LBB/LBB Code/R codes/4-NGS/Ensemble')

Dat=read.table('Read_Count_.txt',sep='\t', header = TRUE)

col1=Dat[,1]
Dat1=Dat[,-1]
row.names(Dat1)=col1


########################

gencode_file = 'gencode.v41.annotation.gtf.gz'

#ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human
#Please download the latest version of annotation file else you can use below code
# gencode_link = paste(
#   'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41',
#   gencode_file,
#   sep = '/'
# )
# download.file(gencode_link, gencode_file, method = 'libcurl') 


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



###################################################
#Extract Expression and add Gene Type and Gene_Symbol  
row1=row.names(Dat1)
indx=which(gene_lengths@listData$ensembl_gene_id %in% row1 )
ens=gene_lengths@listData$ensembl_gene_id[indx]

indx2=which(row1 %in% ens )

Gene_Symbol  <- gene_lengths@listData$gene_name[indx]
Gene_Symbol <-as.data.frame(Gene_Symbol)
Gene_Type  <- gene_lengths@listData$gene_type[indx]
Gene_Type <-as.data.frame(Gene_Type)

full=cbind(Gene_Type,Gene_Symbol,Dat1[indx2,])

write.table(full, "Final_Expression.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(full, "Final_Expression.txt", sep="\t",row.names=T, col.names=T, quote=F) 



