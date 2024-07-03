library(dplyr)

setwd('/Users/macbook/Desktop/Bioinf/Part-c')

# Read the text files into data frames
file1=read.table("count_SRR6191645/counts.txt",sep='\t', header = FALSE)
file2=read.table("count_SRR6191646/counts.txt",sep='\t', header = FALSE)

# set column names
# based on SRR information we know:
#   GSM2808516: 48N; Homo sapiens; RNA-Seq (SRR6191646)
#   GSM2808515: 48C; Homo sapiens; RNA-Seq (SRR6191645)
colnames(file1) <- c("GeneEns", "X48C_COUNT")
colnames(file2) <- c( "GeneEns", "X48N_COUNT")

# number of gene that are not expressed in SRR6191645
number_of_genes_not_expressed_in_X48C <- sum(file1$X48C_COUNT == 0)
print(number_of_genes_not_expressed_in_X48C)
# number of gene that are not expressed in SRR6191646
number_of_genes_not_expressed_in_X48N <- sum(file2$X48N_COUNT == 0)
print(number_of_genes_not_expressed_in_X48N)

# Merge the two data frames by the first column
merged_data <- merge(file1, file2, by = "GeneEns", all = TRUE)

number_of_genes_not_expressed_in_Both <- sum(merged_data$X48N_COUNT == 0 & merged_data$X48C_COUNT == 0)
print(number_of_genes_not_expressed_in_Both)


total_sum_C <- sum(file1$X48C_COUNT, na.rm = TRUE)
print(total_sum_C)
total_sum_N <- sum(file2$X48N_COUNT, na.rm = TRUE)
print(total_sum_N)

# Print the sum
print(total_sum)

write.table(merged_data, "merged_data.txt", sep = "\t", row.names = FALSE, quote = FALSE)
