
## Analysis of Differentially Expressed LncRNAs and mRNA in Colorectal Cancer

This project is based on the study "Differentially Expressed LncRNAs and mRNA Identified by Sequencing Analysis in Colorectal Cancer Patients" (GEO accession number: GSE104836). It involves 20 samples, with 10 normal and 10 cancerous tissues. Paired samples with IDs 48 (GEO accessions GSM2808515 and GSM2808516; SRA run accessions SRR6191645 and SRR6191646) were specifically selected.

### Workflow Summary:

1. **Quality Control and Trimming**:
   - SRA entries were converted to paired Fastq.gz files using `fastq-dump`.
   - Quality control was performed using FastQC to ensure high data quality.

2. **Read Mapping**:
   - Reads were mapped to the reference genome using HISAT2 for alignment.

3. **Building Gene Expression Matrix**:
   - Aligned reads were processed with `htseq-count` to build a gene expression matrix for differential expression analysis.

4. **Differential Gene Expression Analysis**:
   - Conducted to identify differentially expressed genes between normal and cancerous tissues.
