# metadata
Status=c ("Gezond", "Gezond", "Gezond", "Gezond", "Reuma", "Reuma", "Reuma", "Reuma")
Status_table= data.frame(Status)

#desseq dataset
ddsR = DESeqDataSetFromMatrix(countData = echtedata,
                              colData = Status_table,
                              design = ~ Status)
ddsR = DESeq(ddsR)
resultatenR= results(ddsR)
resultatenR

# sorteren
sum(resultatenR$padj < 0.05 & resultatenR$log2FoldChange > 1, na.rm = TRUE)
sum(resultatenR$padj < 0.05 & resultatenR$log2FoldChange < -1, na.rm = TRUE)