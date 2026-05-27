#working directory intellen
setwd("C:/Users/bente/Desktop/school/Casus Reuma/")
getwd()

#instaleren package
install.packages('BiocManager')
BiocManager::install('Rsubread')

#package inladen
library('Rsubread')

#indexeren
buildindex(
  basename = 'ref_human',
  reference = 'GCF_000001405.40_GRCh38.p14_genomic.fna',
  memory = 4000,
  indexSplit = TRUE)

#mappen
align.norm1 = align(index = "ref_human", readfile1 = "SRR4785819_1_subset40k.fastq", readfile2 = "SRR4785819_2_subset40k.fastq", output_file = "norm1.BAM")
align.norm2 = align(index = "ref_human", readfile1 = "SRR4785820_1_subset40k.fastq", readfile2 = "SRR4785820_2_subset40k.fastq", output_file = "norm2.BAM")
align.norm3 = align(index = "ref_human", readfile1 = "SRR4785828_1_subset40k.fastq", readfile2 = "SRR4785828_2_subset40k.fastq", output_file = "norm3.BAM")
align.norm4 = align(index = "ref_human", readfile1 = "SRR4785831_1_subset40k.fastq", readfile2 = "SRR4785831_2_subset40k.fastq", output_file = "norm4.BAM")
align.reuma1 = align(index = "ref_human", readfile1 = "SRR4785979_1_subset40k.fastq", readfile2 = "SRR4785979_2_subset40k.fastq", output_file = "reuma1.BAM")
align.reuma2 = align(index = "ref_human", readfile1 = "SRR4785980_1_subset40k.fastq", readfile2 = "SRR4785980_2_subset40k.fastq", output_file = "reuma2.BAM")
align.reuma3 = align(index = "ref_human", readfile1 = "SRR4785986_1_subset40k.fastq", readfile2 = "SRR4785986_2_subset40k.fastq", output_file = "reuma3.BAM")
align.reuma4 = align(index = "ref_human", readfile1 = "SRR4785988_1_subset40k.fastq", readfile2 = "SRR4785988_2_subset40k.fastq", output_file = "reuma4.BAM")

#vector maken
alles= c("norm1.BAM","norm2.BAM","norm3.BAM", "norm4.BAM", "reuma1.BAM", "reuma2.BAM","reuma3.BAM", "reuma4.BAM")

# count matrix maken
count_matrix <- featureCounts(
  files = alles,
  annot.ext = "genomic.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE, 
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)

#counts van de countmatrix opslaan onder een nieuwe naam
counts <- count_matrix$counts

# de counts van de countmatrix nieuwe kolomnamen geven
colnames(counts) <- c("norm1", "norm2", "norm3", "norm4", "reuma1", "reuma2", "reuma3", "reuma4")
head(counts)

#opslaan
write.csv(counts, "subsetcasus1_countmatrix.csv")

#gekregen countmatrix inlezen
echtedata=read.table("count_matrix_RA.txt", header = TRUE, row.names = 1)

#packages inladen
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("KEGGREST")
library(KEGGREST)
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
BiocManager::install("pathview")
library(pathview)

# colnames geven
colnames(echtedata) <- c("norm1", "norm2", "norm3", "norm4", "reuma1", "reuma2", "reuma3", "reuma4")

#rijnamen geven
rownames(Status_table)=c("norm1", "norm2", "norm3", "norm4", "reuma1", "reuma2", "reuma3", "reuma4")

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

#tabel opslaan
write.table(resultatenR, file = 'ResultatenReuma.csv', row.names = TRUE, col.names = TRUE)

# sorteren
sum(resultatenR$padj < 0.05 & resultatenR$log2FoldChange > 1, na.rm = TRUE)
sum(resultatenR$padj < 0.05 & resultatenR$log2FoldChange < -1, na.rm = TRUE)

#welke genen zijn intresant
hoogste_fold_changeR <- resultatenR[order(resultatenR$log2FoldChange, decreasing = TRUE), ]
laagste_fold_changeR <- resultatenR[order(resultatenR$log2FoldChange, decreasing = FALSE), ]
laagste_p_waardeR <- resultatenR[order(resultatenR$padj, decreasing = FALSE), ]

hoogste_fold_changeR
laagste_fold_changeR
laagste_p_waardeR

#volcano plot maken
EnhancedVolcano(resultatenR,
                lab = rownames(resultatenR),
                x = 'log2FoldChange',
                y = 'padj')

#figuur opslaan
dev.copy(png, 'VolcanoplotCasus.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

#packages inladen voor GO analyse
BiocManager::install("goseq")
library(goseq)
BiocManager::install("geneLenDataBase")
library(geneLenDataBase)
BiocManager::install("org.Dm.eg.db")
library(org.Dm.eg.db)
library("magrittr")
library(dplyr)

# tabel maken met de resultaten van de DESeq
All= rownames(resultatenR)
ALLtabel= as.data.frame(resultatenR)
ALLtabel

#Deg gefilterd op een log foldchange kleiner en gelijk aan -1 en groter of gelijk aan 1 en een adjusted p kleiner dan 0.05
DEG= ALLtabel %>%
  filter(padj<0.05)%>%
  filter(log2FoldChange<=-1| log2FoldChange >=1)
DEG
DEGG=rownames(DEG)
class(DEGG)

#kijken welke van All in DEGG staan of andersom
gene.vector=as.integer(All%in%DEGG)
gene.vector

#packages inladen
BiocManager::install("dplyr")
library(dplyr)
library(stringr)
library(goseq)

#human genome optis bekijken, hg19 is gekozen
supportedOrganisms() %>% filter(str_detect(Genome, "hg"))
names(gene.vector)=All
pwf=nullp(gene.vector,"hg19","geneSymbol")

#Go analyse
GO.wall=goseq(pwf,"hg19", "geneSymbol")

#wat is enriched
class(GO.wall)
head(GO.wall)
nrow(GO.wall)

#kijken welke enriched statistisch significant is
enriched.GO=GO.wall$category[GO.wall$over_represented_pvalue<.05]

#kijken hoeveel GO termen er overblijven
class(enriched.GO)
head(enriched.GO)
length(enriched.GO)

# packages inladen voor het maken van een dotplot
library(ggplot2)
library(dplyr)
BiocManager::install("pathview")
library(pathview)

#Selecteer de top 10 meest significante pathways en bereken -log10(p-waarde)
top10_GO <- GO.wall %>%
  arrange(over_represented_pvalue) %>%
  head(10) %>%
  mutate(log_p = -log10(over_represented_pvalue),
         # Zorgt ervoor dat de plot gesorteerd blijft op significantie
         term = reorder(term, log_p))

#Maak de dotplot
ggplot(top10_GO, aes(x = log_p, y = category)) +
  geom_point(aes(size = numDEInCat, color = log_p)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Top 10 Verrijkte GO Termen",
    x = expression(-log[10] * "(p-value)"),
    y = "GO Term",
    size = "Aantal DE genen",
    color = "-log10(p)"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))

#bar plot is beter
#Maak de bar chart
ggplot(top10_GO, aes(x = numInCat/numDEInCat, y =  reorder(category, -numInCat/numDEInCat), fill = ontology)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Top 10 Meest Significante GO Termen",
    x = "Verhouding genexpressie Reuma en Gezond",
    y = "GO Term",
    fill = "Ontologie"
  ) +
  theme_minimal()

#pathway visualisatie
ALLtabel[1]=NULL
ALLtabel[2:5]=NULL

pathview(
  gene.data = ALLtabel,
  pathway.id = "04672",  
  species = "hsa",          
  gene.idtype = "SYMBOL",     
  limit = list(gene = 5))

