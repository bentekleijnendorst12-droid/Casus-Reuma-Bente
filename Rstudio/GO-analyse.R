#Go analyse
GO.wall=goseq(pwf,"hg19", "geneSymbol")

#wat is enriched
class(GO.wall)
head(GO.wall)
nrow(GO.wall)
