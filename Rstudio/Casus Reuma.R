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
align.eth1 = align(index = "ref_human", readfile1 = "SRR8394576_ethanol_12h_1.fasta.gz", output_file = "eth1.BAM")
align.eth2 = align(index = "ref_human", readfile1 = "SRR8394577_ethanol_12h_2.fasta.gz", output_file = "eth2.BAM")
align.eth3 = align(index = "ref_human", readfile1 = "SRR8394578_ethanol_12h_3.fasta.gz", output_file = "eth3.BAM")
align.con1 = align(index = "ref_human", readfile1 = "SRR8394612_control_12h_1.fasta.gz", output_file = "con1.BAM")
align.con2 = align(index = "ref_human", readfile1 = "SRR8394613_control_12h_2.fasta.gz", output_file = "con2.BAM")
align.con3 = align(index = "ref_human", readfile1 = "SRR8394614_control_12h_3.fasta.gz", output_file = "con3.BAM")