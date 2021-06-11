rm(list=ls())

#source("https://bioconductor.org/biocLite.R")
#biocLite("tximport")

library(tximport)
library(GenomicFeatures)
#import gtf
txdb <-makeTxDbFromGFF("L:/Lab-Sack/komudi/SequencingData/SequencingData_PBMCs/genome_transcript/gencode.v28.annotation.gtf")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k,  columns = "TXNAME", keytype = "GENEID")
tx2genepre <- df[,2:1]

tx2genepre$tr_id=gsub("\\..*","",tx2genepre$TXNAME)

tx2gene<-tx2genepre[,c(3,2)]
names(tx2gene)[1]<-"TXNAME"
head(tx2gene)
samples<-read.table(file="salmon.samples.txt", header = FALSE)
files<-file.path("L:","Lab-Sack","komudi","SequencingData","SequencingData_PBMCs","08-salmonV2",samples$V2,"quant.sf")
names(files)<-samples$V1 #add labels
txiquant.salmon <- tximport(files, type = "salmon", 
                       tx2gene = tx2gene, 
                       ignoreTxVersion = TRUE)
write.table(txiquant.salmon,file = "Salmon_geneCount_fastRefed_GROUPED_1_2v2.txt")

test2<-read.table(file = "L:/Lab-Sack/komudi/SequencingData/SequencingData_PBMCs/08-salmonV2/S2.quant/quant.sf", header = TRUE)
head(test1)
names(test1)[1]<-"TXNAME"

matchans<-match(test1$V1,tx2gene$TXNAME)
sum(is.na(matchans))
which(tx2gene == "ENST00000456328.2", arr.ind = TRUE)