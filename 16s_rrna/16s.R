# Non-drake script

library(tidyverse)
library(here)
library(dada2)
library(drake)

set_path <- here("../data/")
files <- list.files(set_path, pattern = ".fastq")
fnFs <- sort(list.files(set_path, pattern="1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(set_path, pattern="2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
filtFs <- file.path(set_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(set_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,30), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
head(out)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Taxonomy

taxa <- assignTaxonomy(seqtab.nochim, here("../taxa/silva_nr_v138_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa, here("../taxa/silva_nr_v138_train_set.fa.gz"))
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
view(taxa.print)

# genoplotR 
# plotly
# upsetR



taxa_list <- as.data.frame(taxa.print)

aceto <- taxa_list %>% 
  filter(Family == "Acetobacteraceae")

colnames(taxa.print)
