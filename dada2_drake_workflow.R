# Drake script

source("package.R")
source("functions.R")

dada2analysis <- drake_plan(
  set_path = here("../data/"),
  files = list.files(set_path, pattern = ".fastq.gz"),
  fnFs = sort(list.files(set_path, pattern="1.fastq.gz", full.names = TRUE)),
  fnRs = sort(list.files(set_path, pattern="2.fastq.gz", full.names = TRUE)),
  sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1),
  filtFs = file.path(set_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")),
  filtRs = file.path(set_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")),
  out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,30), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE),
  errF = learnErrors(filtFs, multithread=TRUE),
  errR = learnErrors(filtRs, multithread=TRUE),
  derepFs = derepFastq(filtFs, verbose=TRUE) %>% 
    setNames(sample.names),
  derepRs = derepFastq(filtRs, verbose=TRUE) %>% 
    setNames(sample.names),
  dadaFs = dada(derepFs, err=errF, multithread=TRUE),
  dadaRs = dada(derepRs, err=errR, multithread=TRUE),
  mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE),
  seqtab = makeSequenceTable(mergers),
  seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE),
  getN = function(x) sum(getUniques(x)),
  track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) %>% 
    set_colnames(c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")) %>% 
    set_rownames(sample.names),
  taxa = assignTaxonomy(seqtab.nochim, here("../taxa/silva_nr_v138_train_set.fa.gz"), multithread=TRUE),
  taxa_sp = addSpecies(taxa, here("../taxa/silva_species_assignment_v138.fa.gz")),
  taxa.print = set_rownames(taxa_sp, NULL)
)

make(dada2analysis)
loadd(set_path, files, fnFs, fnRs, sample.names, filtFs, filtRs, out, errF, errR, derepFs, derepRs,dadaFs, dadaRs, mergers, seqtab, seqtab.nochim, getN, track, taxa, taxa_sp, taxa.print)

#### phyloseq?

library(phyloseq)

plot_richness(taxa)

# Analysis tools
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
head(out)
plotErrors(errF, nominalQ=TRUE)
dadaFs[[1]]
head(mergers[[1]])
dim(seqtab)
table(nchar(getSequences(seqtab)))
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
head(track)

vis_drake_graph(dada2analysis)