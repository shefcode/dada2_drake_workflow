# Analyzing 16S Microbiome Insight data

source("packages.R")

# Pull metadata for the 11 samples for "filtered_bee_hind_gut_pellet" [Need to incorporate into drake]
load_meta <- read.csv("meta56_66.csv")
rownames(load_meta) <- load_meta$MBI.Sample.ID
meta <- select(load_meta, -c(MBI.Sample.ID))

# drake script
dada2analysis <- drake_plan(
  set_path = here("raw_sequences/samples/filtered_bee_hind_gut_pellet"),
  files = list.files(set_path, pattern = "fastq.gz"),
  fnFs = sort(list.files(set_path, pattern = "R1_001.fastq.gz", full.names = TRUE)),
  fnRs = sort(list.files(set_path, pattern = "R2_001.fastq.gz", full.names = TRUE)),
  sample.names = sapply(strsplit(basename(fnFs), "_"), `[`, 1),
  filtFs = file.path(set_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")),
  filtRs = file.path(set_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")),
  out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                      maxN=0, maxEE=c(2,30), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=FALSE),
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
  taxa = assignTaxonomy(seqtab.nochim, here("taxa/silva_nr_v138_train_set.fa.gz"), multithread=TRUE),
  taxa_sp = addSpecies(taxa, here("taxa/silva_species_assignment_v138.fa.gz")),
  taxa.print = set_rownames(taxa_sp, NULL),
  ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa)),
  dna = Biostrings::DNAStringSet(taxa_names(ps)),
  dna_n = setNames(dna, taxa_names(ps)),
  ps_m = merge_phyloseq(ps, dna_n),
  ps_t = setNames(taxa_names(ps_m), paste0("ASV", seq(ntaxa(ps_m))))
)

make(dada2analysis)

loadd(set_path, files, fnFs, fnRs, sample.names, filtFs, filtRs, out, errF, errR, derepFs, derepRs, dadaFs, dadaRs,
      mergers, seqtab, seqtab.nochim, getN, track, taxa, taxa_sp, taxa.print, ps, dna, dna_n, ps_m, ps_t)
vis_drake_graph(dada2analysis)

# Addition of phloseq
# following the DADA2 tutorial

theme_set(theme_bw())

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_data(meta), tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps) # dna_n
ps <- merge_phyloseq(ps, dna) # ps_m
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # ps_t
ps #

plot_richness(ps, measures=c("Chao1", "Shannon", "Simpson")) # Plot alpha diversity [Works]

# Ordinate
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")

# Bar Plot
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill = "Family")
                                    
########### astrobiomike tutorial  [Taxonomic Summaries] --- In Progress
phyla_counts_tab <- otu_table(tax_glom(ps, taxrank="Phylum")) 
phyla_counts_tab
phyla_tax_vec <- as.vector(tax_table(tax_glom(ps, taxrank="Phylum"))[,2]) 
colnames(phyla_counts_tab) <- as.vector(phyla_tax_vec)
unclassified_tax_counts <- colSums(seqtab.nochim) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)
phyla_and_unidentified_counts_tab
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]
class_counts_tab <- otu_table(tax_glom(ps, taxrank="Class")) 
class_tax_phy_tab <- tax_table(tax_glom(ps, taxrank="Class")) 
phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("Phylum"=phy_tmp_vec, "Class"=class_tmp_vec, row.names = rows_tmp)
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$Phylum == "Proteobacteria", "Class"])
colnames(class_counts_tab) <- as.vector(class_tax_tab$Class) 
proteo_class_counts_tab <- class_counts_tab[colnames(class_counts_tab) %in% proteo_classes_vec, ] 
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[colnames(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)
