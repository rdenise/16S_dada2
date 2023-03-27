library(dada2)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(tidyr)


# Filename parsing

# Filtering
setwd("/data/san/data1/users/remi/16S_analysis")
pathF <- "results/databases/reads_trimmed/try_reads/F" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "results/databases/reads_trimmed/try_reads/R" # CHANGE ME ...
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              trimLeft = 0, truncLen=250, maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)


filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, "results/databases/reads_trimmed/try_reads/seqtab.rds") # CHANGE ME to where you want sequence table saved


st.all <- readRDS("results/databases/reads_trimmed/try_reads/seqtab.rds")
#####st.all <- mergeSequenceTables(st1, st2, st3)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/data/san/data0/databases/silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE) # https://zenodo.org/record/4587955
#taxa <- addSpecies(tax, "Documents/PhD/Databases/silva_species_assignment_v132.fa.gz")
# Write to disk
saveRDS(seqtab, "results/databases/reads_trimmed/try_reads/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "results/databases/reads_trimmed/try_reads/tax_final.rds") # CHANGE ME ...

seqtab <- readRDS("results/databases/reads_trimmed/try_reads/seqtab_final.rds")
tax   <- readRDS("results/databases/reads_trimmed/try_reads/tax_final.rds")
samples.out <- rownames(seqtab)
samdf <- data.frame(ID = samples.out)
row.names(samdf) = samples.out

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

count_table_tax  = t(rbind(ps@otu_table, t(ps@tax_table)))

count_table_one_word = data.frame(count_table_tax[,1:(nrow(ps@otu_table))])

View(head(count_table_one_word))

count_table_one_word$taxonomy = paste0(count_table_tax[,"Kingdom"], ";" , 
                                       count_table_tax[,"Phylum" ], ";" , 
                                       count_table_tax[,"Class"  ], ";" , 
                                       count_table_tax[,"Order"  ], ";" , 
                                       count_table_tax[,"Family" ], ";" , 
                                       count_table_tax[,"Genus"  ], ";")

write.csv(x = count_table_one_word, quote = F, file = "results/databases/reads_trimmed/try_reads/count_table_from_dada2.csv")

#Genus_count_table

agg <- aggregate_taxa(ps, level = "Genus",  verbose = TRUE) #Change me of you want a different level for your count table, Genus is reccomended.

colap <- unite(as.data.frame(agg@tax_table@.Data), newCol, -unique) 

genus_table <- agg@otu_table@.Data

View(head(genus_table))

row.names(genus_table) <- colap$newCol

write.csv(x = genus_table, quote = F, file = "results/databases/reads_trimmed/try_reads/genus_table_from_dada2.csv")





seqtab <- readRDS("results/databases/reads_trimmed/try_reads/seqtab_final.rds")
tax   <- readRDS("results/databases/reads_trimmed/try_reads/tax_final.rds")

samples.out <- rownames(seqtab)
samdf <- data.frame(ID = samples.out)
row.names(samdf) = samples.out

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

pip = t(ps@otu_table@.Data)

pip_otu_table = pip
row.names(pip_otu_table) = c(paste("otu", 1:nrow(pip_otu_table), sep = ""))

OTU_table = data.frame("OTU"=rownames(pip_otu_table),pip_otu_table)
OTU_table$taxonomy = paste(unlist(t(ps@tax_table@.Data)[5,]), unlist(t(ps@tax_table@.Data)[6,]), sep = "|")
#OTU_table = OTU_table[apply(OTU_table[,-c(1, ncol(OTU_table))] == 0, 1, sum) <= ((ncol(OTU_table) -2)  * 0.9), ] #Remove features with prevalence < 10%


write.table(OTU_table, 
            file = "results/databases/reads_trimmed/try_reads/PICRUSt2_otu_table_with_taxonomy.csv", row.names=FALSE, quote = F, sep = ",", )

write.table(OTU_table[,-ncol(OTU_table)], 
            file = "results/databases/reads_trimmed/try_reads/PICRUSt2_otu_table.csv", row.names=FALSE, quote = F, sep = ",", )


repseqs <- c(paste(">otu", 1:nrow(pip_otu_table), "\n", row.names(pip), "\n", sep = ""))
write.table(unname(c(paste(">otu", 
                           1:nrow(pip_otu_table), 
                           "\n", row.names(pip), 
                           sep = ""))), 
            file = "~/results/databases/reads_trimmed/try_reads/PICRUSt2_representative.csv", quote = F, row.names = F, col.names = F)