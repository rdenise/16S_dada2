##########################################################################

library(dada2)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(tidyr)

##########################################################################

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log.file<-file(snakemake@log[[1]],open="wt")
sink(log.file)
sink(log.file, append=TRUE,type="message")

##########################################################################

cat("##########\n")
cat("# Taking care of the first two output\n")
cat("##########\n")

# Prepare arguments (no matter the order)
seqtab = readRDS(snakemake@input[["seqtab"]])
tax = readRDS(snakemake@input[["taxa_table"]])

##########################################################################

# Change the rownames of the sample data to match the seqtab
samples.out <- rownames(seqtab)

# Create a sample data frame
samdf <- data.frame(ID = samples.out)

# Set the rownames of the sample data frame to match the seqtab
row.names(samdf) = samples.out

# Create a phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

# Create a count table with taxonomy
count_table_tax  = t(rbind(ps@otu_table, t(ps@tax_table)))

# Create a count table with taxonomy and only one word per taxonomic level
count_table_one_word = data.frame(count_table_tax[,1:(nrow(ps@otu_table))])

# Add the taxonomy to the count table
count_table_one_word$taxonomy = paste0(count_table_tax[,"Kingdom"], ";" , 
                                       count_table_tax[,"Phylum" ], ";" , 
                                       count_table_tax[,"Class"  ], ";" , 
                                       count_table_tax[,"Order"  ], ";" , 
                                       count_table_tax[,"Family" ], ";" , 
                                       count_table_tax[,"Genus"  ], ";")

# Write the count table to a csv file
write.csv(x = count_table_one_word, 
        quote = F, 
        file = snakemake@output[["count_table"]])

# Genus_count_table
# Change me of you want a different level for your count table, Genus is reccomended.
agg <- aggregate_taxa(ps, level = "Genus",  verbose = TRUE) 

# Create a count table with taxonomy
colap <- unite(as.data.frame(agg@tax_table@.Data), newCol, -unique) 

# Extract the count table
genus_table <- agg@otu_table@.Data

# Add the taxonomy to the count table
row.names(genus_table) <- colap$newCol

# Write the count table to a csv file
write.csv(x = genus_table, 
        quote = F, 
        file = snakemake@output[["genus_table"]])


cat("##########\n")
cat("# Taking care of the three last output\n")
cat("##########\n")

# Change the rownames of the sample data to match the seqtab
samples.out <- rownames(seqtab)

# Create a sample data frame
samdf <- data.frame(ID = samples.out)

# Set the rownames of the sample data frame to match the seqtab
row.names(samdf) = samples.out

# Create a phyloseq object
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

# Create a count table with taxonomy
pip = t(ps@otu_table@.Data)

# Copy the taxonomy to the count table
pip_otu_table = pip

# Add the otu names to the count table
row.names(pip_otu_table) = c(paste("otu", 1:nrow(pip_otu_table), sep = ""))

# Create otu table with taxonomy
OTU_table = data.frame("OTU"=rownames(pip_otu_table),pip_otu_table)

# Add the taxonomy to the count table
OTU_table$taxonomy = paste(unlist(t(ps@tax_table@.Data)[5,]), unlist(t(ps@tax_table@.Data)[6,]), sep = "|")
# OTU_table = OTU_table[apply(OTU_table[,-c(1, ncol(OTU_table))] == 0, 1, sum) <= ((ncol(OTU_table) -2)  * 0.9), ] #Remove features with prevalence < 10%


write.table(OTU_table, 
            file = snakemake@output[["otu_table"]], 
            row.names=FALSE, quote = F, sep = ",")

write.table(OTU_table[,-ncol(OTU_table)], 
            file = snakemake@output[["otu_table_no_taxonomy"]], 
            row.names=FALSE, quote = F, sep = ",")

# Create a fasta file with the representative sequences
repseqs <- c(paste(">otu", 1:nrow(pip_otu_table), "\n", row.names(pip), "\n", sep = ""))

write.table(unname(c(paste(">otu", 
                           1:nrow(pip_otu_table), 
                           "\n", row.names(pip), 
                           sep = ""))), 
            file = snakemake@output[["repseqs"]], 
            quote = F, row.names = F, col.names = F)

# Proper syntax to close the connection for the log file
# but could be optional for Snakemake wrapper
sink(type="message")
sink()