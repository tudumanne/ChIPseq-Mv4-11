#Reference modified from https://github.com/mevers/frag_align_rDNA

#load libraries
library(Rsamtools)
library(tidyverse)

#read in data - rDNA reference and BAM file
rDNA <- readDNAStringSet("KY962518.fa")
bam <- BamFile("unmasked_150bp_overlap_k1000_sorted.bam")

#generate rDNA to genome read map based on BAM reads
df_map <- bam %>%
  scanBam() %>%
  pluck(1) %>%
  keep(names(.) %in% c("qname", "rname", "pos", "qwidth")) %>%
  as_tibble() %>%
  separate(qname, c("tmp", "from", "to", "tmp1", "tmp2"), sep = ";") %>%
  select(-tmp, -tmp2, -tmp2) %>%
  mutate(from_chr = "rDNA", from_start = str_remove(from, "KY962518.1_start="), from_end = str_remove(to, "end=")) %>%
  rename(to_chr = rname, to_start = pos, to_end = qwidth) %>%
  mutate(to_end = to_end + to_start) %>%
  select(from_chr, from_start, from_end, to_chr, to_start, to_end) %>%
  as.data.frame()
df_map$from_start=as.integer(df_map$from_start)
df_map$from_end=as.integer(df_map$from_end)

#export the mapped regions as a .bed file
df_map_table = as.data.frame(df_map)
df_map_table$to_chr = as.numeric(df_map_table$to_chr)

bed <- df_map_table[,c('to_chr', 'to_start', 'to_end')]
colnames(bed) <- c('chrom', 'chromStart', 'chromEnd')
write.table(bed, "regions.bed")
