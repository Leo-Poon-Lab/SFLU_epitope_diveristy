library(tidyverse)
library(Biostrings)
library(readxl)
library(writexl)
library(ggpval)

# NP immunodominant epitope 
# amino acid position 366 - 374
# nucleotide position 1096 - 1122
nt_pos <- 1096:1122

# read metadata
df_meta <- read_excel("../data/sample_key_simple.xlsx")
df_meta <- df_meta %>% filter(Analysis=="Y")
df_meta$sample <- paste0(df_meta$CODE, "-NP")
df_meta$sample[grepl("^6E", df_meta$sample)] <- gsub("-NP", "", df_meta$sample[grepl("^6E", df_meta$sample)])
df_meta <- df_meta %>% arrange(TRANSFER,CODE)

# read the vcf data
files_NP_vcf_all <- list.files("../results/IRMA-results/", "A_NP.vcf", recursive = T)
files_NP_vcf <- files_NP_vcf_all[grepl("^6E", files_NP_vcf_all)]
files_NP_vcf <- c(files_NP_vcf, files_NP_vcf_all[grepl("-NP", files_NP_vcf_all)])
samples_all <- gsub("/.+", "", files_NP_vcf)

check_analysis <- samples_all %in% df_meta$sample
files_NP_vcf <- files_NP_vcf[check_analysis]
samples_all <- samples_all[check_analysis]

# estimate diversity at nucleotide level, and at haplotype level, considering every reads spanning the whole region.
system("epitope_diversity --help")
files_diversity <- sapply(files_NP_vcf, function(file_vcf_i){
	file_bam_i <- gsub("vcf", "bam", file_vcf_i)
	outfile_i <- gsub("vcf", "diversity", file_vcf_i)
	
	system(paste0("epitope_diversity -v -f ../results/IRMA-results/", file_bam_i, " -p ../data/np366.gff -o ../results/IRMA-results/", outfile_i))

	outfile_i
})

df_diversity <- lapply(files_diversity, function(x){
	tmp <- read_tsv(paste0("../results/IRMA-results/", x))
	tmp$sample <- gsub("/.+", "", x)
	tmp
})
df_diversity <- left_join(bind_rows(df_diversity), df_meta, "sample")
wilcox.test(df_diversity$Shannon_entropy[df_diversity$TRANSFER=="virus T cells"], df_diversity$Shannon_entropy[df_diversity$TRANSFER=="vaccine T cells"])
# t.test(df_diversity$Shannon_entropy[df_diversity$TRANSFER=="virus T cells"], df_diversity$Shannon_entropy[df_diversity$TRANSFER=="vaccine T cells"])
p <- ggplot(df_diversity, aes(x=TRANSFER, y=Shannon_entropy)) + geom_boxplot() + geom_jitter(alpha=0.8)+ylab("Shannon entropy (Haplotype)") + theme(axis.title.x = element_blank())
add_pval(p, pairs = list(c(2,3)), test='wilcox.test', alternative='two.sided')
ggsave("../results/boxplot_haplotype_Shannon_entropy.pdf", width = 6, height =6)

wilcox.test(df_diversity$population_nucleotide_diversity[df_diversity$TRANSFER=="virus T cells"], df_diversity$population_nucleotide_diversity[df_diversity$TRANSFER=="vaccine T cells"])
# t.test(df_diversity$population_nucleotide_diversity[df_diversity$TRANSFER=="virus T cells"], df_diversity$population_nucleotide_diversity[df_diversity$TRANSFER=="vaccine T cells"])
p <- ggplot(df_diversity, aes(x=TRANSFER, y=population_nucleotide_diversity)) + geom_boxplot() + geom_jitter(alpha=0.8)+ylab("Nucleotide diversity (Haplotype)") + theme(axis.title.x = element_blank())
add_pval(p, pairs = list(c(2,3)), test='wilcox.test', alternative='two.sided')
ggsave("../results/boxplot_haplotype_nucleotide_diversity.pdf", width = 6, height =6)

# write out the results in excel
df_diversity %>% write_xlsx("../results/data_haplotype.xlsx")
