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
files_NP_vcf <- list.files("../results/IRMA-results/", "A_NP.vcf", recursive = T)
files_NP_vcf <- files_NP_vcf_all[grepl("^6E", files_NP_vcf_all)]
files_NP_vcf <- c(files_NP_vcf, files_NP_vcf_all[grepl("-NP", files_NP_vcf_all)])
samples_all <- gsub("/.+", "", files_NP_vcf)

check_analysis <- samples_all %in% df_meta$sample
files_NP_vcf <- files_NP_vcf[check_analysis]
samples_all <- samples_all[check_analysis]

# df_vcf <- lapply(files_NP_vcf, function(x){
# 	# x <- files_NP_vcf[1]
# 	full_path_i <- paste0("../results/IRMA-results/", x)
# 	tmp <- read_tsv(full_path_i,comment = "##", col_types = cols(.default = "c"))
# 	tmp$sample <- gsub("/.+", "", x)
# 	tmp
# })
# df_vcf <- bind_rows(df_vcf)
# df_vcf_epitope <- df_vcf %>% mutate(POS=as.numeric(POS)) %>% filter(POS<=1122 & POS>=1096)
# unique(df_vcf_epitope$sample)

# read all alleles data
files_NP_allAlleles_all <- list.files("../results/IRMA-results/", "A_NP-allAlleles.txt", recursive = T)
files_NP_allAlleles <- files_NP_allAlleles_all[grepl("^6E", files_NP_allAlleles_all)]
files_NP_allAlleles <- c(files_NP_allAlleles, files_NP_allAlleles_all[grepl("-NP", files_NP_allAlleles_all)])
samples_all <- gsub("/.+", "", files_NP_allAlleles)

check_analysis <- samples_all %in% df_meta$sample
files_NP_allAlleles <- files_NP_allAlleles[check_analysis]
samples_all <- samples_all[check_analysis]

df_allAlleles <- lapply(files_NP_allAlleles, function(x){
	# x <- files_NP_allAlleles[1]
	full_path_i <- paste0("../results/IRMA-results/", x)
	tmp <- read_tsv(full_path_i,comment = "##", col_types = cols(.default = "c"))
	tmp$sample <- gsub("/.+", "", x)
	tmp
})
df_allAlleles <- bind_rows(df_allAlleles)
df_allAlleles_epitope <- df_allAlleles %>% mutate(POS=as.numeric(Position)) %>% filter(POS<=1122 & POS>=1096)

# estimate diversity at nucleotide level, and at locus level, then take average within a region.
## Reference: https://academic.oup.com/ve/article/5/1/vey041/5304643
## Shannon entropy
df_meta$shannon_entropy_locus <- sapply(df_meta$sample, function(sample_i){
	df_allAlleles_tmp <- df_allAlleles_epitope %>% filter(sample==sample_i)
	if(nrow(df_allAlleles_tmp)==0) {return(0)} else {
		hs_locus_i <- sapply(nt_pos, function(nt_pos_i){
			if(!nt_pos_i %in% df_allAlleles_tmp$POS) {return(0)} else {
				p_locus_i <- df_allAlleles_tmp %>% filter(POS==nt_pos_i) %>% .$Frequency %>% as.numeric()
				stopifnot(sum(p_locus_i)>0.99)
				hs_locus_i <- -sum(p_locus_i*log2(p_locus_i))
			}
		})
		mean(hs_locus_i)
	}
})
wilcox.test(df_meta$shannon_entropy_locus[df_meta$TRANSFER=="virus T cells"], df_meta$shannon_entropy_locus[df_meta$TRANSFER=="vaccine T cells"])
# t.test(df_meta$shannon_entropy_locus[df_meta$TRANSFER=="virus T cells"], df_meta$shannon_entropy_locus[df_meta$TRANSFER=="vaccine T cells"])
p <- ggplot(df_meta, aes(x=TRANSFER, y=shannon_entropy_locus)) + geom_boxplot() + geom_jitter(alpha=0.8)+ylab("Shannon entropy (Locus)") + theme(axis.title.x = element_blank())
add_pval(p, pairs = list(c(2,3)), test='wilcox.test', alternative='two.sided')
ggsave("../results/boxplot_locus_Shannon_entropy.pdf", width = 6, height =6)

## nucleotide diversity
df_meta$nucleotide_diversity_locus <- sapply(df_meta$sample, function(sample_i){
	df_allAlleles_tmp <- df_allAlleles_epitope %>% filter(sample==sample_i)
	if(nrow(df_allAlleles_tmp)==0) {return(0)} else {
		hpi_locus_i <- sapply(nt_pos, function(nt_pos_i){
			if(!nt_pos_i %in% df_allAlleles_tmp$POS) {return(0)} else {
				p_locus_i <- df_allAlleles_tmp %>% filter(POS==nt_pos_i) %>% .$Frequency %>% as.numeric()
				stopifnot(sum(p_locus_i)>0.99)
				
				dp_locus_i <- df_allAlleles_tmp %>% filter(POS==nt_pos_i) %>% .$Total %>% .[1] %>% as.numeric()
				freq_locus_i = dp_locus_i * p_locus_i

				hpi_locus_i <- 1-((sum(freq_locus_i*(freq_locus_i-1))/(dp_locus_i*(dp_locus_i-1))))
			}
		})
		mean(hpi_locus_i)
	}
})
wilcox.test(df_meta$nucleotide_diversity_locus[df_meta$TRANSFER=="virus T cells"], df_meta$nucleotide_diversity_locus[df_meta$TRANSFER=="vaccine T cells"])
# t.test(df_meta$nucleotide_diversity_locus[df_meta$TRANSFER=="virus T cells"], df_meta$nucleotide_diversity_locus[df_meta$TRANSFER=="vaccine T cells"])
p <- ggplot(df_meta, aes(x=TRANSFER, y=nucleotide_diversity_locus)) + geom_boxplot() + geom_jitter(alpha=0.8)+ylab("Nucleotide diversity (Locus)") + theme(axis.title.x = element_blank())
add_pval(p, pairs = list(c(2,3)), test='wilcox.test', alternative='two.sided')
ggsave("../results/boxplot_locus_nucleotide_diversity.pdf", width = 6, height =6)

## Fst
## nucleotide diversity between groups, and within groups.
pairs_all <- combn(df_meta$sample[df_meta$TRANSFER!="NO T CELLS"], 2)
df_pairs_all <- as_tibble(t(pairs_all))
names(df_pairs_all) <- c("sample1","sample2")
df_pairs_all$group1 <- factor(df_pairs_all$sample1, df_meta$sample[df_meta$TRANSFER!="NO T CELLS"], labels=df_meta$TRANSFER[df_meta$TRANSFER!="NO T CELLS"])
df_pairs_all$group2 <- factor(df_pairs_all$sample2, df_meta$sample[df_meta$TRANSFER!="NO T CELLS"], labels=df_meta$TRANSFER[df_meta$TRANSFER!="NO T CELLS"])
df_pairs_all$type <- ifelse(df_pairs_all$group1==df_pairs_all$group2, "within", "between")

df_pairs_all$nucleotide_diversity_locus <- sapply(seq_len(nrow(df_pairs_all)), function(i){
	# i=1
	sample1_i <- df_pairs_all$sample1[i]
	sample2_i <- df_pairs_all$sample2[i]
	df_allAlleles_tmp <- df_allAlleles_epitope %>% filter(sample %in% c(sample1_i, sample2_i))

	if(nrow(df_allAlleles_tmp)==0) {return(0)} else {
		hpi_locus_i <- sapply(nt_pos, function(nt_pos_i){
			if(!nt_pos_i %in% df_allAlleles_tmp$POS) {return(0)} else {
				p_locus_i <- df_allAlleles_tmp %>% filter(POS==nt_pos_i) %>% .$Frequency %>% as.numeric()
				# print(p_locus_i)
				stopifnot(sum(p_locus_i)>1.99)
				
				dp_locus_i <- df_allAlleles_tmp %>% filter(POS==nt_pos_i) %>% .$Total %>% as.numeric()
				freq_locus_i <- dp_locus_i * p_locus_i
				dp_sum_locus_i <- sum(unique(dp_locus_i))
				
				stopifnot(abs(sum(freq_locus_i)-dp_sum_locus_i)<1)

				hpi_locus_i <- 1-((sum(freq_locus_i*(freq_locus_i-1))/(dp_sum_locus_i*(dp_sum_locus_i-1))))
			}
		})
		mean(hpi_locus_i)
	}
})
wilcox.test(df_pairs_all$nucleotide_diversity_locus[df_pairs_all$type=="within"], df_pairs_all$nucleotide_diversity_locus[df_pairs_all$type=="between"])
# t.test(df_pairs_all$nucleotide_diversity_locus[df_pairs_all$type=="within"], df_pairs_all$nucleotide_diversity_locus[df_pairs_all$type=="between"])
p <- ggplot(df_pairs_all, aes(x=type, y=nucleotide_diversity_locus)) + geom_boxplot() + geom_jitter(alpha=0.8)+ylab("Nucleotide diversity (Locus)") + theme(axis.title.x = element_blank())
add_pval(p, pairs = list(c(1,2)), test='wilcox.test', alternative='two.sided')
ggsave("../results/boxplot_locus_fst.pdf", width = 6, height =6)

# write out the results in excel
df_meta %>% write_xlsx("../results/data_locus.xlsx")
df_pairs_all %>% write_xlsx("../results/data_locus_fst.xlsx")
