# Please note that the IRMA can only be successfully started on Linux.

library(tidyverse)
library(readxl)

# 2022-09
# data <- read_excel("../data/sample key - NGS - Updated.xlsx")
# data <- data[data[,4]=="Y",]

files_fq_1 <- list.files("..//data/Data files for HK/Raw data/", "R1_001")
files_fq_2 <- gsub("_R1_", "_R2_", files_fq_1)
setwd("../data/Data files for HK/Raw data")

sapply(seq_along(files_fq_1), function(i){
# sapply(1:9, function(i){
	# i=34
	files_fqgz_1_i <- files_fq_1[i]
	files_fqgz_2_i <- files_fq_2[i]
	sample_name <- gsub("_S\\d+_L001_R1_001.fastq.gz$", "", files_fqgz_1_i)
	sample_name <- gsub(".+//", "", sample_name)
	# system(paste0("gzip -d ", files_fqgz_1_i))
	# system(paste0("gzip -d ", files_fqgz_2_i))

	# files_fq_1_i <- gsub("\\.gz", "", files_fqgz_1_i)
	# files_fq_2_i <- gsub("\\.gz", "", files_fqgz_2_i)

	system(paste0("IRMA FLU '", files_fqgz_1_i, "' '", files_fqgz_2_i, "' ", sample_name))
	# system(paste0("IRMA FLU '", files_fq_1_i, "' '", files_fq_2_i, "' ", sample_name))

	# system(paste0("rm ", files_fq_1_i))
	# system(paste0("rm ", files_fq_2_i))

})

results_dirs <- list.dirs(recursive = F)

sapply(results_dirs, function(x) {
	dest_path <- gsub("\\.", "../../../results/", x)
	system(paste0("mv ", x, " ", dest_path))
})

# 2023-05-29

# 2022-09
# data <- read_excel("../data/sample key - NGS - Updated.xlsx")
# data <- data[data[,4]=="Y",]

files_fq_1 <- list.files("..//data/2023 Yi-Mo Raw data/", "R1_001")
files_fq_2 <- gsub("_R1_", "_R2_", files_fq_1)
setwd("../data/2023 Yi-Mo Raw data")

sapply(seq_along(files_fq_1), function(i){
# sapply(1:9, function(i){
	# i=1
	files_fqgz_1_i <- files_fq_1[i]
	files_fqgz_2_i <- files_fq_2[i]
	sample_name <- gsub("_S\\d+_L001_R1_001.fastq.gz$", "", files_fqgz_1_i)
	sample_name <- gsub(".+//", "", sample_name)
	# system(paste0("gzip -d ", files_fqgz_1_i))
	# system(paste0("gzip -d ", files_fqgz_2_i))

	# files_fq_1_i <- gsub("\\.gz", "", files_fqgz_1_i)
	# files_fq_2_i <- gsub("\\.gz", "", files_fqgz_2_i)

	system(paste0("IRMA FLU '", files_fqgz_1_i, "' '", files_fqgz_2_i, "' ", sample_name))
	# system(paste0("IRMA FLU '", files_fq_1_i, "' '", files_fq_2_i, "' ", sample_name))

	# system(paste0("rm ", files_fq_1_i))
	# system(paste0("rm ", files_fq_2_i))

})

results_dirs <- list.dirs(recursive = F)

sapply(results_dirs, function(x) {
	dest_path <- gsub("\\.", "../../results/", x)
	system(paste0("mv ", x, " ", dest_path))
})
