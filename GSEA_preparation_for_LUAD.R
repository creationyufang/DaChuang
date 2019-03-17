#用R处理下，将癌组织和对应的癌旁组织的数据分别提取出来分别作为两组的表达矩阵
#以及或者分组文件（cls文件）
#从上述代码，我获得tumor_number个癌组织样本和对应的normal_number个癌旁样本的表达谱数据，并且将Ensembl ID均转化为了Gene symbol
#（避免之后用GSEA时，再用chip做ID转化）；然后可以直接将txt文件作为输入
library(dplyr)
library(stringr)

data <- read.table(file = "TCGA-LUAD.htseq_fpkm.tsv", sep = "\t", header = T, stringsAsFactors = F, check.names = F)

data_df <- data[, -1]
#str_sub(string, start = 1L, end = -1L) (in package stringr)
#in TCGA, 1~9 represent for tumor, 10~29 stand for normal
#"%in%" <- function(x, table) match(x, table, nomatch = 0) > 0
normal_df <- data_df[, str_sub(names(data_df), 14, 15) %in% 11:19]
normal_number <- length(names(normal_df))

tumor_df <- data_df[, grepl(paste(str_sub(names(normal_df), 1, 12), collapse = "|"), names(data_df)) & !names(data_df) %in% names(normal_df)]
tumor_number <- length(names(tumor_df))

idmapping <- read.table(file = "gencode.v22.annotation.gene.probeMap", sep = "\t", header = T, stringsAsFactors = F)
geneid <- data.frame(id = data$Ensembl_ID, stringsAsFactors = F)
geneid2symbol <- left_join(geneid, idmapping, by = "id")

all_df <- cbind(Name = geneid2symbol$gene, DESCRIPTION = "na", tumor_df, normal_df)

group <- c(rep("Tumor", tumor_number), rep("Normal", normal_number))
group <- paste(group, collapse = " ")
group <- c(paste(c(normal_number+tumor_number, 2, 1), collapse = " "), "# Tumor Normal", group)

write.table(file = "LUAD_fpkm.txt", all_df, sep = "\t", col.names = T, row.names = F, quote = F)
write.table(file = "group_LUAD.cls", group, col.names = F, row.names = F, quote = F)
