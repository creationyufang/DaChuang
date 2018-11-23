######转录数据下载
library(SummarizedExperiment)
library(TCGAbiolinks)
#装了两个包

query.exp.hg38 <- GDCquery(project = "TCGA-LUAD", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM")


#创建了一个查询指令，括号内为条件限定
#“project”参数一共有43个选项，这里选了肺腺癌
#“data.category”参数一共有7个选项，这里选了转录组数据
#“data.type”参数这里一共有多少选项未知，这里选了基因表达量
#“workflow.type”参数这里一共有3个选项，分别是“HTSeq - Counts”“HTSeq - FPKM”"HTSeq - FPKM - UQ"
#（HTSeq:高通量测序）这里表示一个用来处理高通量测序数据的python包，可以得到基因水平的counts表达量
# 第一个选项是原始counts数 第二、第三的选项是normalized过后的量 有公式可以由原始数算出来
# 详情见这个网站 解释得很清楚 https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/ 
# ??????????问题：应该是normalized过后的比较好 但是我们还是不知道什么情况下用哪一种

GDCdownload(query.exp.hg38,files.per.chunk = 1）
#这一步开始下载满足上一步条件的数据。括号里分别为文件名和每次允许下载的文件数 
#"files.per.chunk"参数这里选定了1 意味着文件将一个一个下载。函数说明中说一次下载多个可在数据文件过大时减少下载时的问题
# ????????????问题：这里数值具体选什么？因为之前胃腺癌的数据比较小也选了50 这个文件300+M却只选了1
#下载时出现了报错 RStudio建议添加一个参数 method=‘client’ 添加之后就好了
#“method”这个参数一共两个选项 它默认的method是‘api’ 函数说明里描述api速度快 但是不稳定 很有可能要重复操作

exp.hg38 <- GDCprepare(query = query.exp.hg38)
# 将下载好的query转换成一个SummerizedExperiment的文件，这个以rda为后缀的文件是一个总结性文件，
# 有了它，我们可以不再需要之前下载的raw数据，所以后面的remove.files.prepared可以选择True，
# 这样会把之前下载的大量文件删除，当然也可以留着不删除（即default）。

exp.hg38.values <- assay(exp.hg38) 
#把数据转化成矩阵的形式
rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
#用文件中的“external_gene_name”定义矩阵的行名
write.csv(exp.hg38.values,file = "LUAD_exp_hg38.csv")
#输出文件为csv格式

#######获取肺腺癌临床数据及导出的命令
#下面的比较简单 就不解释了
clinical_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
colnames(clinical_LUAD)
write.table(clinical_LUAD,file="clinical_LUAD",sep = "\t",quote=T,col.names = T)
#临床数据整理
clinical_names=c("Tumor_Sample_Barcode","FAB_classification","days_to_last_followup","Overall_Survival_Status")

clinical_LUAD_m1 = DataFrame(clinical_LUAD$submitter_id,
                             clinical_LUAD$tumor_stage,
                             clinical_LUAD$days_to_last_follow_up,
                             clinical_LUAD$days_to_death)

############TCGA胃腺癌基因数据下载（4个转录组信息文件）
ibrary(TCGAbiolinks)
  library(SummarizedExperiment)
  query.exp.hg38 <- GDCquery(project = “TCGA-STAD”, 
  data.category = “Transcriptome Profiling”, 
  data.type = “Gene Expression Quantification”, 
  workflow.type = “HTSeq - FPKM”)
  GDCdownload(query.exp.hg38,files.per.chunk = 50)
  exp.hg38 <- GDCprepare(query = query.exp.hg38)
  exp.hg38.values <- assay(exp.hg38)
  rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
  write.csv(exp.hg38.values,file = “stad_exp_hg38_FPKM.csv”)
 #上面这一段与之前类似 不解释
 
  query.exp.hg38 <- GDCquery(project = “TCGA-STAD”, 
  data.category = “Transcriptome Profiling”, 
  data.type = “Gene Expression Quantification”, 
  workflow.type = “HTSeq - FPKM-UQ”)
  GDCdownload(query.exp.hg38,files.per.chunk = 50)
  exp.hg38 <- GDCprepare(query = query.exp.hg38)
  exp.hg38.values <- assay(exp.hg38)
  rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
  write.csv(exp.hg38.values,file = “stad_exp_hg38_FPKM_UQ.csv”)
 #上面这一段与之前类似 不解释
 
  query.exp.hg38 <- GDCquery(project = “TCGA-STAD”, 
  data.category = “Transcriptome Profiling”, 
  data.type = “Gene Expression Quantification”, 
  workflow.type = “HTSeq - Counts”)
  GDCdownload(query.exp.hg38,files.per.chunk = 50)
  exp.hg38 <- GDCprepare(query = query.exp.hg38)
  exp.hg38.values <- assay(exp.hg38)
  rownames(exp.hg38.values) <- values(exp.hg38)$external_gene_name
  write.csv(exp.hg38.values,file = “stad_exp_hg38_htseq_counts.csv”)
  #上面这一段与之前类似 不解释 
  
    query.exp <- GDCquery(project = “TCGA-STAD”, 
  legacy = TRUE,
   data.category = “Gene expression”,
   data.type = “Gene expression quantification”,
   platform = “Illumina HiSeq”, 
  file.type = “results”,
   experimental.strategy = “RNA-Seq”)
   #“legacy”参数的意思不知道是不是我查到的这个
   #“GDC Data Portal 中的数据是最新经过统一标准整理的，但有些数据还未开放，而 GDC Legacy Archive 中的数据是所有未经处理的数据，更全面”
   #“data.category”设置成了基因表达
   #“file.type”文件形式 但是不知道具体什么意思
  GDCdownload(query.exp,files.per.chunk = 50)
  stad.hg19 <- GDCprepare(query = query.exp)
  stad.hg19.exp <- assay(stad.hg19)
  rownames(stad.hg19.exp) <- values(stad.hg19)$gene_id
  write.csv(stad.hg19.exp,file = “stad_exp_hg19_rsem.csv”)
  #上面这一段与之前类似 不解释
