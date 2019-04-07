#reference:https://blog.csdn.net/weixin_40466280/article/details/80377317
#生存分析与Rhttps://blog.csdn.net/weixin_40466280/article/details/80377317

#?how to change data/dataframe to values:
#w=lung[2,]
#w=w[1,]

#?how to get specific rows in a dataframe
#reference:https://blog.csdn.net/gpwner/article/details/69048465
temp=data.frame(a=c(1,2,3,4),b=c("a","b","c","d"),d=c(6,7,8,9))
temp[temp$b %in% c("b","c"),]

#?how to delete duplicated ones
a <- a[!duplicated(a)]

#1.library packages
library(survival)
library(survminer)

#2.inport data & select specific columns
 #2.1 import_data
lusc=read.csv("lusc_exp_hg38_htseq_FPKM_uq.csv",header = T)
colnames(lusc)[1]="genesymble"
clinical=read.table("clinical.tsv",header = T,sep="\t")
 #2.2 motify_colnames_of_lusc_and_pick_out_status/timetodeath_of_clinical
temp=colnames(lusc)
temp=chartr(".","-",temp)
colnames(lusc)=temp

submitter_id=substr(temp[-1],1,12)
lusc_temp=lusc
colnames(lusc_temp)=c("genesymble",submitter_id)
sign=substr(temp[-1],14,15)
head(sign)
lusc_cancer=lusc_temp[,c(TRUE,sign=="01")]
rm(lusc_temp)

clinical=clinical[clinical$project_id=="TCGA-LUSC",c("submitter_id","vital_status","days_to_death")]

 #2.3 select_interset_patientID_&_creat_dataframe
id_common=intersect(clinical$submitter_id,colnames(lusc_cancer))
lusc_=lusc_cancer[,c("genesymble",id_common)]
clinical_=clinical[clinical$submitter_id %in% id_common,]
clinical=clinical_
lusc=lusc_

 #2.4 transfer_alive/dead_to_1/2
temp=clinical$vital_status=="alive"
clinical$vital_status[temp]=1
clinical$vital_status[!temp]=2
rm(temp)

#View(head(clinical_))
#View(head(lusc_))
#write.csv(lusc_,file="lusc_.csv",col.names = T,row.names = F,sep = ",",quote=F)
#write.csv(clinical_,file="clinical_.csv",row.names = F,quote=F)
#lusc=read.csv("lusc_.csv",header =T,stringsAsFactors = F)
#clinical=read.csv("clinical_.csv",header=T,stringsAsFactors = F)
 
 #2.4 (optional)combine_expression_and_clinical_data
#temp=clinical_$vital_status
#temp=as.numeric(temp)
#temp=c(colnames(clinical_)[2],temp)
#temp2=clinical_$days_to_death
#temp2=as.numeric(temp2)
#temp2=c(colnames(clinical_)[3],temp2)
#data=rbind(temp,temp2,lusc_)

#3. determine target genes
genelist=c("CCNB2","PLK1","RRM2","UBE2C")
genelist=c("DLC1","CSDC2","COL13A1","CLDNS","CLIP")
#,"CFD","CASS4","CAPN8","C14orf180","ARHGAP44","ANK2","ACSM5"

##4. create cox object
clinical$vital_status=as.numeric(clinical$vital_status)
# 1 for alive, 2 for dead
clinical$days_to_death=as.numeric(clinical$days_to_death)

fit.surv=Surv(clinical$days_to_death,clinical$vital_status)
expr=lusc[genelist,]
#expr=lusc[lusc$genesymble %in% genelist,]
colnames(expr)=chartr(".","-",colnames(expr))

#rownames(expr)=expr[,1]
#expr=expr[,-1]
expr=t(expr)
expr=data.frame(expr)

##expr=data.frame(as.numeric(as.character(expr[,1])),as.numeric(as.character(expr[,2])),as.numeric(as.character(expr[,3])))
##colnames(expr)=genelist[-4]
head(expr)
expr=round(expr)
#temp=expr$ARHGAP44
res.cox=coxph(Surv(clinical$days_to_death,clinical$vital_status)~expr$C1orf112+expr$SCYL3+expr$KEAP1+expr$MTMR7+expr$HECW1)
res.cox
cox.zph(res.cox)
temp=cox.zph(res.cox)
ggcoxzph(temp)

#5. Single Factor Cox
lusc<- lusc[!duplicated(lusc[,1]),]
rownames(lusc)=lusc[,1]
lusc=lusc[,-1]

z=c()
coefficients=c()
score=c()
geneid=c()
numb=c()
HR=c()

#length(lusc[,1])
for(i in length(lusc[,1])){
  temp=lusc[i,]
  temp=t(temp)
  temp=as.vector(temp)
  data=data.frame(temp,clinical$days_to_death,clinical$vital_status)
  res=coxph(Surv(data$clinical.days_to_death,data$clinical.vital_status)~data$temp)
  geneid[i]=rownames(lusc)[i]
  #wald=as.numeric(res$wald.test)
  #z[i]=wald
  coe=as.numeric(res$coefficients)
  coefficients[i]=coe
  HR[i]=exp(coefficients[i])
  #score=as.numeric(res$score)
  #score[i]=score
  numb[i]=i
}
single_cox=data.frame(numb,geneid,coefficients,HR)
single_cox=single_cox[order(single_cox$HR,decreasing = T),]
write.csv(single_cox,file="single_cox.csv",row.names = F,quote = F)
