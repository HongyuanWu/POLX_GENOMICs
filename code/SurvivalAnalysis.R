# load library
suppressMessages(library("survival"))
suppressMessages(library("papeR"))
suppressMessages(library("cowplot"))
suppressMessages(library("VennDiagram"))
suppressMessages(library("grid"))
source("function.R")

# load data 
options(stringsAsFactors = F)
mut.dat<-read.delim("../data/ICI_data_mutations_extended.txt",header=1)
sample.dat<-read.delim("../data/ICI_sample_clinical.tsv",header=1)
pat.dat<-read.delim("../data/ICI_Patient_clinica.tsv",header=1)
clinical.dat<-merge(sample.dat,pat.dat,by="PATIENT_ID")



get_gene_clinical_matrix<-function(gene,mut.dat,clinical.dat){
  sub.mut.data<-mut.dat[mut.dat$Hugo_Symbol%in%focus.gene,]
  pos.sample<-sub.mut.data$Tumor_Sample_Barcode
  pos.patient.id<-substr(pos.sample, 1, nchar(pos.sample) -8)
  tempdata<-clinical.dat
  # Refine MSI
  tempdata$MSI_TYPE[tempdata$MSI_TYPE%in%c("Do not report","Not Available")]=NA
  tempdata$MSI_TYPE[tempdata$MSI_TYPE%in%c("Indeterminate","Stable")]="MSI-L/stable"
  # get status whether has the mutation of specific genes
  tempdata$GeneStatus<-as.factor(ifelse(tempdata$PATIENT_ID%in%pos.patient.id,"Mut","WT"))
  # remove unmatched samples 
  tempdata<-tempdata[tempdata$SOMATIC_STATUS!="Unmatched",]
  # top 20% TMB sample  of each tumor type (adjust for tumor type) 
  split.cancer <- split(tempdata, tempdata$CANCER_TYPE)
  topnames<-lapply(split.cancer, function(x) {return(top_n(x, round(nrow(x)*0.2), TMB_SCORE)$PATIENT_ID)})
  top.names<-unlist(topnames,use.names = F) 
  tempdata$TMB_binary<-ifelse(tempdata$PATIENT_ID %in% top.names,"TMB-H","TMB-L")
  return(tempdata)
  
}

# Single Gene survival analysis ----
focus.gene<-c("POLD1")
plot.data<-get_gene_clinical_matrix(focus.gene,mut.dat,clinical.dat)
surv<-Surv(plot.data$OS_MONTHS,plot.data$OS_STATUS=="DECEASED")
fit<- survfit(surv ~ plot.data$GeneStatus)
p1<-SurvivalCurvePlot(fit, plot.data, pval = T, genename=focus.gene,panel = pal_aaas()(2))

focus.gene<-c("POLE")
plot.data<-get_gene_clinical_matrix(focus.gene,mut.dat,clinical.dat)
surv<-Surv(plot.data$OS_MONTHS,plot.data$OS_STATUS=="DECEASED")
fit<- survfit(surv ~ plot.data$GeneStatus)
p2<-SurvivalCurvePlot(fit, plot.data, pval = T, genename=focus.gene,panel = pal_aaas()(2))

# 
focus.gene<-c("POLE","POLD1")
plot.data<-get_gene_clinical_matrix(focus.gene,mut.dat,clinical.dat)
surv<-Surv(plot.data$OS_MONTHS,plot.data$OS_STATUS=="DECEASED")
fit<- survfit(surv ~ plot.data$GeneStatus)
p3<-SurvivalCurvePlot(fit, plot.data, pval = T, genename="POLE/POLD1",panel = pal_aaas()(2))

temp.list<-list(p1,p2,p3)
survplotlist<-list()
for(i in 1:length(temp.list)){
  surv_plot<-temp.list[[i]]
  survplotlist[[i]] = plot_grid(surv_plot$plot, surv_plot$table, ncol = 1, align = 'v',rel_heights = c(2/3, 1/3))
}

pdf(file="../results/Figure_1DEF.pdf",width = 20,height = 7)
cowplot::plot_grid(plotlist = survplotlist,nrow=1)
dev.off()


# ven plot for calculation overlap be MSS-H and POLE/POLD mutation----------
table(plot.data$MSI_TYPE,plot.data$GeneStatus)
plot.ven.data<-plot.data[,c("PATIENT_ID","MSI_TYPE","GeneStatus")]
plot.ven.data<-plot.data[complete.cases(plot.ven.data),]
list.all<-plot.ven.data$PATIENT_ID
list.msi<-plot.ven.data$PATIENT_ID[which(plot.ven.data$MSI_TYPE%in%"Instable")]
list.mutation<-plot.ven.data$PATIENT_ID[plot.ven.data$GeneStatus=="Mut"]
pdf(file="../results/Figure_1H.pdf",width =8,height = 8)
tempven<-venn.diagram(list(`ALL sample with MSS status`=list.all,
                            `MSI-H`=list.msi,
                            `POLE/POLD1 mutated sample`=list.mutation),
                       fill = pal_lancet()(3),
                       col = "transparent",
                       alpha = 0.5,
                       cex = 2.5,
                       filename = NULL)
grid_draw(tempven)
dev.off()



# Survival Analysis of  POLE/POLD1 Mutation and MSI status----
focus.gene<-c("POLE","POLD1")
plot.data<-get_gene_clinical_matrix(focus.gene,mut.dat,clinical.dat)
plot.data.tmb.POLE<-plot.data[,c("TMB_SCORE","GeneStatus","OS_MONTHS","OS_STATUS","PATIENT_ID","MSI_TYPE","TMB_binary")]
plot.data.tmb.POLE<-plot.data.tmb.POLE[!is.na(plot.data.tmb.POLE$MSI_TYPE),]
plot.data.tmb.POLE$TMB_POLX<-NA

for(i in 1:nrow(plot.data.tmb.POLE )){
    if(plot.data.tmb.POLE$MSI_TYPE[i]=="Instable"){
        plot.data.tmb.POLE$TMB_POLX[i] <- "Instable"
    }else if((plot.data.tmb.POLE$GeneStatus[i]=="Mut")&&(plot.data.tmb.POLE$MSI_TYPE[i]=="MSI-L/stable")){
        plot.data.tmb.POLE$TMB_POLX[i] <- "Mut/MSI-L"
    }else if((plot.data.tmb.POLE$GeneStatus[i]=="WT")&&(plot.data.tmb.POLE$MSI_TYPE[i]=="MSI-L/stable")){
        plot.data.tmb.POLE$TMB_POLX[i] <- "WT/MSI-L"
    }
}

# check data
table(plot.data.tmb.POLE$TMB_POLX)
# set orders
plot.data.tmb.POLE$TMB_POLX<-factor(plot.data.tmb.POLE$TMB_POLX,levels = c("Instable","WT/MSI-L","Mut/MSI-L"))
surv<-Surv(plot.data.tmb.POLE$OS_MONTHS,plot.data.tmb.POLE$OS_STATUS=="DECEASED")
fit<- survfit(surv ~ plot.data.tmb.POLE$TMB_POLX)
p5<-SurvivalCurvePlot(fit, plot.data.tmb.POLE, pval = T, genename="POLE/POLD1 mutated in Non-MSI-H subgroup",panel = pal_aaas()(3),
                      label=c("MSI-H","Mut+nonMSI-H","WT+nonMSi-H"))
# p5
dev.off()   
pdf(file="../results/Figure_1H.pdf",width =10,height = 8)
print(p5)
dev.off()


#Multivariable Cox regression data analysis ----
#interest varable 
focus.gene<-c("POLE","POLD1")
# int.variable<-c("CANCER_TYPE","MSI_TYPE","SEX","GeneStatus","TMB_SCORE","AGE_AT_SEQ_REPORT")
int.variable<-c("CANCER_TYPE","MSI_TYPE","GeneStatus")
plot.data2<-get_gene_clinical_matrix(focus.gene,mut.dat,clinical.dat)
plot.data2<-plot.data2[,c(int.variable,"OS_MONTHS","OS_STATUS")]
plot.data2<-plot.data2[complete.cases(plot.data2),]



# unicox regression 
surv<-Surv(plot.data2$OS_MONTHS,plot.data2$OS_STATUS=="DECEASED")
unicox.dat<-batchUnivarCOXfun(surv,plot.data2[,int.variable])
filtered.variable<-row.names(unicox.dat[unicox.dat[,7]<0.05,])

# multiple regression analysis
fmla <- as.formula(paste("surv ~ ", paste(filtered.variable, collapse= "+")))
cox.fit <- coxph(fmla, data=plot.data2)
cox.fit2.write <- prettify(summary(cox.fit),digits=4)
write.csv(cox.fit2.write,file="../results/Cox_regression_multivariable_analysis.csv")
# survival analysis