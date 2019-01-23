suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(ggthemes))
source("function.R")


#-------------------------------------------------------------
# plot general mutation frequencies of POLE or POLD1 genes with TMB data #
#------------------------------------------------------------
# Load Data----

bioportal.data<-read.delim("../data/bioportalData.tsv",header = T)
bioportal.data.agre<-aggregate(bioportal.data[,2:6],list(bioportal.data$GeneralCancerType),sum)
colnames(bioportal.data.agre)[1]="tumor_type"
bioportal.data.agre$tumor_type<-as.character(bioportal.data.agre$tumor_type)
# get pole total rate 
sum(bioportal.data.agre$POLE+bioportal.data.agre$POLE.POLD1)/sum(bioportal.data.agre$Total)
# get pold total rate 
sum(bioportal.data.agre$POLD1+bioportal.data.agre$POLE.POLD1)/sum(bioportal.data.agre$Total)
#write.csv(bioportal.data.agre,file = "preprocessed.table.csv")


# get cancer type with sample number  more than 30
bioportal.data.agre<-bioportal.data.agre[bioportal.data.agre$Total>30,]
# get cancer type with pole or pold mutations more than 5 
bioportal.data.agre<-bioportal.data.agre[(bioportal.data.agre$POLE+bioportal.data.agre$POLE.POLD1+bioportal.data.agre$POLD1)>5,]
bioportal.data.agre$totalFre<-round((bioportal.data.agre$POLE+bioportal.data.agre$POLE.POLD1+bioportal.data.agre$POLD1)*100/bioportal.data.agre$Total,digits = 2)
bioportal.data.agre<-bioportal.data.agre[order(bioportal.data.agre$totalFre,decreasing = T),]
bioportal.data.agre$allMutated<-bioportal.data.agre$POLE+bioportal.data.agre$POLE.POLD1+bioportal.data.agre$POLD1
bioportal.data.agre$POLD1<-round(bioportal.data.agre$POLD1*100/bioportal.data.agre$Total,digits = 2)
bioportal.data.agre$POLE<-round(bioportal.data.agre$POLE*100/bioportal.data.agre$Total,digits = 2)
bioportal.data.agre$POLE.POLD1<-round(bioportal.data.agre$POLE.POLD1*100/bioportal.data.agre$Total,digits = 2)

bp1 <-bar_plot_function(bioportal.data.agre)


pdf("../results/Figure_1B.pdf",width = 12,height = 7)
print(bp1)
dev.off()


#----------------------------------------#
# plot POLE/POLD1 mutation with TMB data #
#---------------------------------------#
#Load data ----
msk.tmb.dat<-read.delim("../data/MSK_tmb_data.tsv",header=1)
#head(msk.tmb.dat)
MSKCC.panel.size<-1.22# 1.22MB for this panel 
msk.tmb.dat$TMB<-msk.tmb.dat$snv/MSKCC.panel.size
msk.tmb.dat$POLXstatus<-ifelse(msk.tmb.dat$pole_mutations+msk.tmb.dat$pold1_mutations>=1 ,"POLE/POLD1 Mut","WT")
plot_tmb_data_polX<-msk.tmb.dat[,c("TMB","POLXstatus","GeneralTumorType")]

#get cancer types have enough values >5
ava.cancer.names<-colnames(table(plot_tmb_data_polX[,-1]))[table(plot_tmb_data_polX[,-1])[1,]>=5]
plot_tmb_data_polX<-plot_tmb_data_polX[which(plot_tmb_data_polX$GeneralTumorType%in%ava.cancer.names),]
plot_tmb_data_polX$GeneralTumorType<-as.character(plot_tmb_data_polX$GeneralTumorType)
# sort by frequency
confusionTable<-data.frame(Mut=table(plot_tmb_data_polX$GeneralTumorType,plot_tmb_data_polX$POLXstatus)[,1],Total=table(plot_tmb_data_polX$GeneralTumorType,plot_tmb_data_polX$POLXstatus)[,2])
confusionTable$fre<-confusionTable[,1]/confusionTable[,2]
confusionTable<-confusionTable[order(confusionTable$fre,decreasing = T),]

# add sample number into plot 
confusionTable$GeneralTumorType<-row.names(confusionTable)
confusionTable$newName<-with(confusionTable,paste(GeneralTumorType,"\n (",Mut,"/",Total-Mut,")",sep = ""))

plot_tmb_data_polX_merge<-merge(plot_tmb_data_polX,confusionTable,by="GeneralTumorType")
plot_tmb_data_polX_merge$GeneralTumorType<-plot_tmb_data_polX_merge$newName

plot_tmb_data_polX_merge$GeneralTumorType<-factor(plot_tmb_data_polX_merge$GeneralTumorType,levels = confusionTable$newName)

pbox1<-plot_paired_boxplot(plot_tmb_data_polX_merge)


# pbox2<-plot_paired_boxplot(plot_tmb_data_pold1,type="POLD1")

pdf("../results/Figure_1C.pdf",width = 12,height = 7)
print(pbox1)
dev.off()
