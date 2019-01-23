suppressMessages(library(ggsci))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(tidyr))
suppressMessages( library(survminer))
SurvivalCurvePlot <- function(fit, PSD, pvalue, genename,panel=c("#FED439FF","#709AE1FF"),label=c()){
  md<-surv_median(fit)

  if(length(label)>=2){
    lb<-paste(label,"(",md$median,"month)",sep="")
    BestCurve <- ggsurvplot(fit, data = PSD, 
                            pval = pvalue,
                            xlab = "Time(Months)",
                            legend.labs=lb,
                            conf.int = FALSE,
                            legend = c(0.7, 0.88),
                            palette =panel,
                            surv.median.line="hv",
                            risk.table = TRUE, 
                            risk.table.fontsize = 6,
                            risk.table.y.text.col = TRUE, 
                            risk.table.height = 0.30, 
                            legend.title = paste("Group(Median)"),
                            title=genename,
                            pval.size = 7,
                            font.main = c(20, "bold", "black"),
                            font.x = c(18, "plain", "black"), # feature of curve xlab title
                            font.y = c(18, "plain", "black"), # feature of curve ylab title
                            font.tickslab = c(18, "plain", "black"))
  }else{
    BestCurve <- ggsurvplot(fit, data = PSD, 
                            pval = pvalue,
                            xlab = "Time(Months)",
                            # legend.labs=lb,
                            conf.int = FALSE,
                            legend = c(0.7, 0.88),
                            palette =panel,
                            surv.median.line="hv",
                            risk.table = TRUE, 
                            risk.table.fontsize = 6,
                            risk.table.y.text.col = TRUE, 
                            risk.table.height = 0.30, 
                            legend.title = paste("Group(Median)"),
                            title=genename,
                            pval.size = 7,
                            font.main = c(20, "bold", "black"),
                            font.x = c(18, "plain", "black"), # feature of curve xlab title
                            font.y = c(18, "plain", "black"), # feature of curve ylab title
                            font.tickslab = c(18, "plain", "black"))
  }
  # 

  BestCurve$plot <- BestCurve$plot+ theme(legend.text = element_text(size = 16, color = "black", face = "plain"),
                                          axis.title = element_text(color = "black", hjust=0.5, vjust=0.8, size=20, face="plain"),
                                          legend.title = element_text(angle=0, hjust=0.8, vjust=0.8, size=18, face="plain"))
  BestCurve$table<- BestCurve$table + theme(axis.title.y = element_blank(),
                                            axis.text.x = element_text(color = "black", hjust=0.5, vjust=0.8, size=18, face="plain"),
                                            
                                            )
                                            
                                            
  return(BestCurve)
}

batchUnivarCOXfun<-function(surv,data){
  df<-data.frame(name=c("coef","Hazard Ratio","CI (lower)","CI (upper)","se(coef)","z","Pr(>|Z|)","C-index","C-index-se","testCph.p"))
  pro=ncol(data)/10000
  j=1;
  k=1;
  for (i in 1:ncol(data)) {
    fml=as.formula(paste("surv ~ ",names(data)[i]))
    fit<-coxph(surv ~ data[,i])
    p<-data.frame(summary(fit)$coefficients)[1,5]
    test=cox.zph(coxph(surv ~ data[,i]))
    testp= test$table[3]
    #     x<-c(x,y)
    #     coefficients=summary(a)$coefficients
    y=prettify(summary(fit))
    y=y[,c(-1,-8,-9)]
    x=c(t(y)[,1],p,summary(fit)$concordance,testp)
    df=data.frame(df,x)
    if(i==pro*k){
      print(paste(j*k,"items finished"))
      k=k+1
    }
  }
  row.names(df)<-df[,1]
  df<-t(df[,-1])
  row.names(df)=names(data)
  return(df)
}

bar_plot_function<-function(df){
  
  df$tumor_type<-paste(df$tumor_type,"(",df$allMutated,"/",df$Total,")",sep="")
  # plot frame 
  plot.frame<-df[,c(1:4,7)]
  plot.frame.melt<-reshape2::melt(plot.frame,id.vars=c("tumor_type","totalFre"))
  names(plot.frame.melt)=c("Cancer","totalFre","Mutation","Frequency")
  
  plot.frame.melt$Cancer<-factor(plot.frame.melt$Cancer,levels = df$tumor_type)
  plot.frame.melt$Mutation<-factor(plot.frame.melt$Mutation,levels = c("POLE.POLD1","POLE","POLD1"))
  # adjust text position 
  text.data <- ddply(plot.frame.melt, .(Cancer),
                       transform, pos = cumsum(Frequency) - (0.5 * Frequency))
  # text.data$tag<-plot.frame.melt$Frequency>1
  # 
  # ggplot2
 p<- ggplot(plot.frame.melt,aes(x = Cancer, y = Frequency,fill = Mutation)) + 
    geom_bar(position = "stack",stat = "identity")+ylim(0,24)+theme(legend.position=c(0.5,0.8),legend.box = "horizontal")
 
 # geom_text(
 #   aes(label = n, y = n - 1), # << move each label down by 1 unit
 #   position = position_stack(), 
 #   color = "white", fontface = "bold", size = 8
 # )  
 # 
     p<-p+geom_text(data=text.data,aes(x=Cancer,y=pos,label=Frequency),size=3,colour="white")
    p<-p+annotate("text", x = df$tumor_type, y = df$totalFre+0.1, label = paste(df$totalFre,"%",sep=""),angle = 90, hjust = -.05, size =4)
  p<-p  + theme(plot.subtitle = element_text(vjust = 1), 
                  plot.caption = element_text(vjust = 1), 
                  axis.line = element_line(linetype = "dashed"), 
                  axis.ticks = element_line(size = 0.7), 
                  axis.text = element_text(face = "bold", hjust = 0.95), 
                  axis.text.x = element_text(vjust = 0.2, angle = 90,size =10), 
                  panel.background = element_rect(fill = NA))
  
  p<-p+scale_fill_jama()
  p
  
  

  
  return(p)
}



# plot paired box plot 

plot_paired_boxplot<-function(df){
  
  colnames(df)[3]<-"Status"
  pbox<-ggpubr::ggboxplot(df, x = "GeneralTumorType", y = "TMB",
                          palette = "lancet",fill="Status",width=0.5)
  pbox<-pbox +ggpubr::stat_compare_means(aes(group = Status), label = "p.signif")
  pbox<-pbox+theme(plot.subtitle = element_text(vjust = 1), 
             plot.caption = element_text(vjust = 1), 
             axis.line = element_line(linetype = "solid"), 
             axis.ticks = element_line(size = 0.7), 
             axis.text = element_text(face = "bold", hjust = 0.95), 
             axis.text.x = element_text(vjust = 0.2, angle = 90,size =10), 
             panel.background = element_rect(fill = NA))
  return(pbox)
  
  
}
