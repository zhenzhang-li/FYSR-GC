library(data.table)
library(randomForest)
library(ggplot2)
library(reshape2)
library(pROC)
library(dplyr)
library(robustbase)
library(cowplot)

s1 = "SBS3.SBS.geneCumulativeContributionAbundance.txt"
s10 = "SBS6.SBS.geneCumulativeContributionAbundance.txt"
s1 = as.data.frame(fread(s1, header=F, sep="\t"))
s1$V1 = paste(s1$V1,'_SBS18',sep='')
s10 = as.data.frame(fread(s10, header=F, sep="\t"))
s10$V1 = paste(s10$V1,'_SBS44',sep='')
s10 = s10[-1, ]
data = data.frame(rbind(s1, s10))

GroupT = as.data.frame(fread("Group.txt"))
NamesT = GroupT$SampleID
index = NULL
for(i in 1:length(NamesT))
{
  ix = which(data[1,] == NamesT[i])
  if(length(ix)>0)
  {
    index = c(index, ix)
  }    
}
dataT = data[, c(1, index)] 
dataT[1,1] = 'SampleID'
write.table(dataT, file="Test.matrix.txt",quote = F, row.names = F, col.names = F, sep="\t")


load("train1.rf.RData")

#----- predice resule of test.set by markers 
conf2 <- read.table('Group.txt',header=T,sep='\t',row.names=1,check.names=F)
colnames(conf2) <- c('state')
dat2 <- read.table('Test.matrix.txt',header=T,sep='\t',row.names=1,check.names=F)
dat2 <- data.frame(t(dat2))
set.seed(999)
p.test<-predict(train1.rf, dat2, type='prob') 
Dnamess <- intersect(rownames(as.matrix(p.test)),rownames(as.matrix(conf2)))
pre.test <- as.data.frame(cbind(pre.value=as.matrix(p.test[Dnamess,2]),as.matrix(conf2)[Dnamess,]))
colnames(pre.test) <- c('pre.value', 'state')
write.table(pre.test,'cross_validation.marker.predict.in.test.txt',sep='\t',quote=F)
conf2 <- conf2[Dnamess, ]
number_ticks <- function(n) {function(limits) pretty(limits, n)} ## for axis tick number

##### plot figure 6 #####
pre.data <- read.table('cross_validation.marker.predict.in.test.txt',head=T,check.names=F)
rownames(pre.data) <- NULL
colnames(pre.data) <- c('pre.value','state')
pre.data$state <- factor(pre.data$state,levels=c('0','1'))
colors <- c('blue','red')
pre.data <- pre.data[order(pre.data$pre.value),]
p6 <- ggplot(pre.data) +
  geom_point(aes(as.numeric(reorder(as.numeric(row.names(pre.data)),pre.value)),pre.value,color=state),size=2) +
  geom_hline(yintercept=0.5,linetype='dashed',color='grey') +
  labs(x='Samples',y='Probability of Disease',color='') +
  scale_color_manual(values= colors) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(0,1,0.5)) +
  scale_x_continuous(breaks=number_ticks(5)) +
  theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',angle=90,hjust=0.5,size=8),
        axis.title = element_text(color='black',size=8),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'),
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.text = element_text(size=8),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.key.size = unit(0.1,'in'),
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
        panel.border = element_rect(colour = 'black', fill=NA),
        panel.grid = element_blank(),
        panel.background = element_blank()
  )

#----- test.ROC -----
outcome.test = conf2
outcome.test <- sub('0','0',outcome.test)
outcome.test <- sub('1','1',outcome.test)
##### plot figure 7 #####
test.roc1 <- roc(outcome.test,p.test[,2],percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
sens.ci <- as.data.frame(ci.se(test.roc1, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c('sp','se.low','se.median','se.high'))
write.table(sens.ci,'test_ci_result.txt',sep='\t',quote=F,col.names=T,row.names=F)
test.ci <- read.table('test_ci_result.txt',head=T,check.names=F)
p7 <- ggroc(test.roc1,color='darkred',legacy.axes=TRUE)+
  annotate(geom='text',x=0.5,y=0.15,label=paste('AUC =',round(test.roc1$ci[2],4)),size=3,hjust=0) +
  annotate(geom='text',x=0.5,y=0.05,label=paste('95% CI:',round(test.roc1$ci[1],4),'-',round(test.roc1$ci[3],4)),size=3,hjust=0) +
  geom_abline(intercept = 0, slope = 1, color = 'gray') +
  geom_ribbon(data=test.ci,aes(x=1-sp,ymin=se.low,ymax=se.high),fill='red',alpha=0.2) +
  scale_x_continuous(breaks=seq(0, 1, 0.5)) +
  scale_y_continuous(breaks=seq(0, 1, 0.5)) +
  theme(axis.text.x = element_text(color='black',size=10),
        axis.text.y = element_text(color='black',size=10,angle=90,hjust=0.5),
        axis.title = element_text(color='black',size=10),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'),
        legend.position = 'none',
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
        panel.border = element_rect(colour = 'black', fill=NA),
        panel.grid = element_blank(),
        panel.background = element_blank()
  )

pdf('BGI.RUC.pdf',width=12,height=6)
print(plot_grid(p6,p7,ncol=2,align='hv'))
dev.off()
