rfcv1 <- function (trainx, trainy, cv.fold = 5, scale = "log", step = 0.5,mtry = function(p) max(1, floor(sqrt(p))), recursive = FALSE,...)
{
  classRF <- is.factor(trainy)
  n <- nrow(trainx)
  p <- ncol(trainx)
  if (scale == "log") 
  {
    k <- floor(log(p, base = 1/step))
    n.var <- round(p * step^(0:(k - 1)))
    same <- diff(n.var) == 0
    if (any(same))
      n.var <- n.var[-which(same)]
    if (!1 %in% n.var)
      n.var <- c(n.var, 1)
  }else {
    n.var <- seq(from = p, to = 1, by = step)
  }
  k <- length(n.var)
  cv.pred <- vector(k, mode = "list")
  for (i in 1:k) cv.pred[[i]] <- rep(0,length(trainy))
  if (classRF) 
  {
    f <- trainy
  }else {
    f <- factor(rep(1:5, length = length(trainy))[order(order(trainy))])
  }
  nlvl <- table(f)
  idx <- numeric(n)
  for (i in 1:length(nlvl)) 
  {
    idx[which(f == levels(f)[i])] <- sample(rep(1:cv.fold,length = nlvl[i]))
  }
  res=list()
  accuracy=list()
  gini=list()
  for (i in 1:cv.fold) 
  {
    all.rf <- randomForest(trainx[idx != i, , drop = FALSE],trainy[idx != i],importance = TRUE)
    aa = predict(all.rf,trainx[idx == i, , drop = FALSE],type="prob")
    cv.pred[[1]][idx == i] <- as.numeric(aa[,2])
    impvar <- (1:p)[order(all.rf$importance[, 3], decreasing = TRUE)]
    acc_value <- as.numeric(importance(all.rf)[, 3])
    gini_value <- as.numeric(importance(all.rf)[, 4])
    res[[i]]=impvar
    accuracy[[i]]= acc_value
    gini[[i]] = gini_value
    for (j in 2:k) 
    {
      imp.idx <- impvar[1:n.var[j]]
      sub.rf <- randomForest(trainx[idx != i, imp.idx,drop = FALSE], trainy[idx != i])
      bb <- predict(sub.rf,trainx[idx ==i,imp.idx, drop = FALSE],type="prob")
      cv.pred[[j]][idx == i] <- as.numeric(bb[,2])
      if (recursive) 
      {
        impvar <- (1:length(imp.idx))[order(sub.rf$importance[,3], decreasing = TRUE)]
        acc_value <- as.numeric(importance(sub.ref)[, 3])
        gini_value <- as.numeric(importance(sub.ref)[, 4])
      }
      NULL
    }
    NULL
  }
  if (classRF) 
  {
    error.cv <- sapply(cv.pred, function(x) mean(ifelse(x>0.5,1,0) != as.numeric(as.character(trainy))))  
    #error.cv <- sapply(cv.pred, function(x) mean(factor(ifelse(x>0.5,1,0))!=trainy))
  }else {   
    #error.cv <- sapply(cv.pred, function(x) mean((trainy - x)^2))
    error.cv <- sapply(cv.pred, function(x) mean((as.numeric(as.character(trainy)) - x)^2))
  }
  names(error.cv) <- names(cv.pred) <- n.var
  list(n.var = n.var, error.cv = error.cv, predicted = cv.pred,res=res, acc=accuracy, gini=gini)
}


library(data.table)
library(randomForest)
library(ggplot2)
library(reshape2)
library(pROC)
library(dplyr)
library(robustbase)
library(cowplot)

s1 = "SBS1.SBS.geneCumulativeContributionAbundance.txt"
s10 = "SBS10.SBS.geneCumulativeContributionAbundance.txt"
s1 = as.data.frame(fread(s1, header=F, sep="	"))
s1$V1 = paste(s1$V1,'_SBS18',sep='')
s10 = as.data.frame(fread(s10, header=F, sep="	"))
s10$V1 = paste(s10$V1,'_SBS44',sep='')
s10 = s10[-1, ]
data = data.frame(rbind(s1, s10))
da = as.data.frame(fread("Group.txt"))  

aafeq = 0.001
ssfeq = 0
cv.fold = 10
trc = 0.95
nrep = 3
path = "result"
if(trc<1)
{	  
  trcinx0 = trcinx1 = NULL  
  inx0 = which(da$Group == 0)
  inx1 = which(da$Group == 1)
  trcinx0 = sample(1:length(inx0), floor(trc*length(inx0))) 
  trcinx1 = sample(1:length(inx1), floor(trc*length(inx1)))  
  inx0tr = inx0[trcinx0]
  inx1tr = inx1[trcinx1]
  gr1 = da[c(inx0tr,inx1tr),]
}else{
  gr1 = da
}

tmpgr1 = NULL
for(i in 1:nrep)
{
  tmpgr1 = rbind(tmpgr1, gr1)
}
write.table(tmpgr1, file=paste(path, "/Group1.txt", sep = ""), quote = F, sep="	", col.names = T,row.names = F)

Group = as.data.frame(fread(paste(path, "/Group1.txt", sep = "")))
Names = Group$SampleID
da1 = t(as.data.frame(fread(paste(path, "/Group1.txt", sep = ""),header = F, sep="	")))
index = NULL
for(i in 1:length(Names))
{
  ix = which(data[1,] == Names[i])
  if(length(ix)>0)
  {
    index = c(index, ix)
  }    
}
dataa = data[, c(1, index)] 
dataa[1,1] = 'SampleID'
write.table(dataa, file=paste(path, "/matrix.txt", sep = ""),quote = F, row.names = F, col.names = F, sep="	")

GroupT = as.data.frame(fread(paste(path, "/Group2.txt", sep = "")))
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
write.table(dataT, file=paste(path, "/Test.matrix.txt", sep = ""),quote = F, row.names = F, col.names = F, sep="	")

GroupT = as.data.frame(fread(paste(path, "/Group.txt", sep = "")))
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
write.table(dataT, file=paste(path, "/TRain.original.matrix.txt", sep = ""),quote = F, row.names = F, col.names = F, sep="	")

da2 = as.data.frame(fread(paste(path, "/matrix.txt", sep = ""),header = F, sep="	"))
Names = as.character(da2[1, ])
index = NULL
for(i in 1:length(Names))
{
  ix = which(da1[1,] == Names[i])
  if(length(ix)>0)
  {
    index = c(index, ix[1])
  }    
}
da1 = da1[, index]
data = data.frame(rbind(da1, da2))
data[2,1] = 'Type'
data = data[-3, ]
tmp = data[1:2, ]
tmp1 = data[-c(1:2), ]
meanT = apply(tmp1[,2:ncol(tmp1)], 1, function(x){x=as.numeric(as.character(x)); return(mean(x))})
ix = meanT >= aafeq
tmp1 = tmp1[ix, ]
Num = apply(tmp1[,2:ncol(tmp1)], 1, function(x){x=as.numeric(as.character(x)); return(sum(x>0))})
ix = Num/(ncol(tmp1)-1) >= ssfeq
tmp1 = tmp1[ix, ]
data = rbind(tmp,tmp1)
write.table(data, file=paste(path, "/metaphlan2.RF.txt", sep = ""),quote = F, row.names = F, col.names = F, sep="	")

infile = paste(path, "/metaphlan2.RF.txt", sep = "")
set.seed(999)
dat1 <- read.table(infile, head=T,sep="	",row.names=1,check.names=FALSE,  stringsAsFactors = FALSE)
conf<-dat1[1,]
outcome = as.factor(as.character(conf))
outcome <-as.factor(outcome)  
X <- as.data.frame(t(dat1[-1,]))
X$outcome = outcome  
result <- replicate(5, rfcv1(X[,-ncol(X)], X$outcome, cv.fold=cv.fold, step=0.9), simplify=FALSE)
save(result, file = paste(path,'/result.RData',sep=""))
error.cv <- sapply(result, '[[', 'error.cv')
error.cv.cbm <- as.data.frame(cbind(rowMeans(error.cv), error.cv))
cutoff <- min(error.cv.cbm[,1])+sd(error.cv.cbm[,1])
cutoff.num <- nrow(error.cv.cbm[error.cv.cbm[,1]<=cutoff,])
optimal.set.feature.num <- as.numeric(rownames(error.cv.cbm[error.cv.cbm[,1]<=cutoff,])[cutoff.num]) # obtain marker number    
error.cv.data <- cbind(error.cv,rowMeans(error.cv))
colnames(error.cv.data) <- c('E1','E2','E3','E4','E5','Mean') 
error.cv.data <- melt(error.cv.data)  
number_ticks <- function(n) {function(limits) pretty(limits, n)} ## for axis tick number 

#----- pick out markers from the crossvalidation result -----
k = 1
nl = length(result)
nl1 = length(result[[1]]$res)
nl2 = length(result[[1]]$res[[1]])
b <- matrix(0, ncol=nl2, nrow=nl*nl1)  ## ncol is gene, genus or mlg number, nrow = 5*10
accuracy <- matrix(0,ncol=nl2,nrow=nl*nl1) ## for plot feature importance
gini <- matrix(0,ncol=nl2,nrow=nl*nl1)

for(i in 1:nl)
{
  for(j in 1:nl1)
  {
    b[k,] <- result[[i]]$res[[j]]
    accuracy[k,] <- result[[i]]$acc[[j]]
    gini[k,] <- result[[i]]$gini[[j]]
    k = k+1
  }
}
mlg.list <- b[,1:10]
list <- c()
k = 1
for(i in 1:10)
{
  for(j in 1:50)
  {
    list[k] <- mlg.list[j,i]
    k = k+1
  }
}  
mlg.sort <- as.matrix(table(list))
mlg.sort <- mlg.sort[rev(order(mlg.sort[,1])),]
pick <- as.numeric(names(head(mlg.sort,optimal.set.feature.num)))

########################
# update by lw
optimal.set.feature.num <- length(pick)
p1 <- ggplot(error.cv.data,aes(x=Var1,y=value,color=Var2)) +
  geom_line() +
  geom_vline(xintercept = optimal.set.feature.num, color='red') +
  coord_trans(x="log2") +
  scale_x_continuous(breaks=c(1,2,5,10,20,50),labels=c(c(1,2,5,10,20,50))) +
  labs(x='Number of variables',y='CV Error') +
  scale_color_manual(values= c('gray','gray','gray','gray','gray','black')) +
  theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',size=8,angle=90,hjust=0.5),
        axis.title = element_text(color='black',size=8),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'),
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
        legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill=NA),
        panel.grid = element_blank(),
        panel.background = element_blank()
  )
########################

tmp = X[,-ncol(X)]
mlg.pick <- colnames(tmp)[pick]
write.table(mlg.pick,paste(path,'/cross_validation_pick.txt',sep=""),sep='	',quote=F) 
colnames(accuracy) <- rownames(dat1)[-1]
colnames(gini) <- rownames(dat1)[-1]  
median_acc <- as.matrix(rowMedians(t(accuracy)))
median_gini <- as.matrix(rowMedians(t(gini)))  
importance.data <- cbind(median_acc,median_gini)
colnames(importance.data) <- c('MeanDecreaseAccuracy','MeanDecreaseGini')
write.table(importance.data,paste(path,'/feature.importance.xls',sep=""),sep='	',quote=F,col.names=NA)

#----- predict result of train.set by markers -----
train1 <- X[,c(pick,nrow(dat1))] ##
set.seed(999)
train1 <- data.frame(train1)
train1.rf <- randomForest(outcome~., data = train1,importance = TRUE)
save(train1.rf, file = paste(path,'/train1.rf.RData',sep=""))
train1.pre <- predict(train1.rf,type='prob')
p.train <- train1.pre[,2]
rowNNs = rownames(as.matrix(p.train))
rowNNs = sapply(strsplit(rowNNs, "\\."),"[[", 1)
urowNNs = unique(rowNNs)
Tnamess <- intersect(urowNNs, unique(Group$SampleID))
combine <- matrix(0, length(Tnamess), 3)
for(i in 1:length(Tnamess))
{
  combine[i, 2] = Tnamess[i]
  ix = rowNNs == Tnamess[i]
  combine[i, 1] = mean(as.numeric(as.character(p.train[ix])))
  ix = Group$SampleID == Tnamess[i]
  combine[i, 3] = Group$Group[ix][1]
}
colnames(combine) <- c('predict.value','SampleID','Group')
write.table(combine,paste(path,'/cross_validation.marker.predict.in.train.txt',sep=""),sep='	',quote=F)
###### plot figure 2 #####
temp.data <- read.table(paste(path,'/cross_validation.marker.predict.in.train.txt',sep=""),head=T,check.names=F)
temp.data$Group <- factor(temp.data$Group,levels=c('0','1'))
colors <- c('blue','red')

p2 <- ggplot(data=temp.data,aes(x=Group,y=predict.value)) +
  stat_boxplot(geom = 'errorbar',width=0.3,color=colors) +
  geom_boxplot(aes(color=Group),lwd=0.5,outlier.shape=1,notch=F) + 
  geom_hline(yintercept=0.5,linetype='dashed',color='grey') +
  labs(x='',y='Probability of Disease',fill= '') +
  scale_color_manual(values= colors) +
  scale_y_continuous(limits = c(0, 1),breaks=seq(0, 1, 0.5)) + 
  theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',size=8,angle=90,hjust=0.5),
        axis.title = element_text(color='black',size=8),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'),
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
        legend.position = 'none',
        panel.border = element_rect(colour = 'black', fill=NA),
        panel.grid = element_blank(),
        panel.background = element_blank()
  )
#----- ROC in train -----
##### plot figure 3 #####
train.roc <- roc(outcome,p.train,percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,direction = "<")
sens.ci <- as.data.frame(ci.se(train.roc, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c('sp','se.low','se.median','se.high'))
write.table(sens.ci,paste(path,'/train_ci_result.txt',sep=""),sep='	',quote=F,col.names=T,row.names=F)
train.ci <- read.table(paste(path,'/train_ci_result.txt',sep=""),head=T,check.names=F)

p3 <- ggroc(train.roc,color='darkred',legacy.axes=TRUE)+
  annotate(geom='text',x=0.5,y=0.15,label=paste('AUC =',round(train.roc$ci[2],4)),size=3,hjust=0) +
  annotate(geom='text',x=0.5,y=0.05,label=paste('95% CI:',round(train.roc$ci[1],4),'-',round(train.roc$ci[3],4)),size=3,hjust=0) +
  geom_abline(intercept = 0, slope = 1, color = 'gray') +
  geom_ribbon(data=train.ci,aes(x=1-sp,ymin=se.low,ymax=se.high,fill='red'),fill='red',alpha=0.2) +
  scale_x_continuous(breaks=seq(0, 1, 0.5)) +
  scale_y_continuous(breaks=seq(0, 1, 0.5)) +
  theme(axis.text.x = element_text(color='black',size=8),
        axis.text.y = element_text(color='black',size=8,angle=90,hjust=0.5),
        axis.title = element_text(color='black',size=8),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'),
        legend.position = 'none',
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
        panel.border = element_rect(colour = 'black', fill=NA),
        panel.grid = element_blank(),
        panel.background = element_blank()
  )

###########################  update by lw
#----- predice resule of train.original.set by markers ----- 
conf2 <- read.table(paste(path,'/Group.txt',sep=""),header=T,sep='	',row.names=1,check.names=F)
colnames(conf2) <- c('state')
dat2 <- read.table(paste(path,'/TRain.original.matrix.txt',sep=""),header=T,sep='	',row.names=1,check.names=F)
dat2 <- data.frame(t(dat2))
set.seed(999)
p.test<-predict(train1.rf, dat2, type='prob') 
Dnamess <- intersect(rownames(as.matrix(p.test)),rownames(as.matrix(conf2)))
pre.test <- as.data.frame(cbind(pre.value=as.matrix(p.test[Dnamess,2]),as.matrix(conf2)[Dnamess,]))
colnames(pre.test) <- c('pre.value', 'state')
write.table(pre.test,paste(path,'/cross_validation.marker.predict.in.train.original.txt',sep=""),sep='	',quote=F)
conf2 <- conf2[Dnamess, ]

##### plot figure 4 #####
train.original.data <- read.table(paste(path,'/cross_validation.marker.predict.in.train.original.txt',sep=""),head=T,check.names=F)
rownames(train.original.data) <- NULL
colnames(train.original.data) <- c('pre.value','state')
train.original.data$state <- factor(train.original.data$state,levels=c('0','1'))
colors <- c('blue','red')
train.original.data <- train.original.data[order(train.original.data$pre.value),]
p4 <- ggplot(train.original.data) +
  geom_point(aes(as.numeric(reorder(as.numeric(row.names(train.original.data)),pre.value)),pre.value,color=state),size=2) +
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
        #  	legend.key.width = unit(0.2, 'in'),
        #		legend.key.height = unit(0.2, 'in'),
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
##### plot figure 5 #####
test.roc <- roc(outcome.test,p.test[,2],percent=FALSE,partial.auc.correct=TRUE,ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,plot=F,auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE)
sens.ci <- as.data.frame(ci.se(test.roc, specificities=seq(0, 1, 0.05)))
sens.ci <- setNames(cbind(rownames(sens.ci),sens.ci,row.names=NULL),c('sp','se.low','se.median','se.high'))
write.table(sens.ci,paste(path,'/train.original_ci_result.txt',sep=""),sep='	',quote=F,col.names=T,row.names=F)
test.ci <- read.table(paste(path,'/train.original_ci_result.txt',sep=""),head=T,check.names=F)
p5 <- ggroc(test.roc,color='darkred',legacy.axes=TRUE)+
  annotate(geom='text',x=0.5,y=0.15,label=paste('AUC =',round(test.roc$ci[2],4)),size=3,hjust=0) +
  annotate(geom='text',x=0.5,y=0.05,label=paste('95% CI:',round(test.roc$ci[1],4),'-',round(test.roc$ci[3],4)),size=3,hjust=0) +
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
###############################	


#----- predice resule of test.set by markers 
conf2 <- read.table(paste(path,'/Group2.txt',sep=""),header=T,sep='	',row.names=1,check.names=F)
colnames(conf2) <- c('state')
dat2 <- read.table(paste(path,'/Test.matrix.txt',sep=""),header=T,sep='	',row.names=1,check.names=F)
dat2 <- data.frame(t(dat2))
set.seed(999)
p.test<-predict(train1.rf, dat2, type='prob') 
Dnamess <- intersect(rownames(as.matrix(p.test)),rownames(as.matrix(conf2)))
pre.test <- as.data.frame(cbind(pre.value=as.matrix(p.test[Dnamess,2]),as.matrix(conf2)[Dnamess,]))
colnames(pre.test) <- c('pre.value', 'state')
write.table(pre.test,paste(path,'/cross_validation.marker.predict.in.test.txt',sep=""),sep='	',quote=F)
conf2 <- conf2[Dnamess, ]

##### plot figure 6 #####
pre.data <- read.table(paste(path,'/cross_validation.marker.predict.in.test.txt',sep=""),head=T,check.names=F)
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
        #  	legend.key.width = unit(0.2, 'in'),
        #		legend.key.height = unit(0.2, 'in'),
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
write.table(sens.ci,paste(path,'/test_ci_result.txt',sep=""),sep='	',quote=F,col.names=T,row.names=F)
test.ci <- read.table(paste(path,'/test_ci_result.txt',sep=""),head=T,check.names=F)
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

pdf(paste(path,'/BGI.RUC.pdf',sep=""),width=25,height=6)
print(plot_grid(p1,p2,p3,p4,p5,p6,p7,ncol=7,align='hv'))
dev.off()

#----- plot importance of top n features -----
kl = length(pick)
imp.data <- read.table(paste(path,'/feature.importance.xls',sep=""),header=T,check.names=F)
imp.data <- setNames(cbind(rownames(imp.data),imp.data,row.names=NULL),c('Species','MeanDecreaseAccuracy','MeanDecreaseGini'))
imp.data <- arrange(imp.data,desc(MeanDecreaseAccuracy))
imp.data <- head(imp.data,n=kl)
number_ticks <- function(n) {function(limits) pretty(limits, n)} ## for axis tick number

p8 <- ggplot(imp.data,aes(x=reorder(Species,MeanDecreaseAccuracy),y=MeanDecreaseAccuracy)) +
  geom_bar(stat = 'identity',width=0.7) + 
  labs(y='MDA',x='',fill='',color='') +
  scale_y_continuous(breaks=number_ticks(3)) +
  coord_flip() +
  theme(axis.text.x=element_text(color='black',size=10),
        axis.text.y=element_text(color='black',size=9),
        axis.title.x=element_text(color='black',size=10),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'), 
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
        panel.border = element_rect(colour = 'black', fill=NA),
        panel.grid = element_blank(),
        panel.background = element_blank()
  )

imp.data <- arrange(imp.data,desc(MeanDecreaseGini))
imp.data <- head(imp.data,n=kl)
p9 <- ggplot(imp.data,aes(x=reorder(Species,MeanDecreaseGini),y=MeanDecreaseGini)) +
  geom_bar(stat = 'identity',width=0.7) + 
  labs(y='MDG',x='',fill='',color='') +
  scale_y_continuous(breaks=number_ticks(3)) +
  coord_flip() +
  theme(axis.text.x=element_text(color='black',size=10),
        axis.text.y=element_text(color='black',size=9),
        axis.title.x=element_text(color='black',size=10),
        axis.line = element_line(color='black'),
        axis.ticks = element_line(color='black'), 
        plot.margin = unit(c(0.351, 0.05, 0.05, 0.051), 'in'),
        panel.border = element_rect(colour = 'black', fill=NA),
        panel.grid = element_blank(),
        panel.background = element_blank()
  )

pdf(paste(path,'/feature.importance.pdf',sep=""),width=16,height=10)
print(plot_grid(p8,p9,ncol=2,align='hv'))
dev.off()  
dataTXT = t(c(path,train.roc$ci[2],test.roc$ci[2],test.roc1$ci[2]))
write.table(dataTXT, file=paste(path,"/Final.AUC.txt",sep=""),quote=F,col.names=F,row.names=F,sep="	")
