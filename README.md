---
title: "Phytohormone profiling in an evolutionary framework"
author: Vojtěch Schmidt, Roman Skokan, Katarina Kurtović, Stanislav Vosolsobě, Roberta
  Filepová, Samuel Haluška, Thomas Depaepe, Anthony Pil, Petre Dobrev, Václav Motyka,
  Dominique Van Der Straeten, Jan Petrášek
date: "`r Sys.Date()`"
output:
  html_document:
    df_print: paged
  pdf_document: default
---

# Phytohormone profiling in an evolutionary framework

### Vojtěch Schmidt, Roman Skokan, Katarina Kurtović, Stanislav Vosolsobě, Roberta Filepová, Samuel Haluška, Thomas Depaepe, Anthony Pil, Petre Dobrev, Václav Motyka, Dominique Van Der Straeten, Jan Petrášek
    

# Statistical analysis
  
## Full R source code written by Stanislav Vosolsobě
  

### Required libraries

```{r}
library(emmeans)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(plotfunctions)
```



## ETHYLENE PRODUCTION
  
```{r, results='hide'}

eth <- read.table("ethylene",header = T)
sample <- paste(eth$Strain,eth$Type,sep = " ")

pdf("eth_dist.pdf",width=6,height=6)
par(las=2,mar=c(16,4,2,2))
boxplot(eth$Ethylene~sample,xlab="",ylab="Ethylene production",col=c("gray","seagreen"))
dev.off()
pdf("eth_dist_trans.pdf",width=6,height=6)
par(las=2,mar=c(16,4,2,2))
boxplot(log(eth$Ethylene)~sample,xlab="",ylab="Ethylene production",col=c("gray","seagreen"))
dev.off()
```

Ethylene production data are non-normally distributed, but can be normalised by log transformation and simple ANOVA can be used for statistical analysis.

### ANOVA

```{r}
m1 <- aov(log(eth$Ethylene)~eth$Strain*eth$Type)
summary(m1)
```

### Confidence intervals

```{r}
par(las=1,mar=c(6,20,2,2))
em1 <- emmeans(m1, pairwise ~ Type | Strain)
summary(em1)
plot(em1,comparisons=T,xlab="Ethylene production (log scale)",ylab="Sample type")
pval <- summary(em1)$contrasts$p.value

pdf("eth_dist_stat.pdf",width=7,height=7)
par(las=2,mar=c(16,4,2,2))
bp <- boxplot(eth$Ethylene~sample,xlab="",ylab="Ethylene production",plot=F)
boxplot(eth$Ethylene~sample,xlab="",ylab="Ethylene production",ylim=c(0,65+max(c(bp$stats[5,],bp$out))),xaxt="n",col=alpha(c("gray","seagreen"),f=0.7),border=c("gray","seagreen"))
axis(side=1,at=1:(2*length(pval)),labels=lapply(strsplit(bp$names,split=" "),`[`,2))
axis(side=1,at=(1:length(pval))*2-0.5,labels=lapply(strsplit(bp$names,split=" "),`[`,1)[(1:length(pval))*2],line = 3,lwd = 0)
for(i in 1:length(pval)){
  lab="n.s."
  if(pval[i]<0.05) lab = "*"
  if(pval[i]<0.01) lab = "**"
  if(pval[i]<0.001) lab = "***"
  y0 <- bp$stats[5,2*i]
  text(x = 2*i-0.5,y = y0+60,labels = lab,adj = 0.5)
  text(x = 2*i-0.5,y = y0+40,labels = formatC(pval[i],digits = 3),adj = 0.5,cex=0.75)
  segments(2*i-1,y0+25,2*i,y0+25,col="grey")
  segments(2*i-1,bp$stats[5,2*i-1]+10,2*i-1,y0+25,col="grey")
  segments(2*i,bp$stats[5,2*i]+10,2*i,y0+25,col="grey")
}
dev.off()
```

## PHYTOHORMONE CONTENT

### Estimation of means, sd and relative errors for all samples

```{r}
fh <- read.table("hormones",header = T)
sample <- as.factor(paste(fh$Strain,fh$Type,fh$Stage,sep = " "))



rerS <- matrix(NA,ncol=ncol(fh)-6,nrow=nlevels(sample))
rerT <- matrix(NA,ncol=ncol(fh)-6,nrow=nlevels(sample))
colnames(rerS) <- colnames(fh)[7:ncol(fh)]
rownames(rerS) <- unique(sample)
colnames(rerT) <- colnames(fh)[7:ncol(fh)]
rownames(rerT) <- unique(sample)

meanS <- matrix(NA,ncol=ncol(fh)-6,nrow=nlevels(sample))
sdS <- matrix(NA,ncol=ncol(fh)-6,nrow=nlevels(sample))
colnames(meanS) <- colnames(fh)[7:ncol(fh)]
rownames(meanS) <- unique(sample)
colnames(sdS) <- colnames(fh)[7:ncol(fh)]
rownames(sdS) <- unique(sample)
total <- matrix(NA,ncol=4*(ncol(fh)-6),nrow=nlevels(sample))

labs <- matrix("",ncol=4*(ncol(fh)-6),2)
labs[1,1+4*(0:(ncol(fh)-7))]<-colnames(fh)[-c(1:6)]
labs[2,]<-rep(c("Mean","Sd","Relative error (biological)","Relative error (technical)"),ncol(fh)-6)

for(i in 1:nlevels(sample)){
  print(i)
  sfh <- fh[sample==unique(sample)[i],]
  print(unique(sample)[i])
  for(h in 1:(ncol(fh)-6)){
    if(!is.na(mean(sfh[,6+h],na.rm=T))){
      #print("OK")
      if(nlevels(as.factor(sfh$Rep[!is.na(sfh[,6+h])]))>1){
        lmt <-aov(sfh[,6+h]~as.factor(sfh$Rep))
        cf <- coef(lmt)
        tm <- c(cf[1],cf[2:length(cf)]+cf[1])
        rerS[i,h] <- sd(tm)/mean(tm)
        rerT[i,h] <- sqrt(unlist(summary(lmt))[6])/mean(tm)
      }else{
        rerS[i,h] <- sd(sfh[,6+h],na.rm=T)/mean(sfh[,6+h],na.rm=T)
      }
      meanS[i,h] <- mean(sfh[,6+h],na.rm=T)
      sdS[i,h] <- sd(sfh[,6+h],na.rm=T)
      total[i,4*(h-1)+c(1,2,3,4)] <- c(meanS[i,h], sdS[i,h], rerS[i,h], rerT[i,h])
    }
  }
}

pdf("fh_sample_heat.pdf",width=12,height=28)
par(las=2,mar=c(2,18,7,2))
image(z=t(rerS)[,nrow(rerS):1],xaxt="n",yaxt="n")
axis(side = 2,at=(1:nlevels(sample)-1)/(nlevels(sample)-1),labels = rev(unique(sample)))
axis(side = 3,at=(1:(ncol(fh)-6)-1)/((ncol(fh)-6)-1),labels = colnames(fh)[7:ncol(fh)])
box()
dev.off()

pdf("fh_tech_heat.pdf",width=12,height=28)
par(las=2,mar=c(2,18,7,2))
image(z=t(rerT)[,nrow(rerT):1],xaxt="n",yaxt="n")
axis(side = 2,at=(1:nlevels(sample)-1)/(nlevels(sample)-1),labels = rev(unique(sample)))
axis(side = 3,at=(1:(ncol(fh)-6)-1)/((ncol(fh)-6)-1),labels = colnames(fh)[7:ncol(fh)])
box()
dev.off()

pdf("legend.pdf",width = 1.5,height = 2)
par(mar=c(0.5,3,0.5,0.1),las=1)
lgd <- matrix(0:200/100,nrow = 1)
image(lgd,,xaxt="n",yaxt="n")
axis(side = 2,at=0:4/10*5/2,labels = format(0:4/10*5,nsmall = 1))
box()
dev.off()


pdf("fh_sample_boxplot.pdf",width=10,height=7)
par(las=2,mar=c(6,4,2,2))
boxplot(rerS,ylab="Relative error of mean")
dev.off()

pdf("fh_tech_boxplot.pdf",width=10,height=7)
par(las=2,mar=c(6,4,2,2))
boxplot((rerT),ylab="Relative error of mean")
dev.off()

write.csv(rerS,"sample.csv")
write.csv(rerT,"techni.csv")
write.csv(meanS,"means.csv")
write.csv(sdS,"sds.csv")
write.csv(total,"total.csv")
write.csv(labs,"labs.csv")
```


### AXENIC vs. NON-AXENIC COMPARISON BIOMASS

```{r}
fh <- read.table("hormones",header = T)
stage <- subset(fh,subset = (To_axenic.vs.non.axenic=="Y")&(Type=="biomass"))
vars <- (stage$Strain)
spec <- unlist(lapply(strsplit(vars,split = "_"),`[`,1))
axen <- unlist(lapply(strsplit(unlist(lapply(strsplit(vars,split = "_"),`[`,2)),split="-"),`[`,1))
axen[axen!="CCAC"] <- "non-axenic"
axen[axen=="CCAC"] <- "axenic"
stage$Strain <- spec
stage$Stage <- axen

vars <- as.factor(stage$Strain)
pstage <- matrix(NA,ncol=ncol(stage)-6,nrow=nlevels(vars)+2)

colnames(pstage) <- colnames(stage)[7:ncol(stage)]
rownames(pstage) <- c(unique(vars),as.factor(c("Axenicity:Strain","Axenicity")))

for(i in 7:ncol(stage)){
  print(colnames(stage)[i])
  actual <- subset(stage,subset = !is.na(stage[,i]))
  tech <- paste(actual$Strain,actual$Stage,actual$Rep,sep = " ")
  if(length(unique(actual$Stage))<2) next
  if(length(unique(actual$Strain))<2){
    lme1 <- lmer(log(actual[,i])~actual$Stage + (1|tech))
    lme3 <- lmer(log(actual[,i])~(1|tech))
    em1 <-tryCatch({
      emmeans(lme1, ~ Stage)
    },error=function(cond){
      message(cond)
      return=NA
    })
    if(is.na(em1)) next
    pp <- summary(pairs(em1))
    pstage[rownames(pstage)==unique(actual$Strain),i-6] <- log(pp$p.value) - log(0.05)
    pstage[nlevels(vars)+1,i-6] <- log(anova(lme1,lme3)$`Pr(>Chisq)`[2])- log(0.05)
  }else{
    lme1 <- lmer(log(actual[,i])~actual$Strain * actual$Stage + (1|tech))
    lme2 <- lmer(log(actual[,i])~actual$Strain + actual$Stage + (1|tech))
    lme3 <- lmer(log(actual[,i])~actual$Strain  + (1|tech))
    em1 <-tryCatch({
      emmeans(lme1, ~ Stage | Strain)
    },error=function(cond){
      message(cond)
      return=NA
    })
    if(is.na(em1)) next
    pp <- summary(pairs(em1))
    pstage[,i-6] <- log(pp$p.value[match(rownames(pstage), pp$Strain)]) - log(0.05)
    pstage[nlevels(vars)+2,i-6] <- log(anova(lme1,lme2)$`Pr(>Chisq)`[2])- log(0.05)
    pstage[nlevels(vars)+1,i-6] <- log(anova(lme2,lme3)$`Pr(>Chisq)`[2])- log(0.05)
  }
}  

# Color scale
rc1 <- colorRampPalette(colors = c("gray", "gray"), space = "rgb")(5*abs(max(na.rm = T,pstage)))
rc2 <- colorRampPalette(colors = c("papayawhip","darkred","black"), space = "rgb")(5*abs(min(na.rm = T,pstage)))
rampcols_blank <- c(rc1, rc2)

# Heatmap
pdf("axen_biomass_stats.pdf",width=12,height=3)
par(las=2,mar=c(2,12,7,2))
image(z=t(-1*pstage)[,nrow(pstage):1],xaxt="n",yaxt="n",col=rampcols_blank)
axis(side = 2,at=(1:nrow(pstage)-1)/(nrow(pstage)-1),labels = rev(rownames(pstage)))
axis(side = 3,at=(1:(ncol(stage)-6)-1)/((ncol(stage)-6)-1),labels = colnames(stage)[7:ncol(stage)])
box()
dev.off() 

# Legend
pdf("axen_biomass_legend.pdf",width = 1.5,height = 2)
par(mar=c(0.5,4,0.5,0.1),las=1)
rs <- range(-1*pstage,na.rm=T)
lgd <- matrix((10*rs[1]):(10*rs[2])/10,nrow = 1)
image(lgd,xaxt="n",yaxt="n",col=rampcols_blank)
cap <- c(0.05,0.001,1e-6,1e-12,1e-18)
lsg <- (-log(cap)+log(cap[1])-rs[1])/(rs[2]-rs[1])
axis(side = 2,at=lsg,labels = prettyNum(cap,nsmall = 2))
box()
dev.off()

# Table
write.csv(exp(pstage+log(0.05)),file = "axen_biomass.csv")
```

### AXENIC vs. NON-AXENIC COMPARISON MEDIUM

```{r}
fh <- read.table("hormones",header = T)
stage <- subset(fh,subset = (To_axenic.vs.non.axenic=="Y")&(Type=="medium"))
vars <- (stage$Strain)
spec <- unlist(lapply(strsplit(vars,split = "_"),`[`,1))
axen <- unlist(lapply(strsplit(unlist(lapply(strsplit(vars,split = "_"),`[`,2)),split="-"),`[`,1))
axen[axen!="CCAC"] <- "non-axenic"
axen[axen=="CCAC"] <- "axenic"
stage$Strain <- spec
stage$Stage <- axen

vars <- as.factor(stage$Strain)
pstage <- matrix(NA,ncol=ncol(stage)-6,nrow=nlevels(vars)+2)

colnames(pstage) <- colnames(stage)[7:ncol(stage)]
rownames(pstage) <- c(unique(vars),as.factor(c("Axenicity:Strain","Axenicity")))

for(i in 7:ncol(stage)){
  print(colnames(stage)[i])
  actual <- subset(stage,subset = !is.na(stage[,i]))
  tech <- paste(actual$Strain,actual$Stage,actual$Rep,sep = " ")
  if(length(unique(actual$Stage))<2) next
  if(length(unique(actual$Strain))<2){
    lme1 <- lmer(log(actual[,i])~actual$Stage + (1|tech))
    lme3 <- lmer(log(actual[,i])~(1|tech))
    em1 <-tryCatch({
      emmeans(lme1, ~ Stage)
    },error=function(cond){
      message(cond)
      return=NA
    })
    if(is.na(em1)) next
    pp <- summary(pairs(em1))
    pstage[rownames(pstage)==unique(actual$Strain),i-6] <- log(pp$p.value) - log(0.05)
    pstage[nlevels(vars)+1,i-6] <- log(anova(lme1,lme3)$`Pr(>Chisq)`[2])- log(0.05)
  }else{
    if(length(unique(paste(actual$Strain,actual$Stage)))==length(unique(actual$Strain))) next
    lme1 <- lmer(log(actual[,i])~actual$Strain * actual$Stage + (1|tech))
    lme2 <- lmer(log(actual[,i])~actual$Strain + actual$Stage + (1|tech))
    lme3 <- lmer(log(actual[,i])~actual$Strain  + (1|tech))
    em1 <-tryCatch({
      emmeans(lme1, ~ Stage | Strain)
    },error=function(cond){
      message(cond)
      return=NA
    })
    if(is.na(em1)) next
    pp <- summary(pairs(em1))
    pstage[,i-6] <- log(pp$p.value[match(rownames(pstage), pp$Strain)]) - log(0.05)
    pstage[nlevels(vars)+2,i-6] <- log(anova(lme1,lme2)$`Pr(>Chisq)`[2])- log(0.05)
    pstage[nlevels(vars)+1,i-6] <- log(anova(lme2,lme3)$`Pr(>Chisq)`[2])- log(0.05)
  }
}  

# Color scale
rc1 <- colorRampPalette(colors = c("gray", "gray"), space = "rgb")(5*abs(max(na.rm = T,pstage)))
rc2 <- colorRampPalette(colors = c("papayawhip","darkred","black"), space = "rgb")(5*abs(min(na.rm = T,pstage)))
rampcols_blank <- c(rc1, rc2)

# Heatmap
pdf("axen_medium_stats.pdf",width=12,height=3)
par(las=2,mar=c(2,12,7,2))
image(z=t(-1*pstage)[,nrow(pstage):1],xaxt="n",yaxt="n",col=rampcols_blank)
axis(side = 2,at=(1:nrow(pstage)-1)/(nrow(pstage)-1),labels = rev(rownames(pstage)))
axis(side = 3,at=(1:(ncol(stage)-6)-1)/((ncol(stage)-6)-1),labels = colnames(stage)[7:ncol(stage)])
box()
dev.off() 

# Legend
pdf("axen_medium_legend.pdf",width = 1.5,height = 2)
par(mar=c(0.5,4,0.5,0.1),las=1)
rs <- range(-1*pstage,na.rm=T)
lgd <- matrix((10*rs[1]):(10*rs[2])/10,nrow = 1)
image(lgd,xaxt="n",yaxt="n",col=rampcols_blank)
cap <- c(0.05,0.001,1e-6,1e-12,1e-18)
lsg <- (-log(cap)+log(cap[1])-rs[1])/(rs[2]-rs[1])
axis(side = 2,at=lsg,labels = prettyNum(cap,nsmall = 2))
box()
dev.off()

# Table
write.csv(exp(pstage+log(0.05)),file = "axen_medium.csv")
```








## STAGE COMPARISON BIOMASS

```{r}
fh <- read.table("hormones",header = T)
stage <- subset(fh,subset = (To_stages.comparison=="Y")&(Type=="biomass"))
vars <- as.factor(stage$Strain)

pstage <- matrix(NA,ncol=ncol(stage)-6,nrow=nlevels(vars)+2)

colnames(pstage) <- colnames(stage)[7:ncol(stage)]
rownames(pstage) <- c(unique(vars),as.factor(c("Stage:Strain","Stage")))

for(i in 7:ncol(stage)){
  print(colnames(stage)[i])
  actual <- subset(stage,subset = !is.na(stage[,i]))
  tech <- paste(actual$Strain,actual$Stage,actual$Rep,sep = " ")
  if(length(unique(actual$Stage))<2) next
  lme1 <- lmer(log(actual[,i])~actual$Strain * actual$Stage + (1|tech))
  #lme2 <- attr(step(lme1),"model")
  lme2 <- lmer(log(actual[,i])~actual$Strain + actual$Stage + (1|tech))
  lme3 <- lmer(log(actual[,i])~actual$Strain  + (1|tech))
  em1 <-tryCatch({
    emmeans(lme1, ~ Stage | Strain)
    },error=function(cond){
      message(cond)
      return=NA
    })
  if(is.na(em1)) next
  pp <- summary(pairs(em1))
  pstage[,i-6] <- log(pp$p.value[match(rownames(pstage), pp$Strain)]) - log(0.05)
  pstage[nlevels(vars)+2,i-6] <- log(anova(lme1,lme2)$`Pr(>Chisq)`[2])- log(0.05)
  pstage[nlevels(vars)+1,i-6] <- log(anova(lme2,lme3)$`Pr(>Chisq)`[2])- log(0.05)
}  

# Color scale
rc1 <- colorRampPalette(colors = c("gray", "gray"), space = "rgb")(5*abs(max(na.rm = T,pstage)))
rc2 <- colorRampPalette(colors = c("papayawhip","darkred","black"), space = "rgb")(5*abs(min(na.rm = T,pstage)))
rampcols_blank <- c(rc1, rc2)

# Heatmap
pdf("stage_biomass_stats.pdf",width=12,height=6)
par(las=2,mar=c(2,12,7,2))
image(z=t(-1*pstage)[,nrow(pstage):1],xaxt="n",yaxt="n",col=rampcols_blank)
axis(side = 2,at=(1:nrow(pstage)-1)/(nrow(pstage)-1),labels = rev(rownames(pstage)))
axis(side = 3,at=(1:(ncol(stage)-6)-1)/((ncol(stage)-6)-1),labels = colnames(stage)[7:ncol(stage)])
box()
dev.off() 

# Legend
pdf("stage_biomass_legend.pdf",width = 1.5,height = 2)
par(mar=c(0.5,4,0.5,0.1),las=1)
rs <- range(-1*pstage,na.rm=T)
lgd <- matrix((10*rs[1]):(10*rs[2])/10,nrow = 1)
image(lgd,xaxt="n",yaxt="n",col=rampcols_blank)
cap <- c(0.05,0.001,1e-6,1e-12,1e-18)
lsg <- (-log(cap)+log(cap[1])-rs[1])/(rs[2]-rs[1])
axis(side = 2,at=lsg,labels = prettyNum(cap,nsmall = 2))
box()
dev.off()

# Table
write.csv(exp(pstage+log(0.05)),file = "stage_biomass.csv")
```


## STAGE COMPARISON MEDIUM

```{r}
fh <- read.table("hormones",header = T)
stage <- subset(fh,subset = (To_stages.comparison=="Y")&(Type=="medium"))
vars <- as.factor(stage$Strain)

pstage <- matrix(NA,ncol=ncol(stage)-6,nrow=nlevels(vars)+2)

colnames(pstage) <- colnames(stage)[7:ncol(stage)]
rownames(pstage) <- c(unique(vars),as.factor(c("Stage:Strain","Stage")))

for(i in 7:ncol(stage)){
  print(colnames(stage)[i])
  actual <- subset(stage,subset = !is.na(stage[,i]))
  tech <- paste(actual$Strain,actual$Stage,actual$Rep,sep = " ")
  if(length(unique(actual$Stage))<2) next
  lme1 <-lmer(log(actual[,i])~actual$Strain * actual$Stage + (1|tech))
  lme2 <- lmer(log(actual[,i])~actual$Strain + actual$Stage + (1|tech))
  lme3 <- lmer(log(actual[,i])~actual$Strain  + (1|tech))
  em1 <-tryCatch({
    emmeans(lme1, ~ Stage | Strain)
  },error=function(cond){
    message(cond)
    return=NA
  })
  if(is.na(em1)) next
  pp <- summary(pairs(em1))
  pstage[,i-6] <- log(pp$p.value[match(rownames(pstage), pp$Strain)]) - log(0.05)
  pstage[nlevels(vars)+2,i-6] <- log(anova(lme1,lme2)$`Pr(>Chisq)`[2])- log(0.05)
  pstage[nlevels(vars)+1,i-6] <- log(anova(lme2,lme3)$`Pr(>Chisq)`[2])- log(0.05)
}  

# Color scale
rc1 <- colorRampPalette(colors = c("gray", "gray"), space = "rgb")(5*abs(max(na.rm = T,pstage)))
rc2 <- colorRampPalette(colors = c("papayawhip","darkred","black"), space = "rgb")(5*abs(min(na.rm = T,pstage)))
rampcols_blank <- c(rc1, rc2)

# Heatmap
pdf("stage__medium_stats.pdf",width=12,height=6)
par(las=2,mar=c(2,12,7,2))
image(z=t(-1*pstage)[,nrow(pstage):1],xaxt="n",yaxt="n",col=rampcols_blank)
axis(side = 2,at=(1:nrow(pstage)-1)/(nrow(pstage)-1),labels = rev(rownames(pstage)))
axis(side = 3,at=(1:(ncol(stage)-6)-1)/((ncol(stage)-6)-1),labels = colnames(stage)[7:ncol(stage)])
box()
dev.off() 

# Legend
pdf("stage__medium_legend.pdf",width = 1.5,height = 2)
par(mar=c(0.5,4,0.5,0.1),las=1)
rs <- range(-1*pstage,na.rm=T)
lgd <- matrix((10*rs[1]):(10*rs[2])/10,nrow = 1)
image(lgd,xaxt="n",yaxt="n",col=rampcols_blank)
cap <- c(0.05,0.001,1e-6,1e-12,1e-18)
lsg <- (-log(cap)+log(cap[1])-rs[1])/(rs[2]-rs[1])
axis(side = 2,at=lsg,labels = prettyNum(cap,nsmall = 2))
box()
dev.off()

# Table
write.csv(exp(pstage+log(0.05)),file = "stage_medium.csv")
```



