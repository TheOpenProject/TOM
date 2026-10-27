## The purpose of this R script is to 
## 1) Infer tumour growth rate under the control condition
## 2) Construct an ElasticNet regression model to explain the tumour growth rate
## 
## Data: Novartis mouse clinical trial
## 
## Gao et al. (2015)
## High-throughput screening using patient-derived tumor xenografts to predict 
## clinical trial drug response.
## Nat Med. 21(11):1318-25.  doi: 10.1038/nm.3954.
## 
## LICENCE
## TheOpenProject/TOM is licensed under the
## GNU General Public License v3.0
## 
## https://github.com/TheOpenProject/TOM/blob/master/LICENSE
## Permissions of this strong copyleft license are conditioned on 
## making available complete source code of licensed works and 
## modifications, which include larger works using a licensed work, 
## under the same license. Copyright and license notices must be 
## preserved. Contributors provide an express grant of patent rights.
## 
## COPYRIGHT © 2019 by Beyond Consulting Ltd.
## 
## DISCLAIMER
## The content of this presentation may be subject to alterations and 
## updates. Therefore, information expressed in this presentation may not 
## reflect the most up-to-date information, unless otherwise notified by 
## an authorised representative independent of this presentation. No 
## contractual relationship is created by this presentation by any person 
## unless specifically indicated by agreement in writing.
## 
## 
## Beyond Consulting Ltd. Company. Registration no. 09698685 Registered in
## Englandand Wales. Registered Office: 14 Tytherington Park Road, 
## Macclesfield, Cheshire, SK10 2EL. United Kingdom
## 
## 
## CONTACT INFORMATION 
## Please contact Tao You: Email: tao.you@letsgobeyond.co.uk 
## or visit www.letsgobeyond.co.uk.
## Created in RStudio Version 1.1.453 – © 2009-2019 RStudio, Inc. 

library(gdata)
library(nlme)
library(ggplot2)
library(propagate)
library(saemix)
library(reshape2)
library(readxl)
library(plyr)
library(dplyr)
library(DAAG)
library(minpack.lm)
library(lattice)
library(GGally)
# heat map
library(pheatmap)
library(glmnet)


## SECTION 1. Tumour Volume modelling
## =======================================================================
## Load Novartis TV data
setwd('~/Documents/Beyond/Consulting/Novartis Mouse Clinical Trial/')

# dat = read_excel("Novartis PDX data - selected.xlsx",sheet="PCT raw data")
# dat = read_excel("Genomic profiling of PDXs and raw response and curve metrics of PCTs.xlsx",sheet="PCT raw data")
dat = read_excel("RNAseq_fpkm and PCT raw data.xlsx",sheet="PCT raw data")
colnames(dat) = c("MODEL","TYPE","TREATMENT","TV","BW","TIME","DEL.TV","DEL.BW")
dat.ctrl = dat[which(dat$TREATMENT=='untreated'),]
NumberOfAnimals = length(unique(dat.ctrl$MODEL))
dat.ctrl$ID = 0
for (i in 1:NumberOfAnimals){
  dat.ctrl$ID[which(dat.ctrl$MODEL == unique(dat.ctrl$MODEL)[i])] = i
}
dat.ctrl$TRAD = (3*dat.ctrl$TV/(4*pi))^(1/3)
dat.ctrl = as.data.frame(dat.ctrl)

####################################################################
# Nonlinear Mixed-effects Modelling
# Histology
histology = unique(dat.ctrl$TYPE) # "GC"    "CRC"   "BRCA"  "NSCLC" "PDAC"  "CM"
dat.ctrl.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl) # while time course
# dat.ctrl.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl[dat.ctrl$TIME<28,]) # day 0-20
ggplot(data = dat.ctrl.grp) +
  geom_line(aes(x=TIME,y=TV,group=MODEL))
length(unique(dat.ctrl.grp$MODEL))

############################################
# 1. Linear model
linearPopMod = nlme(TV~4/3*pi*(r0+g*TIME)^3,
                    fixed = r0+g~1,
                    random = pdDiag(r0+g~1),
                    data = dat.ctrl.grp, 
                    start = c(3,0.5),
                    method='ML',
                    control = nlmeControl(
                      pnlsMaxIter=10,
                      msMaxIter=100,
                      tolerance=1e-3)
                    )
summary(linearPopMod) # successful inference
growth = data.frame(g=linearPopMod$coefficients$fixed[['g']] + linearPopMod$coefficients$random$MODEL[,2])
growth$MODEL = rownames(growth)
rownames(growth) = NULL
histology = unique(data.frame(MODEL=dat.ctrl$MODEL,TYPE=dat.ctrl$TYPE))
growth=merge(growth,histology,by = "MODEL",all.x = TRUE, all.y = FALSE)
IQR(growth$g) # 0.04976554
quantile(growth$g, probs = seq(0, 1, 0.25))
#          0%         25%         50%         75%        100% 
# -0.08354210  0.01121225  0.03023423  0.06097779  0.26063449 
# examine model predictions based on both fixed + random effects
plot(linearPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,xlim=c(-200,2500),ylim=c(-200,2500),pch=20)
# examine correlation among the random effects
plot(linearPopMod$coefficients$random$MODEL[,1],linearPopMod$coefficients$random$MODEL[,2])
# R Squared
cor(linearPopMod$coefficients$random$MODEL[,1], linearPopMod$coefficients$random$MODEL[,2])^2

# histogram of control PDX growth rates
# natural scale - all histology pooled together
qplot(growth$g,
      geom="histogram",
      xlab = "Growth (mm/day)",
      ylab = "Number of PDX's",
      fill=I("blue"),
      col=I("black"),
      alpha=I(.2),
      main = "Control PDX Growth Rate Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15))
# natural scale - colour by histology
qplot(growth$g,
      geom="histogram",
      xlab = "Growth (mm/day)",
      ylab = "Number of PDX's",
      # fill=I("blue"),
      fill=factor(growth$TYPE), 
      col=I("black"),
      alpha=I(.5),
      main = "Control PDX Growth Rate Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15)) +
  # theme_bw(base_size = 22) +
  guides(fill=guide_legend(title="Histology"))
# log2 scale
qplot(log2(growth$g+1),
      geom="histogram",
      xlab = "log2( Growth (mm/day) +1) ",
      ylab = "Number of PDX's",
      fill=I("blue"), 
      col=I("black"), 
      alpha=I(.2),
      main = "Control PDX Growth Rate Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15))

############################################
# Read PDX genomics information
############################################
# RNAseq.fpkm = read_excel("Genomic profiling of PDXs and raw response and curve metrics of PCTs.xlsx",sheet="RNAseq_fpkm")
RNAseq.fpkm = read_excel("RNAseq_fpkm and PCT raw data.xlsx",sheet="RNAseq_fpkm")
RNAseq.fpkm = as.data.frame(RNAseq.fpkm)
# str(RNAseq.fpkm)
RNAseq.trm = RNAseq.fpkm
row.names(RNAseq.trm) = unlist(RNAseq.trm[,1])
RNAseq.trm = RNAseq.trm[,-1]
# remove genes with zero expression across all PDX's
RNAseq.trm = RNAseq.trm[rowSums(RNAseq.trm)!=0,]
# sort the columns (PDX) by PDX names
RNAseq.trm = RNAseq.trm[,order(colnames(RNAseq.trm))]
# To enable comparison between different samples, the FPKM data need to be transformed into TPM
# convert FPKM (Fragments Per Kilobase Million) to TPM (Transcripts Per Kilobase Million)
# TPMi=( FPKMi / sum(FPKMj ) * 10^6
for (i in 1:ncol(RNAseq.trm)){
  RNAseq.trm[,i] = RNAseq.trm[,i]/sum(RNAseq.trm[,i])*10^6
}
# extract columns containing the controls
length(unique(growth$MODEL))
x = RNAseq.trm[colnames(RNAseq.trm) %in% growth$MODEL]


# # log2 transform the data
# # A) log2(x+1)
# x = log2(x+1)
# # B) log2(x)
# x = log2(x)
# x[x==-Inf]=NA
# x = x[complete.cases(x),] # genes without zero expression in any sample: 9346 genes
# # C) filter x<1 = 0
# x[x<1] = 0
# x = log2(x+1)
# # D) filter x<1 = 0
# x[x<1] = 0
# x = log2(x)
# x[x==-Inf]=NA
# x = x[complete.cases(x),] # genes without zero expression in any sample: 6480 genes
# E) log2(x+0.01)          min(x[x!=0]) = 0.01104447
x = log2(x+0.001)
# # F) filter x<0.1 = 0
# x[x<0.1] = 0
# x = log2(x)
# x[x==-Inf]=NA
x = x[complete.cases(x),] # genes without zero expression in any sample: 6480 genes

hist(colSums(x)/nrow(x))
quantile(colSums(RNAseq.trm[colnames(RNAseq.trm) %in% growth$MODEL]), probs = seq(0, 1, 0.25))
#       0%      25%      50%      75%     100% 
# 277721.4 427665.5 494392.8 593434.7 905430.7
hist(rowSums(x)/ncol(x))
quantile(rowSums(x), probs = seq(0, 1, 0.25))
#           0%          25%          50%          75%         100% 
# 8.444290e-05 2.399568e-01 1.751287e+00 3.680772e+00 1.164479e+01 
hist(unlist(x[1,]),breaks = seq(0,10,0.1))
# a = rowSums(x)/ncol(x)
# a = a[a<1000]
# hist(log2(a+1))

# distribution of total gene expression for each sample
qplot(colSums(RNAseq.trm[colnames(RNAseq.trm) %in% growth$MODEL]),
      geom="histogram",
      # xlab = expression(paste("Total Gene Expression for Each Sample: ", Sigma, "log2(gene_i+1)")),
      xlab = "Total Gene Expression (TPM)",
      ylab = "Number of PDX's",
      fill=I("blue"), 
      col=I("black"), 
      alpha=I(.2),
      main = "Total Gene Expression Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15))

qplot(colSums(x[colnames(x) %in% growth$MODEL]),
      geom="histogram",
      # xlab = expression(paste("Total Gene Expression for Each Sample: ", Sigma, "log2(gene_i+1)")),
      xlab = "Total Gene Expression (copy number)",
      ylab = "Number of PDX's",
      fill=I("blue"), 
      col=I("black"), 
      alpha=I(.2),
      main = "Total Gene Expression Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15))

# distribution of expression of each gene averaged over all samples
qplot(rowSums(x)/ncol(x),
      geom="histogram",
      # xlab = expression(paste(Sigma[Sample],"log2(mRNA"[Sample_i],"+0.01)/N")),
      xlab = expression(paste(Sigma[Sample],"log2(mRNA"[Sample_i],")/N")),
      ylab = "Number of Genes",
      fill=I("blue"), 
      col=I("black"), 
      alpha=I(.2),
      main = "Gene Expression Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15))

# distribution of expression of each gene in each sample
qplot(unlist(x),
      geom="histogram",
      # xlab = expression(paste("log2(mRNA"[Sample_i],"+0.01)")),
      xlab = expression(paste("log2(mRNA"[Sample_i],")")),
      ylab = "Number of Genes in all PDX's",
      fill=I("blue"), 
      col=I("black"), 
      alpha=I(.2),
      main = "Gene Expression Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15))

# extract growth rates for PDX with genomic information
y = growth[growth$MODEL %in% colnames(RNAseq.trm),]
# sort the rows (PDX) by PDX names
y = y[order(y$MODEL),]
# select PDX's that showed growth, not shrinkage
# y = y[y$g>0,]
# x = x[,y$MODEL]

# ensure x and y are in the same order
colnames(x) == y$MODEL
rownames(t(x)) == y$MODEL

####################################################
############   ElasticNet Regression    ############
####################################################
x = as.matrix(t(x))
y = y[c("g")]
y = as.matrix(y)


# model construction
fit = glmnet(x, y)
plot(fit,label = TRUE)
plot(fit, xvar = "lambda", label = TRUE)
plot(fit, xvar = "dev", label = TRUE)
# Print shows from left to right the number of nonzero coefficients (Df), 
# the percent (of null) deviance explained (%dev) and the value of ?? (Lambda). 
print(fit)
coef(fit)
Lambda = min(fit$lambda)
Lambda
coefficients = as.data.frame(as.matrix(coef(fit, s=Lambda))) # select the biggest gene signature
colnames(coefficients) = "value"
coefficients$genes = rownames(coefficients)
head(coefficients)
coefficients = coefficients[coefficients$value!=0,]
dim(coefficients)[1]
coefficients

# generate a heatmap for these genes
# select genes
myvars <- colnames(x) %in% coefficients$genes
signature.dat <- data.frame(x)[myvars]
# get meta data
signature = cbind(MODEL=row.names(signature.dat), signature.dat)
# generate labels
signature$label = row.names(signature.dat)
# select genes and label
signature = signature[,-1]
# melt data
dat = melt(signature, id.vars = "label")
# summary statistics
dat.stats = ddply(dat, ~ label+variable, summarise,
                  mean = mean(value),
                  SD = sd(value),
                  SEM = sd(value)/sqrt(length(value)))
# cast data for heatmap
dat.heatmap = dcast(dat.stats, variable~label, value.var = "mean")
dat.heatmap.geneNames = dat.heatmap[,1]
dat.heatmap = dat.heatmap[,-1]
rownames(dat.heatmap) = dat.heatmap.geneNames
# Draw heatmaps
pheatmap(dat.heatmap)

####################################################
# Do cross-valication only once
# cv.glmnet is the main function to do cross-validation here, along with various supporting methods such as plotting and prediction.
# set.seed(94352514)
cvfit = cv.glmnet(x, y, parallel = TRUE)
plot(cvfit)
cvfit$lambda.min # view the selected ??'s and the corresponding coefficients 
# coef(fit,s=cvfit$lambda.min) # We can obtain the actual coefficients at one or more ??'s within the range of the sequence - sparse matrix format
# coef(cvfit, s = "lambda.min") # This is equivalent to the previous command - sparse matrix format
geneSignature = as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
colnames(geneSignature) = "value"
geneSignature$genes = rownames(geneSignature)
head(geneSignature)
geneSignature = geneSignature[geneSignature$value!=0,]
dim(geneSignature)[1]
geneSignature

as.matrix(coef(cvfit, s = "lambda.min")) # non-sparse format
coef(cvfit, s = "lambda.1se") # This gives the most regularized model such that error is within one standard error of the minimum
predict(cvfit, newx = x[1:5,], s = "lambda.min") # Predictions can be made based on the fitted cv.glmnet object

####################################################
# Repeat cross-valication many times and select the best result
CVTimes = 10
lambda.mins = rep(0,CVTimes)
for (i in 1:CVTimes){
  assign(paste("cvfit", i, sep = ""), cv.glmnet(x, y, parallel = TRUE))
  lambda.mins[i] = eval(parse(text=paste("cvfit", i, "$lambda.min", sep = "")))
  print(i)
}
lambda.mins
min(lambda.mins)
indexMinLambda = match(min(lambda.mins),lambda.mins)
indexMinLambda
cvfitBest = eval(parse(text=paste("cvfit", indexMinLambda, sep = "")))
plot(cvfitBest)
geneSignature = as.data.frame(as.matrix(coef(cvfitBest, s = "lambda.min")))
colnames(geneSignature) = "value"
geneSignature$genes = rownames(geneSignature)
head(geneSignature)
geneSignature = geneSignature[geneSignature$value!=0,]
dim(geneSignature)[1]
geneSignature

head(geneSignature)
# select genes
myvars <- colnames(x) %in% geneSignature$genes
signature.dat <- data.frame(x)[myvars]

# distribution of expression of each gene averaged over all samples
qplot(colSums(signature.dat)/nrow(signature.dat),
      geom="histogram",
      # xlab = expression(paste(Sigma[Sample],"log2(mRNA"[Sample_i],"+0.01)/N")),
      xlab = expression(paste(Sigma[Sample],"log2(mRNA"[Sample_i],")/N")),
      ylab = "Number of Genes",
      fill=I("blue"), 
      col=I("black"), 
      alpha=I(.2),
      main = "Expression of Gene Signature") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15))


####################################################
# plot data vs. prediction
# plot(y,predict(fit,newx=x,s=cvfit$lambda.min))
plot(y,predict(fit,newx=x,s=min(lambda.mins)),
     xlab="Data",ylab="Pred",main="Nett Tumour Growth Rate (mm/day)",
     xlim=c(-0.1,0.3),ylim=c(-0.1,0.3),pch=20)
abline(0,1)
grid()
# ggplot
DataPred = data.frame(Data = as.matrix(y),
                      Pred = as.matrix(predict(fit,newx=x,s=min(lambda.mins))))
                      # Pred = as.matrix(predict(fit,newx=x,s=cvfit$lambda.min)))
colnames(DataPred) = c("Data","Pred")
ggplot(data=DataPred) +
  geom_point(aes(x = Data, y = Pred)) +
  geom_abline() +
  # xlim(-0.1,0.3) +
  scale_x_continuous(limits = c(-0.1,0.3)) +
  ylim(c(-0.1,0.3)) +
  ggtitle("Nett Tumour Growth Rate (mm/day)") +
  theme_bw(base_size = 15)
cor(DataPred$Data,DataPred$Pred,method = "pearson")^2 # r squared 0.7271345
