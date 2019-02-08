## =======================================================================
## PREAMBLE
## THE OPEN PROJECT (TOP)
## Vision
## Affordable and effective novel therapies discovered and developed based
## on all accessible, relevant data in a timely manner
## 
## 
## Mission
## Pioneer in translational modelling to develop, validate and improve 
## quantitative methods and tools for accurate experimental design to 
## enable robust decision making in drug discovery and development
## 
## 
## Participants
## Open: Any one can join for free to share data, models, codes and ideas
## Transparent: All results are properly documented to help the community
## Meritocracy: Participants are required to demonstrate understanding of
## the code, rules, and culture of the project before being invited to
## join
## Participation: Participation is by invitation only
## 
## 
## TOP was initiated by Tao You at Beyond Consulting Ltd.
## 
## 
## TRANSLATIONAL ONCOLOGY MODELLING GOALS
## Q: How do in vivo tumour models grow?
## This question breaks down into the following 3 questions.
## Goal 1. Which rate laws accurately recapitulate PDX tumour growth?
## Goal 2. Which rate laws accurately predict PDX tumour growth?
## Goal 3. Can in vivo gene expression predict in vivo tumour growth?
## 
## 
## THE OPEN PROJECT - TRANSLATIONAL ONCOLOGY MODELLING WIKI
## A wiki page exists on GitHub to answer some of the questions:
## https://github.com/TheOpenProject/TOM/wiki
## Why should I care?
## Who should contribute to TOP?
## Why contribute to The Open Project (TOP)?
## Vision & mission of TOP
## What is TOP doing?
## Ways to contribute to TOP
## Who contributes to TOP?
## The spirit of TOP
## 
## 
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
## Copyright © 2019 by Beyond Consulting Ltd.
## 
## Disclaimer
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
## Contact Information: 
## Please contact Tao You: Email: tao.you@letsgobeyond.co.uk 
## or visit www.letsgobeyond.co.uk.
## Created in RStudio Version 1.1.453 – © 2009-2019 RStudio, Inc. 
## 
## Beyond Consulting Ltd is a trusted advisor on translational oncology 
## PK/PD modelling. We help drug discovery and development projects to 
## design the right experiment, articulate results with all accessible 
## data, and generate compelling evidence by strengthening the link 
## between preclinical models and clinical data using proven innovative 
## modelling.
## =======================================================================

library(gdata)      # tested in version 2.18.0
library(nlme)       # tested in version 3.1-131
library(ggplot2)    # tested in version 2.2.1
library(propagate)  # tested in version 1.0-6
library(saemix)     # tested in version 2.1
library(deSolve)    # tested in version 1.20
library(reshape2)   # tested in version 1.4.3
library(readxl)     # tested in version 1.0.0
library(plyr)       # tested in version 1.8.4
library(dplyr)      # tested in version 0.7.8
library(DAAG)       # tested in version 1.22
library(minpack.lm) # tested in version 1.2-1
library(lattice)    # tested in version 0.20-35
library(GGally)     # tested in version 1.4.0

## SECTION 1. Tumour Volume modelling
## =========================================================================================================================
## Load Novartis PDX trial Tumour Volume (TV) data
# To run the R script, you would need to obtain the Novartis mouse clinical trial data from
# https://www.nature.com/articles/nm.3954#supplementary-information
# Excel files 
# Supplementary Table 1 
# Genomic profiling of PDXs and raw response and curve metrics of PCTs.
setwd('/appropriate directory to find the following file/')
dat = read_excel("Genomic profiling of PDXs and raw response and curve metrics of PCTs.xlsx",sheet="PCT raw data")
colnames(dat) = c("MODEL","TYPE","TREATMENT","TV","BW","TIME","DEL.TV","DEL.BW")
dat.ctrl = dat[which(dat$TREATMENT=='BKM120'),]
NumberOfAnimals = length(unique(dat.ctrl$MODEL))
dat.ctrl$ID = 0
for (i in 1:NumberOfAnimals){
  dat.ctrl$ID[which(dat.ctrl$MODEL == unique(dat.ctrl$MODEL)[i])] = i
}
dat.ctrl$TRAD = (3*dat.ctrl$TV/(4*pi))^(1/3)
dat.ctrl = as.data.frame(dat.ctrl)

duration = rep(0,NumberOfAnimals)
for (i in 1:NumberOfAnimals){
  duration[i] = max(dat.ctrl$TIME[dat.ctrl$ID == i])
}
DURATION = data.frame(DURATION=duration,TYPE=dat.ctrl$TYPE[dat.ctrl$TIME==0])
# histogram of duration
qplot(DURATION$DURATION,
      geom="histogram",
      binwidth = 5,
      xlab = "days",
      ylab = "Number of Animals",
      fill=I("blue"),
      col=I("black"), 
      alpha=I(.2),
      main = "Experimental Duration Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
           text = element_text(size = 15, family = "Tahoma"),
           axis.title = element_text(face="bold"),
           axis.text.x=element_text(size = 15)) 
# Histogram of duration by histology 
qplot(DURATION$DURATION,
      geom="histogram",
      binwidth = 5,
      xlab = "days",
      ylab = "Number of Animals",
      fill=factor(DURATION$TYPE),
      col=I("black"), 
      alpha=I(.5),
      main = "Experimental Duration Histogram") +
  theme(plot.title = element_text(size = 15, family = "Tahoma", face = "bold"),
        text = element_text(size = 15, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 15)) +
  guides(fill=guide_legend(title="Histology"))

# number of experiments: 224
length(unique(dat.ctrl$MODEL))
# number of histology: 224
length(unique(dat.ctrl$TYPE))
unique(dat.ctrl$TYPE) #  "GC"    "CRC"   "BRCA"  "NSCLC" "PDAC"  "CM"
length(unique(dat.ctrl$MODEL[dat.ctrl$TYPE=='GC']))     # 44
length(unique(dat.ctrl$MODEL[dat.ctrl$TYPE=='CRC']))    # 42
length(unique(dat.ctrl$MODEL[dat.ctrl$TYPE=='BRCA']))   # 39
length(unique(dat.ctrl$MODEL[dat.ctrl$TYPE=='NSCLC']))  # 29
length(unique(dat.ctrl$MODEL[dat.ctrl$TYPE=='PDAC']))   # 37
length(unique(dat.ctrl$MODEL[dat.ctrl$TYPE=='CM']))     # 33

# plot TV: Spaghetti plot
qplot(TIME, TV, 
      data=dat.ctrl,
      log="y",
      group=ID,
      col=factor(TYPE),
      xlim = c(0,100),
      geom='line',
      xlab="Time (Day)",ylab="TV (cubic mm)",main = "Individual TV") +
  labs(col="Histology") +
  theme_bw(base_size = 22) +
  guides(col=guide_legend(title="Histology"))


# plot Tumour Radius (TRAD): Spaghetti plot
qplot(TIME, TRAD, 
      data=dat.ctrl,
      group=ID,
      col=factor(TYPE),
      geom='line',
      xlim = c(0,100),
      xlab="Time (Day)",ylab="TRAD (mm)",main = "Individual TRAD") +
  # xlim(c(0,100)) +
  labs(col="Histology") +
  theme_bw(base_size = 22) +
  guides(col=guide_legend(title="Histology"))




############################################################################################################################
# NONLINEAR MIXED-EFFECTS MODELLING
# Construct grouped data for modelling
histology = unique(dat.ctrl$TYPE) # "GC"    "CRC"   "BRCA"  "NSCLC" "PDAC"  "CM"
dat.ctrl.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl)
dat.ctrl.GC.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl[which(dat.ctrl$TYPE %in% c("GC")),])
dat.ctrl.CRC.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl[which(dat.ctrl$TYPE %in% c("CRC")),])
dat.ctrl.BRCA.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl[which(dat.ctrl$TYPE %in% c("BRCA")),])
dat.ctrl.NSCLC.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl[which(dat.ctrl$TYPE %in% c("NSCLC")),])
dat.ctrl.PDAC.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl[which(dat.ctrl$TYPE %in% c("PDAC")),])
dat.ctrl.CM.grp<-groupedData(TV~TIME|MODEL, data=dat.ctrl[which(dat.ctrl$TYPE %in% c("CM")),])
# Plot individual tumour growth wrt each histology
plot(dat.ctrl.GC.grp, group = dat.ctrl.GC.grp$MODEL)
plot(dat.ctrl.CRC.grp, group = dat.ctrl.CRC.grp$MODEL)
plot(dat.ctrl.BRCA.grp, group = dat.ctrl.BRCA.grp$MODEL)
plot(dat.ctrl.NSCLC.grp, group = dat.ctrl.NSCLC.grp$MODEL)
plot(dat.ctrl.PDAC.grp, group = dat.ctrl.PDAC.grp$MODEL)
plot(dat.ctrl.CM.grp, group = dat.ctrl.CM.grp$MODEL)


############################################################################################################################
# GOAL 1. INFERENCE: WHICH RATE LAWS ACCURATELY RECAPITULATE PDX TUMOUR GROWTH?
# Method: Evaluate all models based on all data
# 1. Linear model
linearPopMod = nlme(TV~4/3*pi*(r0+g*TIME)^3,
                    fixed = r0+g~1,
                    random = pdDiag(r0+g~1), # This differs from linearPopMod1
                    data = dat.ctrl.grp, 
                    start = c(3,0.5),
                    method='ML',
                    control = nlmeControl(
                      pnlsMaxIter=10,
                      msMaxIter=100,
                      tolerance=1e-3)
                    )
linearPopMod1 = nlme(TV~4/3*pi*(r0+g*TIME)^3,
                     fixed = r0+g~1,
                     random = r0+g~1,  # This differs from linearPopMod
                     data = dat.ctrl.grp, 
                     start = c(3,0.5),
                     method='ML',
                     control = nlmeControl(
                       pnlsMaxIter=10,
                       msMaxIter=100,
                       tolerance=1e-3)
                     )
anova(linearPopMod,linearPopMod1) # p-value shows no significant difference
summary(linearPopMod) # logLik = -19776.3
# Fixed effects: r0 + g ~ 1 
#       Value  Std.Error   DF  t-value p-value
# r0 3.604944 0.04298111 3160 83.87275       0
# g  0.041170 0.00324775 3160 12.67640       0
plot(linearPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,
     xlim=c(-200,2500),ylim=c(-200,2500),pch=20,
     main="Linear Model")
cor(dat.ctrl.grp$TV,linearPopMod$fitted[,2]) # R Squared = 0.966
plot(augPred(linearPopMod),  level = 0:1, length.out = 2, ylim = c(0,1e3))

############################################
# generate P value by looking at all groups
anova(linearPopMod)
# look at 95% confidence intervals of parameter inference
intervals(linearPopMod, which = "fixed")     # fixed effects only
intervals(linearPopMod, which = "var-cov")   # random effects only
intervals(linearPopMod, which = "all")       # both effects

# diagnostic plots
# plot residuals with fitted values
plot(linearPopMod)
# plot residuals with data
plot(linearPopMod, resid(.)~TV, abline=0,na.rm=TRUE)
# examine residuals with time
plot(linearPopMod, resid(.)~TIME, abline=0,na.rm=TRUE)
# examine model predictions based on fixed effects only
plot(linearPopMod, fitted(.,0)~TV, abline=c(0,1),na.rm=TRUE)
# examine model predictions based on both fixed + random effects
plot(linearPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE)
# examine correlation among the random effects
plot(linearPopMod$coefficients$random$MODEL[,1],linearPopMod$coefficients$random$MODEL[,2],
     xlab="r0 random effects",
     ylab="g random effects",
     main="Linear Model")
grid()

# Pearson correlation test
# http://www.sthda.com/english/wiki/correlation-test-between-two-variables-in-r
cor.test(x = linearPopMod$coefficients$random$MODEL[,1],
         y = linearPopMod$coefficients$random$MODEL[,2],
         # use = "everything",
         method = "pearson")
# => p-value = 0.4577 The random effects are NOT significantly correlated

# plot individual fitting
plot(augPred(linearPopMod),  level = 0:1, length.out = 2, ylim = c(0,1e3))

# 2. Exponential
expPopMod = nlme(TV~TV0*exp(a0*TIME),
                 fixed = TV0+a0~1,
                 random = pdDiag(TV0+a0~1), # This differs from linearPopMod1
                 data = dat.ctrl.grp, 
                 start = c(150,0.01),
                 method='ML',
                 control = nlmeControl(
                   pnlsMaxIter=10,
                   msMaxIter=100,
                   tolerance=1e-3)
                 )
expPopMod1 = nlme(TV~TV0*exp(a0*TIME),
                  fixed = TV0+a0~1,
                  random = TV0+a0~1, # This differs from linearPopMod
                  data = dat.ctrl.grp, 
                  start = c(150,0.01),
                  method='ML',
                  control = nlmeControl(
                    pnlsMaxIter=10,
                    msMaxIter=100,
                    tolerance=1e-3)
                  )
anova(expPopMod,expPopMod1) # p-value (0.3761) shows no significant difference
summary(expPopMod) # logLik = -19541.76
# Fixed effects: TV0 + a0 ~ 1 
#         Value Std.Error   DF  t-value p-value
# TV0 219.19145  5.397353 3160 40.61091       0
# a0    0.02584  0.002044 3160 12.63863       0
plot(expPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,
     xlim=c(-200,2500),ylim=c(-200,2500),pch=20,
     main="Exponential Model")
cor(dat.ctrl.grp$TV,expPopMod$fitted[,2]) # R Squared = 0.9706353
# examine correlation among the random effects
plot(expPopMod$coefficients$random$MODEL[,1],expPopMod$coefficients$random$MODEL[,2],
     xlab="TV0 random effects",
     ylab="a0 random effects",
     main="Exponential Model")
grid()
cor.test(expPopMod$coefficients$random$MODEL[,1],expPopMod$coefficients$random$MODEL[,2])

# 3. Exponential-linear
expLinearPopMod = nlme(TV~(TIME<=tau)*TV0*exp(a0*TIME) +
                         (TIME>tau)*(exp(a0*tau)*TV0 + exp(tau*a0)*TV0*a0*(TIME-tau)),
                       fixed = TV0+a0+tau~1,
                       random = pdDiag(TV0+a0+tau~1), # This differs from expLinearPopMod1
                       data = dat.ctrl.grp, 
                       start = c(150,0.1,10),
                       method='ML',
                       control = nlmeControl(
                         pnlsMaxIter=10,
                         msMaxIter=100,
                         tolerance=1e-3)
                       )
expLinearPopMod1 = nlme(TV~(TIME<=tau)*TV0*exp(a0*TIME) +
                          (TIME>tau)*(exp(a0*tau)*TV0 + exp(tau*a0)*TV0*a0*(TIME-tau)),
                        fixed = TV0+a0+tau~1,
                        random = pdDiag(TV0+a0~1), # This differs from expLinearPopMod
                        data = dat.ctrl.grp, 
                        start = c(150,0.1,10),
                        method='ML',
                        control = nlmeControl(
                          pnlsMaxIter=10,
                          msMaxIter=100,
                          tolerance=1e-3)
                        )
expLinearPopMod2 = nlme(TV~(TIME<=tau)*TV0*exp(a0*TIME) +
                          (TIME>tau)*(exp(a0*tau)*TV0 + exp(tau*a0)*TV0*a0*(TIME-tau)),
                        fixed = TV0+a0+tau~1,
                        random = TV0+a0~1, # This differs from expLinearPopMod
                        data = dat.ctrl.grp, 
                        start = c(150,0.1,10),
                        method='ML',
                        control = nlmeControl(
                          pnlsMaxIter=10,
                          msMaxIter=100,
                          tolerance=1e-3)
                        )
anova(expLinearPopMod,expLinearPopMod1) # p-value shows no significant difference
anova(expLinearPopMod,expLinearPopMod1,expLinearPopMod2) # p-value shows no significant difference
summary(expLinearPopMod)
summary(expLinearPopMod1)
plot(expLinearPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,
     xlim=c(-200,2500),ylim=c(-200,2500),pch=20,
     main="Exponential-linear Model")
cor(dat.ctrl.grp$TV,expLinearPopMod$fitted[,2]) # R Squared = 0.970636
# examine correlation among the random effects
# TV0 - a0
plot(expLinearPopMod$coefficients$random$MODEL[,1],expLinearPopMod$coefficients$random$MODEL[,2],
     xlab="TV0 random effects",
     ylab="a0 random effects",
     main="Exponential-linear Model")
grid()
cor.test(expLinearPopMod$coefficients$random$MODEL[,1],expLinearPopMod$coefficients$random$MODEL[,2])
# TV0 - tau
plot(expLinearPopMod$coefficients$random$MODEL[,1],expLinearPopMod$coefficients$random$MODEL[,3],
     xlab="TV0 random effects",
     ylab="tau random effects",
     main="Exponential-linear Model")
grid()
cor.test(expLinearPopMod$coefficients$random$MODEL[,1],expLinearPopMod$coefficients$random$MODEL[,3])
# a0 - tau
plot(expLinearPopMod$coefficients$random$MODEL[,2],expLinearPopMod$coefficients$random$MODEL[,3],
     xlab="a0 random effects",
     ylab="tau random effects",
     main="Exponential-linear Model")
grid()
hist(expLinearPopMod$coefficients$random$MODEL[,3])
cor.test(expLinearPopMod$coefficients$random$MODEL[,2],expLinearPopMod$coefficients$random$MODEL[,3])


# 4. logistic growth
# This fitting works
logisticPopMod = nlme(TV~TV0*K/(TV0+(K-TV0)*exp(-a*TIME)),
                      fixed = TV0+K+a~1,
                      random = pdDiag(TV0~1), # This is unique and any other structure wouldn't fit
                      data = dat.ctrl.grp, 
                      start = c(150,2500,0.04),
                      method='ML',
                      control = nlmeControl(
                        pnlsMaxIter=10,
                        msMaxIter=100,
                        tolerance=1e-3)
                      )
# This fitting fails
logisticPopMod1 = nlme(TV~TV0*K/(TV0+(K-TV0)*exp(-a*TIME)),
                       fixed = TV0+K+a~1,
                       random = pdDiag(TV0+a~1), # This is different from logisticPopMod
                       data = dat.ctrl.grp, 
                       start = c(150,2000,0.02),
                       method='ML',
                       control = nlmeControl(
                         pnlsMaxIter=20,
                         msMaxIter=100,
                         maxIter = 100,
                         tolerance=1e-3)
                       )

# Sample in log scale: This fitting fails
logisticPopMod = nlme(TV~exp(TV0)*exp(K)/(exp(TV0)+(exp(K)-exp(TV0))*exp(-exp(a)*TIME)),
                      fixed = TV0+K+a~1,
                      random = pdDiag(TV0+a~1),
                      data = dat.ctrl.grp, 
                      start = c(log(150),log(5000),log(0.1)),
                      method='ML',
                      control = nlmeControl(
                        pnlsMaxIter=10,
                        msMaxIter=100,
                        tolerance=1e-3)
                      )

anova(logisticPopMod)
summary(logisticPopMod)
plot(logisticPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,
     xlim=c(-200,2500),ylim=c(-200,2500),pch=20,
     main="Logistic Model")
cor(dat.ctrl.grp$TV,logisticPopMod$fitted[,2]) # R Squared = 0.9021192
# examine correlation among the random effects
# TV0 - a0
hist(logisticPopMod$coefficients$random$MODEL[,1],xlab="TV0 random effects",
     main="Logistic Model")
grid()

# Fit logistic model using SAEMIX
PDX <- saemixData(name.data = dat.ctrl, name.group = "MODEL",
                  name.predictors = "TIME", name.response = "TV",
                  units = list(x = "Days", y = "uL"))
plot(PDX)
grid()

logistic.model <- function(psi, id, xidep) {
  TIME <- xidep[, 1]
  TV0 <- psi[id, 1]
  K <- psi[id, 2]
  a <- psi[id, 3]
  resp <- TV0*K/(TV0+(K-TV0)*exp(-a*TIME))
  return(resp)
}
PDX.model <- saemixModel(model = logistic.model,
                         description = "Logistic growth", 
                         psi0 = matrix(c(150, 1000, 0.1),ncol = 3, byrow = TRUE, 
                                       dimnames = list(NULL, c("TV0", "K","a"))),
                         covariance.model = matrix(c(1, 0, 0, 
                                                     0, 1, 0,
                                                     0, 0, 1),
                                                   ncol = 3, byrow = 3))
options = saemixControl(seed = 94352514, save = FALSE, save.graphs = FALSE,
                        nbiter.saemix = c(3000,500))
PDX.fit <- saemix(PDX.model, PDX, options)
plot(PDX.fit,plot.type="observations.vs.predictions")
plot(PDX.fit,plot.type="random.effects")
plot(PDX.fit,plot.type="correlations")
cor(dat.ctrl.grp$TV,fitted(PDX.fit,type=c("ipred"))) # R Squared = 0.9803193

# 5. Gompertz model
GompertzPopMod = nlme(TV~TV0*exp(alpha/beta*(1-exp(-beta*TIME))),
                      fixed = TV0+alpha+beta~1,
                      random = pdDiag(TV0+alpha+beta~1),
                      data=dat.ctrl.grp,
                      start=c(TV0 = 150,
                              alpha = 0.01,
                              beta = 0.1),
                      method='ML',
                      control = nlmeControl(
                        pnlsMaxIter=10,
                        msMaxIter=100,
                        tolerance=1e-3)
                      )

GompertzPopMod1 = nlme(TV~TV0*exp(alpha/beta*(1-exp(-beta*TIME))),
                       fixed = TV0+alpha+beta~1,
                       random = TV0+alpha+beta~1,
                       data=dat.ctrl.grp,
                       start=c(TV0 = 150,
                               alpha = 0.01,
                               beta = 0.1),
                       method='ML',
                       control = nlmeControl(
                         pnlsMaxIter=10,
                         msMaxIter=100,
                         tolerance=1e-3)
                       )
# log space to ensure inferred values are positive - this inference fails unfortunately
GompertzPopMod = nlme(TV~TV0*exp(exp(alpha)/exp(beta)*(1-exp(-exp(beta)*TIME))),
                      fixed = TV0+alpha+beta~1,
                      random = pdDiag(TV0+alpha+beta~1),
                      data=dat.ctrl.grp,
                      start=c(TV0 = 150,
                              alpha = log(0.01),
                              beta = log(0.1)),
                      method='ML',
                      control = nlmeControl(
                        pnlsMaxIter=10,
                        msMaxIter=100,
                        tolerance=1e-3)
)


anova(GompertzPopMod,GompertzPopMod1) # p-value shows GompertzPopMod1 is significantly better
summary(GompertzPopMod)
summary(GompertzPopMod1)
plot(GompertzPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,main="Gompertz Base Model")
plot(GompertzPopMod1, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,
     xlim=c(-200,2500),ylim=c(-200,2500),pch=20,
     main="Gompertz Model")
cor(dat.ctrl.grp$TV,GompertzPopMod1$fitted[,2]) # R Squared = 0.9781097
# examine correlation among the random effects
plot(GompertzPopMod$coefficients$random$MODEL[,1],GompertzPopMod$coefficients$random$MODEL[,2])
plot(GompertzPopMod$coefficients$random$MODEL[,1],GompertzPopMod$coefficients$random$MODEL[,3])
plot(GompertzPopMod$coefficients$random$MODEL[,2],GompertzPopMod$coefficients$random$MODEL[,3])
plot(GompertzPopMod1$coefficients$random$MODEL[,1],GompertzPopMod1$coefficients$random$MODEL[,2],
     xlab="TV0 random effect",
     ylab="alpha random effects",
     main="Gompertz Model")
plot(GompertzPopMod1$coefficients$random$MODEL[,1],GompertzPopMod1$coefficients$random$MODEL[,3],
     xlab="TV0 random effect",
     ylab="beta random effects",
     main="Gompertz Model")
plot(GompertzPopMod1$coefficients$random$MODEL[,2],GompertzPopMod1$coefficients$random$MODEL[,3],
     xlab="alpha random effect",
     ylab="beta random effects",
     main="Gompertz Model")
cor.test(GompertzPopMod$coefficients$random$MODEL[,1],GompertzPopMod$coefficients$random$MODEL[,2])
cor.test(GompertzPopMod$coefficients$random$MODEL[,1],GompertzPopMod$coefficients$random$MODEL[,3])
cor.test(GompertzPopMod$coefficients$random$MODEL[,2],GompertzPopMod$coefficients$random$MODEL[,3])
cor.test(GompertzPopMod1$coefficients$random$MODEL[,1],GompertzPopMod1$coefficients$random$MODEL[,2])
cor.test(GompertzPopMod1$coefficients$random$MODEL[,1],GompertzPopMod1$coefficients$random$MODEL[,3])
cor.test(GompertzPopMod1$coefficients$random$MODEL[,2],GompertzPopMod1$coefficients$random$MODEL[,3])
plot(augPred(GompertzPopMod),  level = 0:1, length.out = 2, ylim = c(0,1e3))


# Fit Gompertz model using SAEMIX: Sample in log space to sure inferred values are positive
PDX <- saemixData(name.data = dat.ctrl, name.group = "MODEL",
                  name.predictors = "TIME", name.response = "TV",
                  units = list(x = "Days", y = "uL"))
Gompertz.model <- function(psi, id, xidep) {
  TIME <- xidep[, 1]
  TV0 <- psi[id, 1]
  alpha <- exp(psi[id, 2])
  beta <- exp(psi[id, 3])
  resp <- TV0*exp(alpha/beta*(1-exp(-beta*TIME)))
  return(resp)
}
PDX.model <- saemixModel(model = Gompertz.model,
                         description = "Gompertz growth",
                         psi0 = matrix(c(150, log(0.05), log(0.001)),ncol = 3, byrow = TRUE,
                                       dimnames = list(NULL, c("TV0", "alpha","beta"))),
                         covariance.model = matrix(c(1, 0, 0, 
                                                     0, 1, 0, 
                                                     0, 0, 1),
                                                   ncol = 3, byrow = 3))
options = saemixControl(seed = 94352514, save = FALSE, save.graphs = FALSE,
                        nbiter.saemix = c(3000,500))
PDX.fit <- saemix(PDX.model, PDX, options)
plot(PDX.fit,plot.type="observations.vs.predictions")
grid()
plot(PDX.fit,plot.type="random.effects")
plot(PDX.fit,plot.type="correlations")
cor(dat.ctrl.grp$TV,fitted(PDX.fit,type=c("ipred"))) # R Squared = 0.9689807




############################################################################################################################
# GOAL 1. INFERENCE: DOES TUMOUR GROWTH RATE DEPEND ON HISTOLOGY?
# Method: Evaluate all models based on all data
# 1. Linear model
linearHistologyPopMod = nlme(TV~4/3*pi*(r0+g*TIME+
                                          (TYPE=="CM")*g_CM*TIME+
                                          (TYPE=="CRC")*g_CRC*TIME+
                                          (TYPE=="GC")*g_GC*TIME+
                                          (TYPE=="NSCLC")*g_NSCLC*TIME+
                                          (TYPE=="PDAC")*g_PDAC*TIME)^3,
                             fixed = r0+g+g_CM+g_CRC+g_GC+g_NSCLC+g_PDAC~1,
                             random = pdDiag(r0+g+g_CM+g_CRC+g_GC+g_NSCLC+g_PDAC~1),
                             data = dat.ctrl.grp,
                             start = c(3,0.1,0,0,0,0,0),
                             method='ML',
                             control = nlmeControl(
                               pnlsMaxIter=10,
                               msMaxIter=100,
                               tolerance=1e-3)
                             )
anova(linearPopMod,linearHistologyPopMod) # linearHistologyPopMod is significantly better than linearPopMod
summary(linearHistologyPopMod)
plot(linearHistologyPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,
     xlim=c(-200,2500),ylim=c(-200,2500),pch=20,
     main="Linear Histology Model")
cor(dat.ctrl.grp$TV,linearHistologyPopMod$fitted[,2])^2 # R Squared = 0.9341562


linearHistologyPopMod1 = nlme(TV~4/3*pi*(r0+(TYPE=="BRCA")*g_BRCA*TIME+
                                           (TYPE=="CM")*g_CM*TIME+
                                           (TYPE=="CRC")*g_CRC*TIME+
                                           (TYPE=="GC")*g_GC*TIME+
                                           (TYPE=="NSCLC")*g_NSCLC*TIME+
                                           (TYPE=="PDAC")*g_PDAC*TIME)^3,
                              fixed = r0+g_BRCA+g_CM+g_CRC+g_GC+g_NSCLC+g_PDAC~1,
                              random = pdDiag(r0+g_BRCA+g_CM+g_CRC+g_GC+g_NSCLC+g_PDAC~1),
                              data = dat.ctrl.grp,
                              start = c(3,0.1,0.1,0.1,0.1,0.1,0.1),
                              method='ML',
                              control = nlmeControl(
                                pnlsMaxIter=10,
                                msMaxIter=100,
                                tolerance=1e-3)
                              )
anova(linearHistologyPopMod1) # all parameters are significant
summary(linearHistologyPopMod1)

############################################
# generate P value by looking at all groups
anova(linearHistologyPopMod)
# look at 95% confidence intervals of parameter inference
intervals(linearHistologyPopMod1, which = "fixed")     # fixed effects only
intervals(linearHistologyPopMod1, which = "var-cov")   # random effects only
intervals(linearHistologyPopMod1, which = "all")       # both effects

# fixed effects
fixed.effects = read.csv(text="fixed,lower,est,upper
r0,3.51576268,3.60077505,3.68578742
BRCA,0.01222056,0.02317984,0.03413913
CM,0.04179566,0.05867597,0.07555628
CRC,0.02926307,0.04346153,0.05766000
GC,0.02080859,0.03185626,0.04290393
NSCLC,0.03335803,0.06253474,0.09171144
PDAC,0.02743318,0.03881295,0.05019273")
fixed.effects = fixed.effects[-1,]
ggplot(fixed.effects, aes(x=est, y=fixed)) + 
  geom_point(aes(x=est, y=fixed)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlim(c(0,0.1)) +
  xlab("PDX Growth Rate (mm/day)") +
  ylab("Histology") + 
  ggtitle("PDX growth by histology") +
  theme_bw(base_size = 15) 

# Inference based on Exponential Model
expHistologyPopMod = nlme(TV~TV0*exp(a0*TIME+
                                       (TYPE=="CM")*a_CM*TIME+
                                       (TYPE=="CRC")*a_CRC*TIME+
                                       (TYPE=="GC")*a_GC*TIME+
                                       (TYPE=="NSCLC")*a_NSCLC*TIME+
                                       (TYPE=="PDAC")*a_PDAC*TIME),
                          fixed = TV0+a0+a_CM+a_CRC+a_GC+a_NSCLC+a_PDAC~1,
                          random = pdDiag(TV0+a0+a_CM+a_CRC+a_GC+a_NSCLC+a_PDAC~1),
                          data = dat.ctrl.grp, 
                          start = c(150,0.01,0,0,0,0,0),
                          method='ML',
                          control = nlmeControl(
                            pnlsMaxIter=10,
                            msMaxIter=100,
                            tolerance=1e-3)
                          )
anova(expHistologyPopMod) # all parameters are significant
summary(expHistologyPopMod)

expHistologyPopMod1 = nlme(TV~TV0*exp((TYPE=="BRCA")*a_BRCA*TIME+
                                       (TYPE=="CM")*a_CM*TIME+
                                       (TYPE=="CRC")*a_CRC*TIME+
                                       (TYPE=="GC")*a_GC*TIME+
                                       (TYPE=="NSCLC")*a_NSCLC*TIME+
                                       (TYPE=="PDAC")*a_PDAC*TIME),
                          fixed = TV0+a_BRCA+a_CM+a_CRC+a_GC+a_NSCLC+a_PDAC~1,
                          random = pdDiag(TV0+a_BRCA+a_CM+a_CRC+a_GC+a_NSCLC+a_PDAC~1),
                          data = dat.ctrl.grp, 
                          start = c(150,0.01,0.01,0.01,0.01,0.01,0.01),
                          method='ML',
                          control = nlmeControl(
                            pnlsMaxIter=10,
                            msMaxIter=100,
                            tolerance=1e-3)
                          )

intervals(expHistologyPopMod1, which = "all")       # both effects
# fixed effects
fixed.effects = read.csv(text="fixed,lower,est,upper
TV0,2.081270e+02,218.75723913,229.38745717
BRCA,6.295684e-03,0.01350908,0.02072248
CM,2.839123e-02,0.03907083,0.04975043
CRC,1.926188e-02,0.02743182,0.03560176
GC,1.272322e-02,0.01952255,0.02632188
NSCLC,1.945404e-02,0.03864721,0.05784039
PDAC,1.792234e-02,0.02434136,0.03076037")
fixed.effects = fixed.effects[-1,]
ggplot(fixed.effects, aes(x=est, y=fixed)) + 
  geom_point(aes(x=est, y=fixed)) +
  geom_errorbarh(aes(xmin=lower, xmax=upper)) +
  xlim(c(0,0.06)) +
  xlab("PDX Growth Constant (/day)") +
  ylab("Histology") + 
  ggtitle("PDX growth by histology") +
  theme_bw(base_size = 15) 






############################################################################################################################
# GOAL 2. PREDICTION: WHICH TUMOUR GROWTH RATE LAWS ACCURATELY PREDICT CONTROL PDX TUMOUR GROWTH?
# Step 1: Create the training (days 0~30) and test (days 31~60) data samples from original data.
# Step 2: For each PDX, train a model using training data and generate predictions
# Step 3: Evaluate prediction accuracy by root-mean-square error (RMSE)
# Step 4: Compare RMSE for different histologies
# 1. Linear Model: TRAD
# 2. Exponential: TV
# 3. Exponential-linear: TV
# 4. logistic growth: TV
# 5. Gompertz: TV
# training and predictions
dat.selected = dat.ctrl[dat.ctrl$ID %in% unique(dat.ctrl$ID[dat.ctrl$TIME>=60]),]
length(unique(dat.selected$ID))
linearRMSE = rep(0,length(unique(dat.selected$ID)))
expRMSE = rep(0,length(unique(dat.selected$ID)))
expLinearRMSE = rep(0,length(unique(dat.selected$ID)))
logisticRMSE = rep(0,length(unique(dat.selected$ID)))
GompertzRMSE = rep(0,length(unique(dat.selected$ID)))
linearRSE = rep(0,length(unique(dat.selected$ID)))
expRSE = rep(0,length(unique(dat.selected$ID)))
expLinearRSE.a0 = rep(0,length(unique(dat.selected$ID)))
expLinearRSE.a1 = rep(0,length(unique(dat.selected$ID)))
logisticRSE.a = rep(0,length(unique(dat.selected$ID)))
logisticRSE.K = rep(0,length(unique(dat.selected$ID)))
GompertzRSE.alpha = rep(0,length(unique(dat.selected$ID)))
GompertzRSE.beta = rep(0,length(unique(dat.selected$ID)))

length(unique(dat.selected$MODEL[dat.selected$TYPE=='GC']))     # 22
length(unique(dat.selected$MODEL[dat.selected$TYPE=='CRC']))    # 15
length(unique(dat.selected$MODEL[dat.selected$TYPE=='BRCA']))   # 19
length(unique(dat.selected$MODEL[dat.selected$TYPE=='NSCLC']))  # 11
length(unique(dat.selected$MODEL[dat.selected$TYPE=='PDAC']))   # 7
length(unique(dat.selected$MODEL[dat.selected$TYPE=='CM']))     # 1

for (i in 1:length(unique(dat.selected$ID))){
  print(i)
  dat.current = dat.selected[dat.selected$ID == unique(dat.selected$ID)[i],]
  # training data: days 0-30
  trainingDat = dat.current[dat.current$TIME<=30,]
  # test data: days 31-60
  testDat = dat.current[dat.current$TIME>30 & dat.current$TIME<=60,]
  
  # 1. linear model: TV
  TRAD0 = trainingDat$TRAD[trainingDat$TIME==0]
  linearMod = nlsLM(TV~4/3*pi*(TRAD0+g*TIME)^3,
                    data=trainingDat,
                    start=list(g=0.1),
                    model = TRUE)
  linearPred <- predict(linearMod, testDat)
  # linearPredRes = linearPred - testDat$TRAD
  linearRMSE[i] = sqrt(mean((linearPred - testDat$TV)^2, na.rm = TRUE))
  linearRSE[i] = summary(linearMod)$coefficients['g','Std. Error']/
    summary(linearMod)$coefficients['g','Estimate']*100
  
  # 2. Exponential: TV
  TV0 = trainingDat$TV[trainingDat$TIME==0]
  expMod = nlsLM(TV~TV0*exp(a0*TIME),
                 data=trainingDat,
                 start=list(a0 = 0.1),
                 # lower=c(1e-5,1e-5),
                 # upper=c(1,1),
                 # algorithm = "port",
                 model = TRUE)
  # summary(expMod)
  expPred <- predict(expMod, testDat)
  expRMSE[i] = sqrt(mean((expPred - testDat$TV)^2, na.rm = TRUE))
  expRSE[i] = summary(expMod)$coefficients['a0','Std. Error']/
    summary(expMod)$coefficients['a0','Estimate']*100
  
  # The 31-day data failed to construct exponential-linear model
  # # 3. Exponential-linear: TV
  # TV0 = trainingDat$TV[trainingDat$TIME==0]
  # expLinearMod = nlsLM(TV~(TIME<=1/a0*log(a1/(a0*TV0)))*TV0*exp(a0*TIME) +
  #                        (TIME>1/a0*log(a1/(a0*TV0)))*(a1/a0 + a1*(TIME-1/a0*log(a1/(a0*TV0)))),
  #                      data=trainingDat,
  #                      start=list(a0 = 0.1,
  #                                 a1 = 10),
  #                      # lower=c(1e-5,1e-5),
  #                      # upper=c(1,1e2),
  #                      # algorithm = "port",
  #                      model = TRUE)
  # expLinearMod = nlsLM(TV~(TIME<=tau)*TV0*exp(a0*TIME) +
  #                        (TIME>tau)*(exp(a0*tau)*TV0 + exp(tau*a0)*TV0*a0*(TIME-tau)),
  #                      data=dat.ctrl[dat.ctrl$ID == i,],
  #                      start=list(a0 = 0.1,
  #                                 tau = 10),
  #                      # lower=c(1e-5,1e-5),
  #                      # upper=c(1,1e2),
  #                      # algorithm = "port",
  #                      model = TRUE)
  # summary(expLinearMod)
  # expLinearPred <- predict(expLinearMod, testDat)
  # expLinearRMSE[i] = sqrt(mean((expLinearPred - testDat$TV)^2, na.rm = TRUE))
  # expLinearRSE.a0[i] = summary(expLinearMod)$coefficients['a0','Std. Error']/
  #   summary(expLinearMod)$coefficients['a0','Estimate']*100
  # expLinearRSE.a1[i] = summary(expLinearMod)$coefficients['a1','Std. Error']/
  #   summary(expLinearMod)$coefficients['a1','Estimate']*100
  
  # 4. logistic growth: TV
  TV0 = trainingDat$TV[trainingDat$TIME==0]
  logisticMod = nlsLM(TV~TV0*K/(TV0+(K-TV0)*exp(-a*TIME)),
                      data=trainingDat,
                      start=list(K=100,
                                 a=0.1),
                      lower=c(0,-1),
                      upper=c(1e5,1),
                      algorithm = "port",
                      model = TRUE)
  # summary(logisticMod)
  logisticPred <- predict(logisticMod, testDat)
  logisticRMSE[i] = sqrt(mean((logisticPred - testDat$TV)^2, na.rm = TRUE))
  logisticRSE.a[i] = summary(logisticMod)$coefficients['a','Std. Error']/
    summary(logisticMod)$coefficients['a','Estimate']*100
  logisticRSE.K[i] = summary(logisticMod)$coefficients['K','Std. Error']/
    summary(logisticMod)$coefficients['K','Estimate']*100
  
  
  # 5. Gompertz: TV
  TV0 = trainingDat$TV[trainingDat$TIME==0]
  GompertzMod = nlsLM(TV~TV0*exp(alpha/beta*(1-exp(-beta*TIME))),
                      data=trainingDat,
                      start=list(alpha = 0.01,
                                 beta = 0.1),
                      lower=c(-1,-1),
                      upper=c(1,1),
                      algorithm = "port",
                      model = TRUE)
  # summary(GompertzMod)
  GompertzPred <- predict(GompertzMod, testDat)
  GompertzRMSE[i] = sqrt(mean((GompertzPred - testDat$TV)^2, na.rm = TRUE))
  GompertzRSE.alpha[i] = summary(GompertzMod)$coefficients['alpha','Std. Error']/
    summary(GompertzMod)$coefficients['alpha','Estimate']*100
  GompertzRSE.beta[i] = summary(GompertzMod)$coefficients['beta','Std. Error']/
    summary(GompertzMod)$coefficients['beta','Estimate']*100
}

result = data.frame(linearRMSE = linearRMSE,
                    expRMSE = expRMSE,
                    logisticRMSE = logisticRMSE,
                    GompertzRMSE = GompertzRMSE,
                    linearRSE = linearRSE,
                    expRSE = expRSE,
                    logisticRSE.a = logisticRSE.a,
                    logisticRSE.K = logisticRSE.K,
                    GompertzRSE.alpha = GompertzRSE.alpha,
                    GompertzRSE.beta = GompertzRSE.beta)
result = cbind(result,unique(dat.selected[c('MODEL','TYPE')]))
# write.csv(result,"PDX Growth Model Fitting Result.csv")  # Save result in csv

unique(dat.selected[c('TYPE')])

hist(linearRMSE/expRMSE)
hist(linearRMSE/logisticRMSE)
hist(linearRMSE/GompertzRMSE)
hist(expRMSE/logisticRMSE)
hist(expRMSE/GompertzRMSE)

hist(log(linearRMSE/expRMSE))
hist(log(linearRMSE/logisticRMSE))
hist(log(linearRMSE/GompertzRMSE))
hist(log(expRMSE/logisticRMSE))
hist(log(expRMSE/GompertzRMSE))

mean(linearRMSE/expRMSE)       # 1.12
mean(linearRMSE/logisticRMSE)  # 1.38
mean(linearRMSE/GompertzRMSE)  # 1.21

# pair plot
pairs(~linearRMSE+expRMSE+logisticRMSE+GompertzRMSE,
      data=result, # result[is.finite(rowSums(result)),],
      xlim=c(0,1e3),
      ylim=c(0,1e3),
      main="RMSE Comparison for Tumour Growth Models")

# pair plot (histology)
super.sym <- trellis.par.get("superpose.symbol")
splom(result[is.finite(result$GompertzRMSE) & result$GompertzRMSE<1e3,][c(
  'linearRMSE','expRMSE','logisticRMSE','GompertzRMSE')],
  groups=result$TYPE,
  panel=panel.superpose,
  xlim=c(0,1e3),
  ylim=c(0,1e3),
  key=list(title="RMSE Comparison for Tumour Growth Models",
           columns=6,
           points=list(pch=super.sym$pch[1:6],
                       col=super.sym$col[1:6]),
           text=list(c("GC","CRC","BRCA","NSCLC","PDAC","CM"))))

# pair plot (histology) / density kernel plot
ggpairs(result[is.finite(result$GompertzRMSE) & result$GompertzRMSE<1e3,]
        [c('linearRMSE','expRMSE','logisticRMSE','GompertzRMSE','TYPE')],
        aes(colour = TYPE, alpha = 0.4)
)

# PDAC: All models perform equally
# CRC: Gompertz fits are very different from any other models; The other models gave similar accuracy
# For CRC data, Gompertz gave 
#   infinite RMSE's: PDX X-1173 and PDX X-5495
#   extremely large RMSE's: PDX X-1855 and PDX X-2846
#   relatively large RMSE: PDX X-1479
# Plot the fitting results to identify PDX's that gave the largest deviations for Gompertz
# RMSE Comparison for Tumour Growth Models: CRC 
pairs(~linearRMSE+expRMSE+logisticRMSE+GompertzRMSE,
      data=result[is.finite(result$GompertzRMSE) & result$TYPE=='CRC',],
      xlim=c(0,3e2),
      ylim=c(0,3e2),
      main="RMSE Comparison for Tumour Growth Models: CRC")
# Observations
# 1) For PDX X-1441, Gompertz fits particularly well but the other models don't
# 2) For PDX X-1479, Gompertz fits particularly poorly but the other models fit well
# Action: identify the two PDX's
result[result$TYPE=='CRC',][c('GompertzRMSE','MODEL')]
# Here we plot the data to find out why
qplot(TIME, TV,
      # Data that fitted well
      # data=dat.selected[dat.selected$TYPE == 'CRC' &
      #   !(dat.selected$MODEL %in% c('X-1173','X-1479','X-1855','X-2846','X-5495')),],
      # Data that did not fit well
      data=dat.selected[dat.selected$MODEL %in% c('X-1173','X-1479','X-1855','X-2846','X-5495'),],
      # log="y",
      group=ID,
      col=factor(MODEL),
      xlim = c(0,60),
      ylim = c(0,500),
      geom='line',
      xlab="Time (Day)",ylab="TV (cubic mm)",main = "Individual TV") +
  labs(col="MODEL") +
  theme_bw(base_size = 22) +
  guides(col=guide_legend(title="MODEL"))


# plot TV: PDAC and CRC
qplot(TIME, TV,
      data=dat.selected[dat.selected$TYPE %in% c('PDAC','CRC'),],
      log="y",
      group=ID,
      col=factor(TYPE),
      xlim = c(0,60),
      geom='line',
      xlab="Time (Day)",ylab="TV (cubic mm)",main = "Individual TV") +
  labs(col="Histology") +
  theme_bw(base_size = 22) +
  guides(col=guide_legend(title="Histology"))

# plot CRC data and model predictions: Gompertz does not fit (well)
CRC_GompertzDoesNotFit = dat.selected[dat.selected$MODEL %in% c('X-1173','X-1479','X-1855','X-2846','X-5495'),]
CRC_GompertzDoesNotFit = CRC_GompertzDoesNotFit[CRC_GompertzDoesNotFit$TIME<=60,] # trim to day 60
CRC_GompertzDoesNotFit$PRED = 0
for (i in 1:length(c('X-1173','X-1479','X-1855','X-2846','X-5495'))){
  print(i)
  dat.current = CRC_GompertzDoesNotFit[
    CRC_GompertzDoesNotFit$MODEL %in% c('X-1173','X-1479','X-1855','X-2846','X-5495')[i],]
  # training data: days 0-20 (the first 3 weeks)
  trainingDat = dat.current[dat.current$TIME<=30,]
  # test data: days 21-41 (the next 3 weeks)
  testDat = dat.current[dat.current$TIME>30 & dat.current$TIME<=60,]
  # 5. Gompertz: TV
  TV0 = trainingDat$TV[trainingDat$TIME==0]
  GompertzMod = nlsLM(TV~TV0*exp(alpha/beta*(1-exp(-beta*TIME))),
                      data=trainingDat,
                      start=list(alpha = 0.01,
                                 beta = 0.1),
                      lower=c(-1,-1),
                      upper=c(1,1),
                      algorithm = "port",
                      model = TRUE)
  # summary(GompertzMod)
  CRC_GompertzDoesNotFit$PRED[CRC_GompertzDoesNotFit$MODEL %in% unique(dat.current$MODEL)] <-
    predict(GompertzMod, dat.current)
}
# plot data (dot) and predictions (line)
ggplot(data=CRC_GompertzDoesNotFit) +
  geom_point(aes(x=TIME,y=TV,col=factor(MODEL))) +
  geom_line(aes(x=TIME,y=PRED,col=factor(MODEL)),size=1) +
  xlim(c(0,60)) +
  ylim(c(0,500)) +
  xlab("Time (Day)") +
  ylab("TV (cubic mm)") + 
  ggtitle("Individual TV") +
  theme_bw(base_size = 22) +
  guides(col=guide_legend(title="Model"))


# plot CRC data and model predictions: Gompertz fits well
CRC_GompertzFitsWell = dat.selected[
  dat.selected$TYPE == 'CRC' & !(dat.selected$MODEL %in% c('X-1173','X-1479','X-1855','X-2846','X-5495')),]
CRC_GompertzFitsWell = CRC_GompertzFitsWell[CRC_GompertzFitsWell$TIME<=60,] # trim to day 60
CRC_GompertzFitsWell$PRED = 0
for (i in 1:length(unique(CRC_GompertzFitsWell$MODEL))){
  print(i)
  dat.current = CRC_GompertzFitsWell[
    CRC_GompertzFitsWell$MODEL == unique(CRC_GompertzFitsWell$MODEL)[i],]
  # training data: days 0-20 (the first 3 weeks)
  trainingDat = dat.current[dat.current$TIME<=30,]
  # test data: days 21-41 (the next 3 weeks)
  testDat = dat.current[dat.current$TIME>30 & dat.current$TIME<=60,]
  # 5. Gompertz: TV
  TV0 = trainingDat$TV[trainingDat$TIME==0]
  GompertzMod = nlsLM(TV~TV0*exp(alpha/beta*(1-exp(-beta*TIME))),
                      data=trainingDat,
                      start=list(alpha = 0.01,
                                 beta = 0.1),
                      lower=c(-1,-1),
                      upper=c(1,1),
                      algorithm = "port",
                      model = TRUE)
  # summary(GompertzMod)
  CRC_GompertzFitsWell$PRED[CRC_GompertzFitsWell$MODEL %in% unique(dat.current$MODEL)] <-
    predict(GompertzMod, dat.current)
}
# plot data (dot) and predictions (line)
ggplot(data=CRC_GompertzFitsWell) +
  geom_point(aes(x=TIME,y=TV,col=factor(MODEL))) +
  geom_line(aes(x=TIME,y=PRED,col=factor(MODEL)),size=1) +
  xlim(c(0,60)) +
  ylim(c(0,500)) +
  xlab("Time (Day)") +
  ylab("TV (cubic mm)") + 
  ggtitle("Individual TV") +
  theme_bw(base_size = 22) +
  guides(col=guide_legend(title="Model"))



############################################################################################################################
# CAN GOMPERTZ MODEL RECAPITULATE THE CRC TIME COURSE AND MAKE GOOD INFERENCE?
# 5. Gompertz
GompertzPopMod = nlme(TV~TV0*exp(alpha/beta*(1-exp(-beta*TIME))),
                      fixed = TV0+alpha+beta~1,
                      random = pdDiag(TV0+alpha+beta~1), # this is different from GompertzPopMod1
                      data=dat.ctrl.CRC.grp,
                      start=c(TV0 = 150,
                              alpha = 0.01,
                              beta = 0.1),
                      method='ML',
                      control = nlmeControl(
                        pnlsMaxIter=10,
                        msMaxIter=100,
                        tolerance=1e-3)
)
GompertzPopMod1 = nlme(TV~TV0*exp(alpha/beta*(1-exp(-beta*TIME))),
                       fixed = TV0+alpha+beta~1,
                       random = TV0+alpha+beta~1, # this is different from GompertzPopMod
                       data=dat.ctrl.CRC.grp,
                       start=c(TV0 = 150,
                               alpha = 0.01,
                               beta = 0.1),
                       method='ML',
                       control = nlmeControl(
                         pnlsMaxIter=10,
                         msMaxIter=100,
                         tolerance=1e-3)
)
anova(GompertzPopMod,GompertzPopMod1) # p-value shows GompertzPopMod1 is significantly better
summary(GompertzPopMod)
plot(GompertzPopMod, fitted(.,1)~TV, abline=c(0,1),na.rm=TRUE,xlim=c(-200,2500),ylim=c(-200,2500),pch=20)
plot(augPred(GompertzPopMod),  level = 0:1, length.out = 2, ylim = c(0,1e3),pch=20)
# Conclusion: 
# 1) Gompertz model is able to recapitulate the CRC time course. 
# 2) It is worth noting that inference for beta was not successful, indicating the model was
# effectively reduced to exponential model
# Fixed effects: TV0 + alpha + beta ~ 1 
#           Value Std.Error  DF   t-value p-value
# TV0   220.11257 11.773994 514 18.694810  0.0000
# alpha   0.02677  0.004688 514  5.711187  0.0000
# beta   -0.00627  0.006702 514 -0.935571  0.3499
# 3) Random effects look OK
# intervals(GompertzPopMod, which = "var-cov")   # random effects only
# Random Effects:
#   Level: MODEL 
#                 lower        est.       upper
# sd(TV0)   56.13932038 71.86781163 92.00293686
# sd(alpha)  0.02064016  0.02779453  0.03742876
# sd(beta)   0.02348742  0.03347588  0.04771212

# Guarantee positive alpha and beta: Sample in the log scale
GompertzPopMod = nlme(TV~TV0*exp(exp(alpha)/exp(beta)*(1-exp(-exp(beta)*TIME))),
                      fixed = TV0+alpha+beta~1,
                      random = pdDiag(TV0+alpha+beta~1),
                      data=dat.ctrl.CRC.grp,
                      start=c(TV0 = 150,
                              alpha = log(0.01),
                              beta = log(0.1)),
                      method='ML',
                      control = nlmeControl(
                        pnlsMaxIter=10,
                        msMaxIter=100,
                        tolerance=1e-3)
)
GompertzPopMod1 = nlme(TV~TV0*exp(exp(alpha)/exp(beta)*(1-exp(-exp(beta)*TIME))),
                       fixed = TV0+alpha+beta~1,
                       random = TV0+alpha+beta~1, # this is different from GompertzPopMod but it wouldn't fit
                       data=dat.ctrl.CRC.grp,
                       start=c(TV0 = 150,
                               alpha = log(0.01),
                               beta = log(0.1)),
                       method='ML',
                       control = nlmeControl(
                         pnlsMaxIter=10,
                         msMaxIter=100,
                         tolerance=1e-3)
)
summary(GompertzPopMod)
