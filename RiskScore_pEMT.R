## ----setup, include=FALSE------------------------------------------------------------------------------
# Load libraries
rm(list=ls())
knitr::opts_chunk$set(echo = F, message = F, warning = F)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(cgdsr)
library(singscore)
library(tidyverse)
library(doBy)
library(survival)
library(survminer)
library(corrplot)
library(Hmisc)
library(ggpubr)
library(gtools)


## ----import dat, include=F-----------------------------------------------------------------------------
# Read Puram et al. gene signature
setwd("~/Desktop/R/3_GitHub/RiskScore_pEMT/")
Puram <- t(read.csv(file = "Puram_pEMT_signature.csv", header = F, sep = ";", stringsAsFactors = F))
pemt.genes <- Puram[,2]
pemt.genes <- gsub(" ", "", pemt.genes, fixed = TRUE)

# Create CGDS object
mycgds = CGDS("https://www.cbioportal.org/")
test(mycgds)

# Get list of cancer studies at server
a <- getCancerStudies(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[105,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[3,1]

# Extract transcriptome from TCGA of Common Puram signatur
gen.prof <- getProfileData(mycgds,pemt.genes,mygeneticprofile,mycaselist)

# Genes from pemt.genes Puram signature found in TCGA patients
table(colnames(gen.prof)%in%pemt.genes) 
found <- which(colnames(gen.prof)%in%pemt.genes)
gen.prof <- gen.prof[,found]
# -> 4 genes of gen.prof not found in pemt.genes. Most likely different naming. For simplicity, 4 genes will be removed from analysis.

# Final check by intersection
length(intersect(colnames(gen.prof),pemt.genes))

# Extracted gen profile data from TCGA complete?
table(complete.cases(gen.prof))


## ----survival data-------------------------------------------------------------------------------------
# Impport TCGA survival data
data.clin <-  readxl::read_xlsx("~/Desktop/R/3_GitHub/RiskScore_pEMT/TCGA_clinical.xlsx")
data.clin$ID <- data.clin$`SAMPLE ID`
data.clin$ID <- chartr('-', '.',data.clin$ID)
data.clin <- orderBy(~ID, data.clin)

# Name ID column to match data
gen.prof$ID <- rownames(gen.prof)
gen.prof <- orderBy(~ID, gen.prof)

"Do genetic profile IDs and survival data IDs match?"
table(data.clin$ID == gen.prof$ID)

# 5 years follow-up cut-off
"Adjust survival data to 5 years follow-up"
table(data.clin$OS_STATUS)
data.clin$OS_MONTHS[data.clin$OS_MONTHS >= 60] <- 60
data.clin$OS_STATUS[data.clin$OS_MONTHS == 60 ] <- "LIVING"
table(data.clin$OS_STATUS)

data.clin$OS_STATUS[data.clin$OS_STATUS == "DECEASED"] <- 1
data.clin$OS_STATUS[data.clin$OS_STATUS == "LIVING"] <- 0
data.clin$OS_STATUS <- as.numeric(data.clin$OS_STATUS)


# Remove HPV+ patients
"Remove HPV+ patients"
table(data.clin$`HPV STATUS`)
ix.HPV <- data.clin$`HPV STATUS`=="HPV-"
data.clin <- data.clin[ix.HPV,]
gen.prof <- gen.prof[ix.HPV,]
#gen.profs <- gen.profs[,colnames(gen.profs) %in% gen.prof$ID]
"After removal"
table(data.clin$`HPV STATUS`)
"Control ID match again"
table(data.clin$ID == gen.prof$ID)
"Control OS status of HPV- and make it numeric"
table(data.clin$OS_STATUS)
"0 = no event/ LIVING, 1 = event/ DECEASED"
table(data.clin$OS_STATUS)


## ----log2 transformation-------------------------------------------------------------------------------
# Randomly select 3 genes
set.seed(1)
n <- sample(1:ncol(gen.prof),3)
gen.prof.log <- log2(gen.prof[-which(colnames(gen.prof)=="ID")]+ 0.00001)
# Gene expr. before and after log2 transformation
rafalib::mypar(3,2)
qqnorm(scale(gen.prof[,n[1]]), main = colnames(gen.prof)[n[1]])
qqline(scale(gen.prof[,n[1]]),col="steelblue",lwd=2)
qqnorm(scale(gen.prof.log[,n[1]]), main = colnames(gen.prof.log)[n[1]])
qqline(scale(gen.prof.log[,n[1]]), col="steelblue",lwd=2)
qqnorm(scale(gen.prof[,n[2]]), main = colnames(gen.prof)[n[2]])
qqline(scale(gen.prof[,n[2]]), col="steelblue",lwd=2)
qqnorm(scale(gen.prof.log[,n[2]]), main = colnames(gen.prof.log)[n[2]])
qqline(scale(gen.prof.log[,n[2]]), col="steelblue",lwd=2)
qqnorm(scale(gen.prof[,n[3]]), main = colnames(gen.prof)[n[3]])
qqline(scale(gen.prof[,n[3]]), col="steelblue",lwd=2)
qqnorm(scale(gen.prof.log[,n[3]]), main = colnames(gen.prof.log)[n[3]])
qqline(scale(gen.prof.log[,n[3]]), col="steelblue",lwd=2)
gen.prof.log$ID <- gen.prof$ID
rafalib::mypar(1,1)


## ----corrmatrix----------------------------------------------------------------------------------------
# Correlation matrix gene expression in TCGA patients
cor.matrix <- rcorr(as.matrix(gen.prof.log[,-ncol(gen.prof.log)]), type = "spearman")
corrplot(cor.matrix$r, 
         method= "shade", type="full", order="hc", tl.cex = .3,
         p.mat = cor.matrix$P, 
         sig.level = 0.05, insig = "blank", tl.col = "black")


## ----Univariate CoxPH table----------------------------------------------------------------------------
"Control ID match again"
table(data.clin$ID==gen.prof.log$ID)
# Run a univariate Cox PH model for every gene
ALL <- data.frame("ID"=data.clin$ID,"Time.OS"=data.clin$OS_MONTHS,
                        "Status.OS"=data.clin$OS_STATUS,
                        gen.prof.log[,colnames(gen.prof.log) %in% pemt.genes])
unvCoxPHtab.ALL <- tibble()
i = 4
for (i in 4:ncol(ALL)) {
  temp <- (coxph(formula = Surv(time = Time.OS,event = Status.OS) ~ ALL[,i], data = ALL) %>% broom::tidy(exp=T) 
           %>% mutate(signif = stars.pval(p.value)))
  unvCoxPHtab.ALL <- rbind(unvCoxPHtab.ALL,temp)
  i = i + 1
}

# Show genes with p-val <= .1
unvCoxPHtab.ALL$term <- colnames(ALL[,4:ncol(ALL)])
unvCoxPHtab.ALL <- unvCoxPHtab.ALL[unvCoxPHtab.ALL$p.value <= 0.1,]
unvCoxPHtab.ALL %>% knitr::kable(digits = c(3,3,3,3,3,3))


## ----RB ALL--------------------------------------------------------------------------------------------
# Run Robust likelihood-based survival modeling 
# Use genes from univ. Cox PH models with p-val <= .1
library(rbsurv)
rb.fit <-rbsurv.default(time=ALL$Time.OS, status=ALL$Status.OS, x=t(ALL[,colnames(ALL) %in% unvCoxPHtab.ALL$term]), gene.ID = rownames(t(ALL[,colnames(ALL) %in% unvCoxPHtab.ALL$term])),set.seed=1)
rb.fit$model
"Genes selected by RB survival modeling"
rb.fit$gene.list 


## ----Visualiztion--------------------------------------------------------------------------------------
# Models visualized with ggforest
cox <- coxph(formula = Surv(time = Time.OS,event = Status.OS) ~ 
                            PLOD2+TPM4+TNFRSF6B+SLC38A5+ANXA5+CXCL14+CD63+P4HA2+C1S, data = ALL)

# Final multivariate Cox PH model
cox.sel <- coxph(formula = Surv(time = Time.OS,event = Status.OS) ~ 
                            SLC38A5+CXCL14+P4HA2+C1S, data = ALL)
ggarrange(ggforest(cox),ggforest(cox.sel), ncol = 1, nrow = 2)


## ----RB Corr-------------------------------------------------------------------------------------------
# Control for collinearity in final Cox PH model
cor.matrix2 <- rcorr(as.matrix(gen.prof.log[,colnames(gen.prof.log) %in% names(cox.sel$coefficients)]), type = "spearman")
corrplot(cor.matrix2$r, 
         method= "shade", type="full", order="hc", tl.cex = 1,
         p.mat = cor.matrix2$P, 
         sig.level = 0.05, insig = "p-value", tl.col = "black")



## ----Risk Score----------------------------------------------------------------------------------------
"Coefficients of genes forming Risk Score"
cox.sel$coefficients
# Simple loop to compute Risk Score
ALL.expr <- ALL[,colnames(ALL) %in% names(cox.sel$coefficients)]
names <- names(cox.sel$coefficients)
i= 1
for(i in 1:length(names)){
  ALL.expr[,names[i]] <- ALL.expr[,names[i]] * cox.sel$coefficients[names[i]]
  i = i+1
}
rm(names)
# Scale Risk Score (Z-scores)
ALL$RiskScore <- scale(rowSums(ALL.expr))
# Binarize Risk Score based on median
ALL$RiskBin <- c("Risk-","Risk+")[(ALL$RiskScore>median(ALL$RiskScore))+1]
# Add Primary Site to data set and rename factors
ALL$PrimSite <- as.factor(data.clin$`PRIMARY SITE`)
levels(ALL$PrimSite) <- list(Hyp="Hypopharynx",Lar="Larynx",OC="Oral Cavity",Oro="Oropharynx")

# Define Status vector 
Status <- as.factor(ALL$Status.OS)
levels(Status) <- list(LIVING="0",DECEASED="1")

# Quick Explorative Data Analysis of computed Risk Score
rafalib::mypar(2,2)
hist(ALL$RiskScore, main="Histogram of Risk Score", xlab = "Risk Score")
qqnorm(ALL$RiskScore, pch = 1, frame = FALSE, main = "QQ-plot of Risk Score")
qqline(scale(ALL$RiskScore), col = "steelblue", lwd = 2)
wilcox <- wilcox.test(ALL$RiskScore[Status=="LIVING"],ALL$RiskScore[Status=="DECEASED"])
plot(Status,ALL$RiskScore,xlab="",ylab="Risk Score",main=paste0("Wilcox test, p= ",signif(wilcox$p.value,3)))
fisher <- fisher.test(table(Status,ALL$RiskBin))
plot(Status,as.factor(ALL$RiskBin),xlab="",ylab="",main=paste0("Binarized: Fisher's exact p= ",signif(fisher$p.value,3)))
rafalib::mypar(1,1)


## ----Surv Risk-----------------------------------------------------------------------------------------
# Cox models with Risk Score before and after binarization
"Continous Cox PH model with Risk Score as feature"
summary(coxph(formula = Surv(time = Time.OS,event = Status.OS) ~ RiskScore,data = ALL))
"Binarized Risk Score visualized with a Kaplan Meier"
(cox.RS <- summary(coxph(formula = Surv(time = Time.OS,event = Status.OS) ~ RiskBin,data = ALL)))


# Extract values of Cox PH model and display in Kaplan Meier
di.cf <- cox.RS$coefficients
di.ci <-cox.RS$conf.int
hr.di <- round(di.cf[2],2) #hazard ratio
ci.di <- paste(round(di.ci[3],2),"-",round(di.ci[4],2),sep="") #95% CI
fit <- survfit(formula = Surv(time = Time.OS,event = Status.OS) ~ RiskBin,data = ALL)
KM_RiskScore <- survminer::ggsurvplot(fit = fit,conf.int = T,pval = F,risk.table = T, 
                                       legend.labs=c("Risk-","Risk+"), palette = c("chartreuse3","brown3"),
                                       surv.plot.height = 2/3, risk.table.height = 1/3)
KM_RiskScore$plot <- KM_RiskScore$plot+ 
  ggplot2::annotate("text", x = 20, y  = 0.1, # x and y coordinates of the text
                    label = paste("logrank p= ",signif(cox.RS$logtest[3],3),"/ HR:",hr.di,"(95% CI:",ci.di,")", sep = " "), 
                    size = 4)
KM_RiskScore


## ----RiskScore PrimSite--------------------------------------------------------------------------------
# Cox PH model of binarized Risk Score and Primary sites
fit <- survfit(formula = Surv(time = Time.OS,event = Status.OS) ~ RiskBin + PrimSite,data = ALL[ALL$PrimSite == "Lar"|ALL$PrimSite == "OC",])
(KM_RiskScore.PrimSite <- survminer::ggsurvplot(fit = fit,conf.int = F,pval = T,risk.table = T, 
                                       legend.labs=c("Risk-/Lar","Risk-/OC","Risk+/Lar","Risk+/OC"),
                                       surv.plot.height = 2/3, risk.table.height = 1/3))

# Pairwise logrank comparison
ps.PrimSite <- pairwise_survdiff(formula = Surv(time = Time.OS,event = Status.OS) ~ RiskBin + PrimSite,data = ALL[ALL$PrimSite == "Lar"|ALL$PrimSite == "OC",])
symnum(ps.PrimSite$p.value, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
   symbols = c("****", "***", "**", "*", "+", "ns"),
   abbr.colnames = F, na = "") %>% knitr::kable()


## ----DFS-----------------------------------------------------------------------------------------------
# Rearrange DFS data
DFS.Status <- data.clin$DFS_STATUS[data.clin$DFS_STATUS!="NA"]
DFS.Status[DFS.Status =="DiseaseFree"] <- 0
DFS.Status[DFS.Status =="Recurred/Progressed"] <- 1
DFS.Status <- as.numeric(DFS.Status)

DFS.Time <- as.numeric(data.clin$DFS_MONTHS[data.clin$DFS_MONTHS!="NA"])
RiskBin.DFS <-  ALL$RiskBin[data.clin$DFS_STATUS!="NA"]
RiskScore <- ALL$RiskScore[data.clin$DFS_STATUS!="NA"]

# 3 years follow-up cut-off
table(DFS.Status)
DFS.Time[DFS.Time >= 36] <- 36
DFS.Status[DFS.Time== 36] <- 0

# Cox PH model of DFS and Kaplan Meier 
DFS <- data.frame(DFS.Time,DFS.Status,RiskBin.DFS)
summary(coxph(formula = Surv(time = DFS.Time,event = DFS.Status) ~ RiskScore))
cox.RS <- summary(coxph(formula = Surv(time = DFS.Time,event = DFS.Status) ~ RiskBin.DFS))
fit <- survfit(formula = Surv(time = DFS.Time,event = DFS.Status) ~ RiskBin.DFS,data = DFS)
di.cf <- cox.RS$coefficients
di.ci <-cox.RS$conf.int
hr.di <- round(di.cf[2],2) 
ci.di <- paste(round(di.ci[3],2),"-",round(di.ci[4],2),sep="") 
KM_RiskScore.DFS <- survminer::ggsurvplot(fit = fit,conf.int = T,pval = F,risk.table = T, 
                                       legend.labs=c("Risk-","Risk+"), palette = c("chartreuse3","brown3"),
                                       surv.plot.height = 2/3, risk.table.height = 1/3)
KM_RiskScore.DFS$plot <- KM_RiskScore.DFS$plot+ 
  ggplot2::annotate("text", x = 15, y  = 0.1, 
                    label = paste("logrank p= ",signif(cox.RS$logtest[3],3),"/ HR:",hr.di,"(95% CI:",ci.di,")", sep = " "), 
                    size = 4)
KM_RiskScore.DFS

