# Genetic association study: 
# Paper: The impact of genetic variations in the serotonergic system on symptom severity and 
# clinical outcome in functional neurological disorders - Weber et al., 2024

# Code written by Samantha Weber & Raquel Cruz. 

# Exemplary code with mock data. 

# This code uses SNPassoc package (https://cran.r-project.org/web/packages/SNPassoc/index.html) as described in Gonzalez et al., 2007 (https://doi.org/10.1093/bioinformatics/btm025)

DataPath = '/Users/samanthaweber/Documents/GitHub/GeneticVariations_FND'
setwd(DataPath)

# Load packages
library(SNPassoc)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(exactRankTests)
library(rstatix)
library(readxl) #for reading exel-files
library(plyr)

#Read in data
data <- read_excel("ExampleData.xlsx")
names(data)

#Create dichotomous variables for group, called AFF; double check if correct
data$outcome<-ifelse(data$outcome=='improved', 0,ifelse(data$outcome=='worse', 1,-99)) 
data$outcome<-ifelse(data$gender=='male', 0,ifelse(data$gender=='female', 1,-99)) 

#check mean age
data%>%
  group_by(outcome) %>%
  get_summary_stats(age, type = "mean_sd")


################################################################################
#### Quality control of data ####

#we normally check per row and per column for missing data. If > 95% is missing in a row/column, then the row/column is excluded. 

#-------------------------------------------------------------------------------
# Quality Control for Genetic Data
#-------------------------------------------------------------------------------

# QC per sample: SNPassoc has a graphic (more down) but no syntax...

names(data) #SDATOS = selected data
data$cmiss<-apply(data[,2:3],1,function(x) length(which(is.na(x))))
data$pmiss<-data$cmiss*100/10
head(data) # check column pmiss and cmiss

# remove colums/rows with more than 90% missing.
data <- data[which(data$pmiss<10),]
dim(data)# nothing to remove

## Quality control using Hardy-Weinberg Equilibrium (HWE) 
# FUNDAMENTAL: CREATE setupSNP:
# we first specify which columns contain the genotyping data, and then what alleles they are without separation

names(data)
myData_S<-setupSNP(data=data,colSNPs=2:3,sep="") #colSNPs = All columns that contain SNPs
myData_S # this is the object that will be referred to in most SNPassoc orders.

plotMissing(myData_S) #visualizes combination of SNPs x Subjects who are missing
class(myData_S)

#---------------------------------------------------------------------------------------
### Hardy-Weinberg Equilibrium Test ###
#---------------------------------------------------------------------------------------
# CALCULATE frequency of genotypes, alleles y do HWE test. 

# summary view of various aspects of SNP quality control

res<-summary(myData_S)# aqui summary da un resultado muy diferente, al hacerse sobre un objeto setupSNP

#---------------------------------------------------------------------------------------
#### ASSOCIATION ANALYSIS FOR ONE SNP ####
#---------------------------------------------------------------------------------------
#TPH2 rs4570625
#TPH1 rs1800532

#---------------------------------------------------------------------------------------
#### PART I MANUSCRIPT ####
#---------------------------------------------------------------------------------------

# Is there an association between genotype (outcome) and symptom severity?
# Symptom severity is a continuous variable, thus we use a linear regression
myData_S$gender<-as.factor(myData_S$gender)
myData_S$outcome<-as.factor(myData_S$outcome)

# Symptom Severity (SS)
asoc_SS<-WGassociation(severity~1,data=myData_S, model=c("codominant","log-additive","dominant","recessive"),genotypingRate=80)
asoc_SS # This will give you the p-values. 

# Add covariates
casoc_SS<-WGassociation(severity~1 + gender + age + bdi,data=myData_S, model=c("codominant","log-additive","dominant","recessive"),genotypingRate=80)
casoc_SS # This will give you the p-values. 

# Now we look at the models in detail
WGstats(asoc_SS)

# Or at individual SNPs. dif = beta of the model
linear_SS<-association(severity~rs4570625,data=myData_S)
linear_SS

clinear_SS<-association(severity~rs4570625 + gender + age + bdi,data=myData_S)
clinear_SS

#or subgroups
linear_SS<-association(severity~rs4570625,subset=(gender=='female'),data=myData_S)
linear_SS


#---------------------------------------------------------------------------------------
#### PART II MANUSCRIPT ####
#---------------------------------------------------------------------------------------

# Outcome: This is a binary variable, thus we use logistic regression

asoc_S1<-WGassociation(outcome~1,data=myData_S, model=c("codominant","log-additive","dominant","recessive"),genotypingRate=80)
asoc_S1
WGstats(asoc_S1)# to see the details of the regression

# Add covariates
casoc_S1<-WGassociation(outcome~1+gender+age+bdi,data=myData_S, model=c("codominant","log-additive","dominant","recessive"),genotypingRate=80)
casoc_S1 # more or less stable when adding sex and age
WGstats(casoc_S1)# para ver las regresiones detalladas

#Check interaction Term --> this I do outside of SNP assoc because it did not seem to be possible inside
#Recessive model
myData_S$birs4570625<-revalue(myData_S$rs4570625,c("G/G"="1", "G/T"="1", "T/T"="0"))
res_birs4570625<-glm(outcome~rs4570625*bdi+gender+age,data=myData_S) # Make here something that makes sense
summary(res_birs4570625)

#---------------------------------------------------------------------------------------
#### PART III MANUSCRIPT ####
#---------------------------------------------------------------------------------------


# TO TAKE THE DATA TO MDR: all the SNPs and in the last column the dependent variable to be analyzed 

data_MDR<-data
names(data_MDR)
myData_MDR<-setupSNP(data=data_MDR,colSNPs=2:3,sep="")

names(myData_MDR)
myData_MDR$outcome<-as.factor(myData_MDR$outcome)
head(myData_MDR)

#This table can be loaded into MDR.MDR can not handle (continuous) covariates. Only outcome and SNPs will be saved.Start with SNPs
# Outcome variable has to be last. 
MDR1<-myData_MDR[,c(10:11,6)]
head(MDR1)
write.table(MDR1, file="SNP_Outcome.txt", row.names=FALSE, quote=FALSE)


# Outcome ----------------------------------------------------------------------
# Multiple linear regression with Findings from MDR ---------------------------

## MODELO REGRESION MULTIPLE - no se hace con SNPassoc PERO usamos algunas funciones de la librer?a para codificar los SNPs antes
## hay que decidir qu? modelo se usa. Por ej. aditivo en todos los SNPS 

# All data
myData_S$ad_rs1800532<-additive(myData_S$rs1800532)
myData_S$ad_rs4570625<-additive(myData_S$rs4570625)

myData_S$d_rs1800532<-dominant(myData_S$rs1800532)
myData_S$d_rs4570625<-dominant(myData_S$rs4570625)

myData_S$r_rs1800532<-recessive(myData_S$rs1800532)
myData_S$r_rs4570625<-recessive(myData_S$rs4570625)

myData_S$cd_rs1800532<-codominant(myData_S$rs1800532)
myData_S$cd_rs4570625<-codominant(myData_S$rs4570625)

head(myData_S)



## MULTIPLE REGRESSION OUTCOME - log-additive model

# putting two SNPs - to make this work better/that it makes sense, we should use much more SNPs
reg_mul<-glm(data=myData_S, outcome~ad_rs1800532+ad_rs4570625, family="binomial")
summary(reg_mul)

# Forward/backwards selection
sel_mul<-step(reg_mul, direction="both")
summary(sel_mul)

# add covariates
creg_mul<-glm(data=myData_S, outcome~ad_rs1800532+ad_rs4570625+age+gender+bdi, family="binomial")
summary(creg_mul)

# Forward/backwards selection
csel_mul<-step(creg_mul, direction="both")
summary(csel_mul)

