# GeneticVariations_FND

Here we provide the explanatory code for the manuscript "The impact of genetic variations in the serotonergic system on symptom severity and clinical outcome in functional neurological disorders - Weber et al., 2024"

The paper contains three analyses. The code here provides example data with the corresponding code in order to facilitate reproducibility of our analyses. For further questions, feel free to contact: samantha.weber@bli.uzh.ch

## 1. Linear regression Symptom severity ~ Genotype + Cov
For the first analysis, a linear regression with the individual SNPs as outcome (independent) variable symptom severity (CGI1) as predictor (dependent) variable was performed using SNPassoc package in R, including age, gender, depression (BDI), trait-anxiety (STAI) and total CTQ scores as covariates. 

## 2. Logistic regression Outcome (good/bad) ~ Genotype + Cov
A logistic regression was performed using SNPassoc package using clinical outcome (dichotomized; improved vs. no change/worsened) as dependent variable, including. For all models, age, gender, depression, trait-anxiety and total CTQ score were included as covariates. 

## 3. Multiple linear regression
The combined effect of multiple SNPs was calculated with the use of multiple logistic/linear regression models using a stepwise selection procedure. This procedure systematically added and removed predictors based on their statistical significance to identify the model that best explains the in order to identify the best interaction model explaining the interaction between SNPs and 1) symptom severity and 2) clinical outcome. 

## 4. MDRI (multifactorial dimensionality reduction)
The combined effect and interaction between SNPs and binary dependent variables (i.e., clinical outcome) were also tested using the multifactorial dimensionality reduction system (MDR; https://sourceforge.net/projects/mdr/files/mdrdt/), in order to define the best interaction models for our analyses. MDR reduces high-dimensional multi-locus genetic data into a single dimension. MDR is optimized to detect gene-gene or gene-environment interactions in (binary) case-control studies. Thus, MDR can be used to identify interactions even in the absence of a main effect of selected genes. A 2-to-3 interaction model was considered, and a 10-fold cross validation (CV) procedure was applied. 


