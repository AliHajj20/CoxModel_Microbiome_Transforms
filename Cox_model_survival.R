# Load the phyloseq object from the saved RData file
load("biom.Rdata")
# Assign the loaded object to a new variable 'phy' and display its summary
phy <- biom
phy
#agglomerate all the features of the same genus in a single variable
phyg <- tax_glom(phy, taxrank = "Genus" )
phyg
# Filter out taxa that are absent (i.e., zero) in more than 95% of the samples
phyg_f = filter_taxa(phyg, function(x) (sum(x == 0)/length(x)) <= 0.95, TRUE)
phyg_f
# Loading coda4microbiome and preparing input data, survival time, status, and covariates
library(coda4microbiome)#this package does log-ratio transformation
# Set data
x = data.frame(t(otu_table(phyg_f)))#extracts the OTU matrix from phyloseq object and transform it to a dataframe
time = as.numeric(sample_data(phyg_f)$T1Dweek)
status = as.numeric(sample_data(phyg_f)$T1D)
covars = data.frame(sample_data(phyg_f)$Sex, # Gender
                    sample_data(phyg_f)$Antibiotics) # Antibiotics

# Fitting a penalized Cox regression model
set.seed(12345)
T1Dmodel <- coda_coxnet(as.matrix(x),#dataframe
                        time,#time to event
                        status,#event occurrence 0/1
                        covar=covars, #data frame with covariates which are variables that might influence the time to event,
                        lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                        alpha = 0.9, nfolds = 5,#elatic net parameter
                        showPlots = FALSE,#no plots
                        coef_threshold = 0.01)#minimum absolute value of the coefficient for a variable to be included in the model
T1Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T1Dmodel$taxa.num
T1Dmodel$taxa.name
T1Dmodel$`log-contrast coefficients`

#signature plot
T1Dmodel$`signature plot`
#survival curve
survcurve <- plot_survcurves(risk.score = T1Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve
#Relative abundance transformation:
library(survival)
set.seed(12345)
T2Dmodel <- abund_coxnet(as.matrix(x),#dataframe
                        time,#time to event
                        status,#event occurrence 0/1
                        covar=covars, #data frame with covariates which are variables that might influence the time to event,
                        lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                        alpha = 0.9, nfolds = 5,#elatic net parameter
                        showPlots = FALSE,#no plots
                        coef_threshold = 0.01)#minimum absolute value of the coefficient for a variable to be included in the model
T2Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T2Dmodel$taxa.num
T2Dmodel$taxa.name

#signature plot
T2Dmodel$`signature plot`
#survival curve
survcurve <- plot_survcurves(risk.score = T2Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve
