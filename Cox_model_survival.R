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
T1Dmodel$`risk_score_plot`

#signature plot
T1Dmodel$`signature plot`
#survival curve
survcurve <- plot_survcurves(risk.score = T1Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve
#A-ilr transformation
library(survival)
set.seed(12345)
T2Dmodel <- abund_coxnet2(phyg_f,#dataframe
                        time,#time to event
                        status,#event occurrence 0/1
                        covar=covars, #data frame with covariates which are variables that might influence the time to event,
                        lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                        alpha = 0.9, nfolds = 5,#elatic net parameter
                        showPlots = FALSE,#no plots
                        coef_threshold = 0.01, transform = "ilr")#minimum absolute value of the coefficient for a variable to be included in the model
T2Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T2Dmodel$taxa.num
T2Dmodel$taxa.name

#risk score
T2Dmodel$`risk_score_plot`
#signature plot
T2Dmodel$`signature_plot`
#survival curve
survcurve_1 <- plot_survcurves(risk.score = T2Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve_1
#B-alr transformation
set.seed(12345)
T3Dmodel <- abund_coxnet2(phyg_f,#dataframe
                          time,#time to event
                          status,#event occurrence 0/1
                          covar=covars, #data frame with covariates which are variables that might influence the time to event,
                          lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                          alpha = 0.9, nfolds = 5,#elatic net parameter
                          showPlots = FALSE,#no plots
                          coef_threshold = 0.01, transform = "alr")#minimum absolute value of the coefficient for a variable to be included in the model
T3Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T3Dmodel$taxa.num
T3Dmodel$taxa.name

#risk score
T3Dmodel$`risk_score_plot`
#signature plot
T3Dmodel$`signature_plot`
#survival curve
survcurve <- plot_survcurves(risk.score = T3Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve
#ilr transformation refers to what exactly?
set.seed(12345)
T4Dmodel <- abund_coxnet2(phyg_f,#dataframe
                          time,#time to event
                          status,#event occurrence 0/1
                          covar=covars, #data frame with covariates which are variables that might influence the time to event,
                          lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                          alpha = 0.9, nfolds = 5,#elatic net parameter
                          showPlots = FALSE,#no plots
                          coef_threshold = 0.01, transform = "ilr")#minimum absolute value of the coefficient for a variable to be included in the model
T4Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T4Dmodel$taxa.num
T4Dmodel$taxa.name

#risk score
T4Dmodel$`risk_score_plot`
#signature plot
T4Dmodel$`signature_plot`
#survival curve
survcurve <- plot_survcurves(risk.score = T2Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve
#relative abundance transformation
set.seed(12345)
T5Dmodel <- abund_coxnet2(phyg_f,#dataframe
                          time,#time to event
                          status,#event occurrence 0/1
                          covar=covars, #data frame with covariates which are variables that might influence the time to event,
                          lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                          alpha = 0.9, nfolds = 5,#elatic net parameter
                          showPlots = FALSE,#no plots
                          coef_threshold = 0.01, transform = "rel")#minimum absolute value of the coefficient for a variable to be included in the model
T5Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T5Dmodel$taxa.num
T5Dmodel$taxa.name

#risk score
T5Dmodel$`risk_score_plot`
#signature plot
T5Dmodel$`signature_plot`
#survival curve
survcurve <- plot_survcurves(risk.score = T5Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve

#clr transformation
set.seed(12345)
T6Dmodel <- abund_coxnet2(phyg_f,#dataframe
                          time,#time to event
                          status,#event occurrence 0/1
                          covar=covars, #data frame with covariates which are variables that might influence the time to event,
                          lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                          alpha = 0.9, nfolds = 5,#elatic net parameter
                          showPlots = FALSE,#no plots
                          coef_threshold = 0.01, transform = "clr")#minimum absolute value of the coefficient for a variable to be included in the model
T6Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T6Dmodel$taxa.num
T6Dmodel$taxa.name

#risk score
T6Dmodel$`risk_score_plot`
#signature plot
T6Dmodel$`signature_plot`
#survival curve
survcurve_6 <- plot_survcurves(risk.score = T6Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve_6
#no transformation
set.seed(12345)
T7Dmodel <- abund_coxnet2(phyg_f,#dataframe
                          time,#time to event
                          status,#event occurrence 0/1
                          covar=covars, #data frame with covariates which are variables that might influence the time to event,
                          lambda = "lambda.1se",#penalization parameter,controls the strength of the penalty applied to the model, controls the overfitting
                          alpha = 0.9, nfolds = 5,#elatic net parameter
                          showPlots = FALSE,#no plots
                          coef_threshold = 0.01, transform = "none")#minimum absolute value of the coefficient for a variable to be included in the model
T4Dmodel

# Extracting selected taxa and log-contrast coefficients from the fitted T1D model
T4Dmodel$taxa.num
T4Dmodel$taxa.name

#risk score
T7Dmodel$`risk_score_plot`
#signature plot
T7Dmodel$`signature_plot`
#survival curve
survcurve <- plot_survcurves(risk.score = T7Dmodel$risk.score,
                             time,
                             status,
                             strata.quantile = 0.5)
survcurve

