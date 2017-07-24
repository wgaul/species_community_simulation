#########################
## Simulation function - 
## create_species
##
## This script defines the 'create_species()' function as part of Willson's
## simulation package
##
## author: Willson Gaul
## created: 9 Jun 2017
## last modified: 24 July 2017
##
#########################

### create_species
# input: parameters defining species probability of occurrence and
#   detection, responses to env. variables, and other charateristics
# output: a 1-row dataframe of parameters for a single species

create_species <- function(occ.formula = NULL, 
                           det.formula = NULL, 
                           species.name = NULL, 
                           prob.occurrence = 1, prob.detection = 1) {
  # Construct a list defining the characteristics of each species
  # ARGS: 
  #   occ.formula: an object of class "formula": a symbolic description of the model to be fitted.  Of the same format as formulas specified to e.g. glm().
  #   det.formula: an object of class "formula".
  #   species.name: a character string
  #   prob.occurrence: a 0 to 1 baseline probability of the species occuring 
  #       at a site
  #   prob.detection: a 0 to 1 probability that the species will be detected 
  #       when present 
  # OUTPUT: a list containing the following elements:
  #   species.name: character string with the name of the species
  #   prob.observed.formula: a double right-hand formula describing covariates
  #       of detection and occupancy in that order for this species.
  #   occ.resp.function: a function giving the species's detection probability 
  #       as a function of the baseline prob.occurrence (intercept) plus the 
  #       coefficients and site-level values for occurrence covariates specified in occ.formula
  #   det.resp.function: a function giving the species's detection probability 
  #       as a function of the baseline prob.detection (intercept) plus the
  #       coefficients and site/visit/observer-level values for detection 
  #       covariates specified in det.formulat
  #   
  is.formula <- function(x) {
    if(class(x) == "formula") T else F
  }
  if(is.null(species.name) || !is.character(species.name)) {
    stop("Provide a character vector to species.name")
  }
  if(!is.numeric(prob.occurrence) || prob.occurrence < 0 || 
     prob.occurrence > 1) {
    stop("Provide a numeric value between 0 and 1 to prob.occurrence")
  }
  if(!is.numeric(prob.detection) || prob.detection < 0 ||
     prob.detection > 1) {
    stop("Provide a numeric value between 0 and 1 to prob.detection")
  }
  if(!is.formula(occ.formula) || !is.formula(det.formula)) {
    stop("Specify a formula for occurrence and detection probabilities.")
  }
  
  # construct function to determine species probability of occurrence
  occ_f <- function(prob = NULL, 
                    occ.formula = NULL) {
    ## TODO: this will require a change in set_species_presence
    # returns the probability of occurrence adjusted by environmental covariates
    # This function is called by 'set_species_presence()'
    # ARGS: prob: the baseline probability of occurrence
    #       occ.coefficients: a list of numeric coefficients for occupancy covs
    #       occ.term.vals = a list of numeric values for environmental covariates
    #           at a single site
    
    logit_prob <- logit(prob)
    
    
    prob <- logistic(logit_prob)
    
    # logit_prob <- logit(prob)
    # for (t in 1:length(occ.coefficients)) {
    #   coef <- occ.coefficients[[t]]
    #   varValue <- occ.term.vals[[t]]
    #   logit_prob <- logit_prob + coef * varValue
    # }
    # prob <- logistic(logit_prob)
    prob
  }
  
  # construct function to determine species probability of detection
  det_f <- function(prob = NULL, 
                     det.coefficients = NULL, 
                     det.term.vals = NULL) {
    ## TODO: this will require a change in ?sample_species()?
    # returns the probability of detection adjusted by covariates
    # this takes as input the values of all coefficients and all variables
    # This function is called by '???'
    # ARGS: prob: the baseline probability of detection
    #       det.coefficients: a list of numeric coefficients for detection covs
    #       det.term.vals = a list of numeric values for detection covariates
    #           for a single survey/sample
    logit_prob <- logit(prob)
    for (t in 1:length(occ.coefficients)) {
      coef <- occ.coefficients[[t]]
      varValue <- occ.term.vals[[t]]
      logit_prob <- logit_prob + coef * varValue
    }
    prob <- logistic(logit_prob)
    prob
  }
  
  # l will be the object returned at the end of the function
  l <- list(species.name = as.character(species.name), 
            base.prob.occurrence = as.numeric(prob.occurrence), 
            base.prob.detection = as.numeric(prob.detection), 
            occ.function = occ_f, 
            det.function = det_f)
  l # return the list of characteristics for this species
}