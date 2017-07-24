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
#########################

### create_species
# input: parameters defining species probability of occurrence and
#   detection, responses to env. variables, and other charateristics
# output: a 1-row dataframe of parameters for a single species

create_species <- function(species.name = NULL, 
                           prob.occurrence = 1, prob.detection = 1, 
                           occ.coefficients = NULL, occ.terms = NULL, 
                           det.coefficients = NULL, det.terms = NULL) {
  # Construct a list defining the characteristics of each species
  # ARGS: 
  #   species.name: a character string
  #   prob.occurrence: a 0 to 1 baseline probability of the species occuring 
  #       at a site
  #   prob.detection: a 0 to 1 probability that the species will be detected 
  #       when present 
  #   occ.coefficients: an optional vector of coefficients for each occ.term 
  #       in the species's occurrence response function. 
  #       This coefficient will be on the logit 
  #       scale, so that this coefficient value can be directly compared to the 
  #       parameter estimate from a logistic regression model.  
  #       e.g. providing a value of 0.03 here and a value of "elevation" to the
  #       resp.terms argument would produce a species for which the log of 
  #       the odds increases by 0.03 for each 1 unit increase in elevation.
  #       Providing c(0.03, 0.01) here and c("elevation", "elevation^2") in 
  #       resp.terms would produce a species for which the probability of 
  #       occurrence at a site 
  #       = logit(prob.occurrence) + 0.03*elevation + 0.01*elevation^2
  #   occ.terms: an optional character vector giving the terms for the species's
  #       response to occurrence covariates, in the form of a character 
  #       string with the name of an environmental variable modified by a 
  #       polynomial term if desired.  e.g. if you have an "elevation" 
  #       variable specified in the call to create_site(), then here you 
  #       could specify c("elevation", "elevation^2) to produce a species 
  #       with a polynomial response to elevation
  #   det.coefficients: an optional vector of coefficients for each det.term 
  #       to include in the species's detection response function. This 
  #       coefficient will be on the logit scale, so that this coefficient 
  #       value can be directly compared to the parameter estimate from a 
  #       logistic regression model.
  #   det.terms: an optional character vector giving the terms for the species's
  #       response to detection covariates, in the form of a character 
  #       string with the name of a detection variable modified by a 
  #       polynomial term if desired.
  #       Note that detection covariates can include site-level covariates 
  #       defined in create_site() and observer-level covariates defined in 
  #       create_observer() (such as a variable indicating how "difficult" the
  #       species is to detect e.g. by people of lower skill level).
  # OUTPUT: a list containing the following elements:
  #   species.name: character string with the name of the species
  #   base.prob.occurrence: the baseline probability of occurrence
  #   base.prob.deteciton: the baseline probability of detection
  #   occ.equation: character string describing the response of occurrence probability to the environmental covariates
  #   det.equation: character string describing response of detection probability to detection covariates
  #   occ.resp.function: a function giving the species's detection probability 
  #       as a function of the baseline prob.occurrence (intercept) plus the 
  #       coefficients and site-level values for occurrence covariates given 
  #       in occ.terms
  #   det.resp.function: a function giving the species's detection probability 
  #       as a function of the baseline prob.detection (intercept) plus the
  #       coefficients and site/visit/observer-level values for detection 
  #       covariates given in det.terms
  #   
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
  if(!is.null(occ.coefficients) && !is.numeric(occ.coefficients)) {
    stop("Any values passed to occ.coefficients must be in the form of a numeric vector")
  }
  if(!is.null(occ.terms) && !is.character(occ.terms)) {
    stop("Any values passed to occ.terms must be in the form of a character vector")
  }
  if(!is.null(det.coefficients) && !is.numeric(det.coefficients)) {
    stop("Any values passed to det.coefficients must be in the form of a numeric vector")
  }
  if(!is.null(det.terms) && !is.character(det.terms)) {
    stop("Any values passed to det.terms must be in the form of a character vector")
  }
  
  # construct function to determine species probability of occurrence
  occ_f <- function(prob = NULL, 
                     occ.coefficients = NULL, 
                     occ.term.vals = NULL) {
    ## TODO: this will require a change in set_species_presence
    # returns the probability of occurrence adjusted by environmental covariates
    # this takes as input the values of all coefficients and all variables
    # This function is called by 'set_species_presence()'
    # ARGS: prob: the baseline probability of occurrence
    #       occ.coefficients: a list of numeric coefficients for occupancy covs
    #       occ.term.vals = a list of numeric values for environmental covariates
    #           at a single site
    logit_prob <- logit(prob)
    for (t in 1:length(occ.coefficients)) {
      coef <- occ.coefficients[[t]]
      varValue <- occ.term.vals[[t]]
      logit_prob <- logit_prob + coef * varValue
    }
    prob <- logistic(logit_prob)
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
  
  # create character strings describing the equations for occ. and det. prob.
  occ_eq <- mapply(FUN = paste, occ.coefficients, occ.terms, 
                   MoreArgs = list(sep = ""))
  occ_eq <- paste(occ_eq, collapse = " + ")
  occ_eq <- paste0("prob.occurrence = base.prob + ", occ_eq)
  det_eq <- mapply(FUN = paste, det.coefficients, det.terms,
                   MoreArgs = list(sep = ""))
  det_eq <- paste(det_eq, collapse = " + ")
  det_eq <- paste0("prob.detection = base.prob + ", det_eq)
  
  l <- list(species.name = as.character(species.name), 
            base.prob.occurrence = as.numeric(prob.occurrence), 
            base.prob.detection = as.numeric(prob.detection), 
            occ.equation = occ_eq, 
            det.equation = det_eq, 
            occ.function = occ_f, 
            det.function = det_f)
  l # return the list of characteristics for this species
}