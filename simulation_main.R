#########################
## Simulation for testing species recording models
##
## This script outlines the structure of the simulation code and provides 
## instructions on how to use and extend the code.
##
## author: Willson Gaul
## created: 9 Jun 2017
## last modified: 20 July 2017
##
#########################

#### main structure
## The following objects will exist:
#
#   - A df of species and their traits
#     -- species as rows, traits as columns
#   - A df of sites
#   - A df of observers
#     -- individual observers as rows, observer traits as columns
#     -- sites as rows, environmental conditions as columns
#   - A df of realized species communities
#     -- sites as rows, species as columns with 1/0 indicating species 
#       presence, not detection
#   - A df of realized samples
#     -- sites as rows, species as columns indicating detection 
#         rather than occurrence
#     -- can have multiple samples (rows) for each site

##### Need functions to:
### these first 3 functions output the rows for the first 3 dfs listed above
# use these to build the universe which creates the realized datasets

#   - create_species
#     -- pass in parameters of species responses to env. variables

#   - create_site
#     -- pass in environmental values for the site

#   - create_observer
#     -- pass in values for observer traits

#### these functions create the individual cells for realized datasets

#   - set_species_presence
#     -- creates a presence/absence value for a species 
#     -- inputs: species parameters, environmental parameters of site

#   - observe_species
#     -- creates a 1/0 observed value for a species
#     -- inputs: species presence/absence value, detection covariates, 
#         observer traits

#   - create site-level species communities
#     -- use the set_species_presence function to get a presence/absence value for
#       each species based on the species's parameters and environmental 
#       parameters at a site
#     -- inputs: a df or list of site environmental parameters, 
#          a df or list of species and their response parameters
#     -- output: df of sites as rows, species as columns, 1 or 0 for each sp.


# setwd("~/Documents/Data_Analysis/UCD/simulation/")
# setwd("~/Dropbox/Link_to_Documents/Data_Analysis/UCD/simulation/")

require(dplyr)

### Function Definitions
logit <- function(p) {
  # ARGS: p: a probability (0 to 1)
  log(p / (1 - p))
}
logistic <- function(x) {
  # ARGS: x: a linear combination of parameters (e.g. prob + coef*varValue)
  1 / (1 + exp(-x))
}

################## Begin Create Functional Units ################
###### Each of these functions creates a single functional unit
###### (e.g. a single species, a single site, a single observer)
###### These are meant to be called by the larger aggregating functions,
###### not directly by the user.  

## Now I am wondering if there should be larger aggregation functions or if
## the user should just call these directly in for loops or something.

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
  
  # l will be the object returned at the end of the function
  l <- list(species.name = species.name, 
            base.prob.occurrence = prob.occurrence, 
            base.prob.detection = prob.detection, 
            occ.function = NULL, 
            det.function = NULL)
  
  # parse occurrence coefficients and terms
  
  
  # parse detection coefficients and terms
  
  
  l # return the list of characteristics for this species
}

create_site <- function(site.name, land.class = NULL, ...) {
  # Create a single site based on user-specified environmental variables
  # ARGS: site.name: a character string with the site name
  #       land.class: an optional character string designating land use class
  #           which will automatically be converted to a dummy variable
  #         ... optional names and values for additional environmental variables
  # output: a 1-row dataframe with the site as row and environmental variables 
  #           as columns
  #
  # pass in as many environmental variables as you want as individual arguments
  # where the name of the variable is the name of the argument, and the
  # value of the variable is the value of the argument.
  # e.g. alt_site(site.name = "st1", altitude = 450, land_class = "forest")
  vars <- list(...)
  site <- as.data.frame(matrix(nrow = 1, 
                               ncol = 1 + length(vars) + length(land.class)))
  colnames(site) <- c("site.name", names(vars))
  site$site.name <- site.name
  
  # add additional environmental variables if there are any
  if (length(vars) > 0) { 
    site[1, 2:(1 + length(vars))] <- vars
  }
  
  # add dummy variable for land use class
  if (!is.null(land.class)) {
    colnames(site)[ncol(site)] <- land.class
    site[1, ncol(site)] <- 1
  }

  site # return 1-row data frame with site name and envirnmental variables
}

### create observer
# input: observer name, names and values for observer traits
# output: a list containing trait names and their values for a single observer

create_observer <- function(observer.name, ...) {
  ob <- list(name = observer.name, ...)

  ###### old data frame way
  # ob <- as.data.frame(matrix(nrow = 1, ncol = 1 + length(vars)))
  # colnames(ob) <- c("observer.name", names(vars))
  # ob$observer.name <- observer.name
  
  ob # return a list with a single observer's name and characteristics
}
######################## End Functional Unit Functions #############

###################### Begin Individual Realizing Functions ##################

set_species_presence <- function(species.row = NULL, site.row = NULL) {
  # Determine whether a species is present at a site
  # This will not be called by the user - the user will call the "create
  # community" function or something like that, which will call this for each
  # species/site cell.
  # ARGS:
  #   species.row: a 1-row data frame for the species
  #   site.row: a 1-row data frame for the site
  #   e.g. set_species_presence(0.05, resp.elevation = 0.001, resp.forest = -1)
  # OUTPUT: a presence/absence value for a species 
  if (is.null(species.row) || is.null(site.row)) {
    stop("You must pass in species and site dataframes with columns detailing
         the environmental variables at the sites and species' responses to 
         those variables.")
  }
  if (nrow(species.row) > 1 || nrow(site.row) > 1) {
    stop("There is more than one row in the species or site dataframes.")
  }
  
  # extract names of species response variables
  # strip 'resp.' prefix from species.row colum names so they match 
  # site environmental variable names
  colnames(species.row)[2:ncol(species.row)] <- gsub("resp.", "", 
                                                     names(species.row)[2:ncol(
                                                       species.row)])
  sp_response_vars <- names(species.row)
  sp_response_vars <- sp_response_vars[-which(
    sp_response_vars %in% c("species.name", 
                            "prob.occurrence", 
                            "prob.detection"))]
  
  # extract names of site environmental variables
  site_env_vars <- names(site.row)
  site_env_vars <- site_env_vars[-which(site_env_vars == "site.name")]
  
  # get baseline probability, which will be modified potentially multiple times
  # by environmental variables below
  prob <- species.row$prob.occurrence[1]
  
  # if there are species response variables, use them to modify probability of
  # occurrence for the species
  if (length(sp_response_vars) > 0) {
    # find variables for which both species response and value at the site
    # are specified.  These can be used to modify probability of occurrence
    # for the species at this site.
    variables <- sp_response_vars[which(sp_response_vars %in% site_env_vars)]
  }
  if (exists("variables") && length(variables) > 0) { # if sites and species share some variables
    for(i in 1:length(variables)) { # loop through usable variables
      # - probability is adjusted by adding the result of
      # (sp response slope * environmental variable value at site) to 
      # the logit(prob) by passing (slope*variableValue) into logistic function
      # - Check for NAs because both the species df and the site df may have
      # rows which have some environmental variables unspecified (and therefore
      # NA).  Those NAs will cause the multiplication below to give a 'prob'
      # of NA, which will make it impossible to determine species presence.
      # - NOTE: using double bracket subsetting because species.row and site.row
      # have class tbl_df which produces a 1x1 df instead of just the value
      # using single bracket subsetting.
      if (!is.na(species.row[[1, variables[i]]]) &&
          !is.na(site.row[[1, variables[i]]])) {
        coef <- species.row[[1, variables[i]]]
        varValue <- site.row[[1, variables[i]]]
        logit_prob <- logit(prob) + coef * varValue
        prob <- logistic(logit_prob)
      }
    }
  }
  
  # stop if probability is outside 0 to 1
  if (prob > 1 || prob < 0) stop("Probability is outside the 0 to 1 range in set_species_presence(). Perhaps logistic transformation didn't work?")

  # generate observation
  value <- rbinom(n = 1, size = 1, prob = prob)
  
  value # return 1 or 0 for present or absent
}

### observe_species
# inputs: species.present - a 1 or 0 indicating species presence or absence
#         species.df - a 1-row data frame containing at least the name and 
#                     probability of detection for this species.  May
#                     also contain other attributes of the species such as 
#                     identification difficulty which might interact with 
#                     observer skill to influence detection
#         observer - a list containing the name and traits of the observer
#         ... - optional detection covariates (wind, day of year, etc)
# outputs: a 1 or 0 indicating whether the species was observerd or not

observe_species <- function(species.present = NULL, species.df = NULL, 
                            observer = NULL, ...) {
  if (is.null(species.present)) {
    stop("Species presence or absence not specified.")}
  if (species.present != 1 & species.present != 0) stop(
    "species.present argument must be a 1 or 0.")
  if (is.null(species.df)) {
    stop("Data frame with species detection probability needed.")}
  if (nrow(species.df) != 1) {
    stop("species.df must be a 1-row dataframe with parameters for a single species.")
  }
  if (is.null(observer)) {
    warning("No observer traits specified.  Unmodified detection probability used.")
  }
  
  # if species is not present it will not be observed.  Right now there is no
  # option for creating "false positive" observations.
  if (species.present == 0) return(0) 
  
  det_vars <- list(...) # capture detection covariates if any

  # set baseline detection probability for this species
  det_prob <- species.df$prob.detection[1]
  
  if (!is.null(observer)) {
    # modify detection probability based on observer skill and identification
    # difficulty of the species
  }
  
  if (length(det_vars) > 0) {
    # modify detection probability based on detection covariates
  }
  
  # stop if probability is outside 0 to 1
  if (det_prob > 1 || det_prob < 0) stop("Probability is outside the 0 to 1 range in observe_species(). Perhaps logistic transformation didn't work?")
  
  # create an observerd/not observed value using the final detection probability
  value <- rbinom(size = 1, n = 1, prob = det_prob) 
  value # return 1 or 0 for detected or not detected
}


####################### end individual realizing functions ##########

######################## Begin dataset realizing functions ############

generate_community_data <- function(species.df = NULL, site.df = NULL) {
  # Create a dataframe of true species presence/absence at sites
  # This is called by the user, and calls the set_species_presence() function
  # ARGS: species.df: a data frame of species with responses to envrionmental 
  #         variables 
  #       site.df: a data frame of sites with env. variable values
  # output: a community data frame with sites as rows, species as columns, and
  #         a 1 or 0 indicating presence or absence of each species at each site
  if (is.null(species.df) || is.null(site.df)) {
    stop("You must provide data frames giving species names and responses to 
         environmental variables and site names and environmental variable values.")
  }
  # initialize results df
  # a row for each site, a column for each species, plus a column for site name
  sp_occurrences <- as.data.frame(matrix(nrow = nrow(site.df), 
                                        ncol = 1 + nrow(species.df)))
  colnames(sp_occurrences) <- c("site.name", species.df$species.name)
  # fill in site names
  sp_occurrences$site.name <- site.df$site.name
  
  for(i in 1:nrow(sp_occurrences)) {
    site <- sp_occurrences$site.name[i] # name of site to fill
    for(j in 2:ncol(sp_occurrences)) {
      sp <- colnames(sp_occurrences)[j] # name of species to fill
      
      # pass the row for this site and the row for this species into 
      # set_species_presence to get a 1/0 value
      sp_occurrences[i, sp] <- set_species_presence(
        species.row = species.df[which(species.df$species.name == sp), ], 
        site.row = site.df[site.df$site.name == site, ])
    }
  }
  sp_occurrences # return the data frame of species occurrences at sites
}

sample_site <- function(sp.occurrences = NULL, species.df = NULL, n = 1, 
                        observers = NULL) {
  # Create a sample of species at a single site based on the 
  # occurrence data produced by generate_community_data()
  #
  # ARGS: sp.occurrences: a 1-row dataframe with site name as the first column,
  #           species as the remaining columns, with 1/0 presence/absence data
  #       species.df: a dataframe defining the traits of each species, including
  #           detection probability, with species as rows and traits as columns
  #       n: the number of samples (or visits) to produce for this site
  #       observers: an optional list containing lists defining the traits 
  #           of observers
  # outputs: a community data frame with site visits as rows, a column for site
  #         name, a column for observer name, and columns for each potential 
  #         species.  Detection or non-detection of each species at each site is
  #         indicated with a 1 or 0
  #
  # The user calls this directly, and it calls observe_species()
  if (is.null(sp.occurrences)) {
    stop("You must provide species occurrence data.")
  }
  if (is.null(species.df)) {
    stop("You must provide a data frame containing at least the species names 
         and their detection probability.")
  }
  if (nrow(sp.occurrences) > 1) {
    stop("More than one row passed into the sp.occurrences argument.")
  }
  if (n < 1) {
    return()
  }
  
  # initialize results df with a row for each sample (visit) and a column
  # for each species plus a column for site name
  samples <- as.data.frame(matrix(nrow = n, ncol = ncol(sp.occurrences)))
  colnames(samples) <- colnames(sp.occurrences)
  
  samples$site.name <- rep(sp.occurrences$site.name[1], times = nrow(samples))
  
  for(i in 1:nrow(samples)) {
    for(j in 2:ncol(samples)) {
      samples[i, j] <- observe_species(
        species.present = sp.occurrences[1, j], 
        species.df = species.df[species.df$species.name == 
                                  colnames(samples)[j], ])
    }
  }
  samples # return the df of n site samples (visits) for this site
}


#######################################################################
######################################################################
##########################################################################
################# testing
