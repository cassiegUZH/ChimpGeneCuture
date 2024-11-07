################################################################################
# BRMS models
################################################################################

# Code for Bayesian regression models to test the effects of genetic data on sharing behaviors
# using "stan" and "brms" packages

rm(list=ls())
library(dplyr)


################################################################################
## GET DYADIC DATA -----

  m.dyads <- read.delim("Data/m.dyads.txt")
  m.dyads <- m.dyads[order(m.dyads$dyad_ID, m.dyads$behavior),]
  
  str(m.dyads)
  
  m.dyads <- filter(m.dyads, !is.na(shared))

# PREPARE VARIABLES: LOG/SQRT CONTINUOUS VARIABLES AND STANDARDISE
  m.dyads$z.log.nepra <- scale(log(m.dyads$nepra_proportion))
  m.dyads$z.log.distance <- scale(log(m.dyads$distance))
  m.dyads$z.sqrt.refugia <- scale(sqrt(m.dyads$mean_refugia_distance))
  m.dyads$z.log.precipitation <- scale(log(m.dyads$mean_annual_precipitation))
  m.dyads$z.obstime <- scale(m.dyads$min_observation_time)
  
# FACTOR BEHAVIOR
  m.dyads$behavior <- as.factor(m.dyads$behavior)
  
# SUBSET DATA TO BEHAVIOURAL CLASSES
  d.nontool <- droplevels(filter(m.dyads, behavioral_class == "non-tool"))
  d.simple <- droplevels(filter(m.dyads, behavioral_class == "simple"))
  d.complex <- droplevels(filter(m.dyads, behavioral_class == "complex"))
  
  
################################################################################
## MODEL FITTING -----
  
# LIBRARIES
  library(rstan)
  library(brms)
  
# PRIORS (WEAK)
  priors_weak <- c(prior(normal(0,1), class=Intercept),
              prior(normal(0,1), class=b),
              prior(exponential(1), class=sd))
  
# FUNCTION TO RUN BAYESIAN MODEL WITH 2,000 ITERATIONS AND 1,000 WARMUP AND 4 CHAINS
  get_brms <- function(formel, dataset, priors, aD){
    bm <- brm(formula = formel, data = dataset,
              family = "bernoulli", prior = priors, 
              warmup = 1000, iter = 2000, chains = 4,
              control=list(adapt_delta =aD))} # Increase this from default = 0.95 to 0.97/0.99 if divergent transitions or Rhat > 1.01
                                              # Only the case for some non-tool models

# RUN MODELS WITH NePRA proportion AS GENETIC PREDICTOR AND SAME HABITAT AS ECOLOGICAL PREDICTOR
  
  # model formula 
  f_nepra_prop <- bf(shared ~ z.log.nepra + same_habitat + z.obstime + (1|dyad_ID) + (1|behavior))
  
  # model with complex behaviors
  m_complex_nepraprop_habitat <- get_brms(formel = f_nepra_prop, dataset = d.complex, 
                                    priors = priors_weak, aD = 0.95)
  
  # model with simple behaviors
  m_simple_nepraprop_habitat <- get_brms(formel = f_nepra_prop, dataset = d.simple, 
                                    priors = priors_weak, aD = 0.95)
  
  # model with nontool behaviors
  m_nontool_nepraprop_habitat <- get_brms(formel = f_nepra_prop, dataset = d.nontool, 
                                   priors = priors_weak, aD = 0.99)
  
  
  # This can be repeated for all combinations of genetic/geographic predictors and ecological predictors
  # note that no genetic and geographic variables are in the same model due to high correlation
  
  # example model formulas:
  f_nepra_b <- bf(shared ~ nepra_binary + z.sqrt.refugia + z.obstime + (1|dyad_ID) + (1|behavior))
  f_distance <- bf(shared ~ z.log.distance + z.log.precipitation + z.obstime + (1|dyad_ID) + (1|behavior))
  # etc.
  
  # note that for models including IBDs, rows including Outamba-Kilimi (No IBD sampling)
  # need to be excluded
  # example
  d.complex.ibd <- d.complex[!is.na(d.complex$ibd_binary),]
  f_ibd <- bf(shared ~ ibd_binary + z.sqrt.refugia + z.obstime + (1|dyad_ID) + (1|behavior))
  
  m_complex_ibd_refugia <- get_brms(formel = f_ibd, dataset = d.complex.ibd, 
                                         priors = priors_weak, aD = 0.95)
  
  
# MODELS TESTING RANDOM SLOPE EFFECTS OF BEHAVIOR
  
  # adapt model formula:
  f_nepra_prop_rs <- bf(shared ~ z.log.nepra + same_habitat + z.obstime + (1|dyad_ID) + (1 + z.log.nepra + same_habitat + z.obstime|behavior))
  
  # model example
  m_complex_nepraprop_habitat_rs <- get_brms(formel = f_nepra_prop_rs, dataset = d.complex, 
                                    priors = priors_weak, aD = 0.95)
  
  # can be repeated for all possible model formulas and behavioural classes
  
# TEST LOGARITHM OF MINIMUM OBSERVATION TIME
  
  # log of obstime variable
  d.complex$z.log.obstime <- scale(log(d.complex$min_observation_time_months))

  # example of formula:
  f_nepra_prop_obstime <- bf(shared ~ z.log.nepra + same_habitat + z.log.obstime + (1|dyad_ID) + (1|behavior))
  
  # model example:
  m_complex_nepraprop_habitat_obstime <- get_brms(formel = f_nepra_prop_obstime, dataset = d.complex, 
                                             priors = priors_weak, aD = 0.95)
  
  
  # again, this is repeated for each model combination
  
# TEST LESS INFORMATIVE WIDE PRIOR
  
  # wide prior
  priors_wide <- c(prior(normal(0,10), class=Intercept),
                   prior(normal(0,10), class=b),
                   prior(exponential(10), class=sd))
  
  
  # formula same as in original model
  f_nepra_prop <- bf(shared ~ z.log.nepra + same_habitat + z.obstime + (1|dyad_ID) + (1|behavior))
  
  # model
  m_complex_nepraprop_habitat_wide <- get_brms(formel = f_nepra_prop, dataset = d.complex, 
                                                  priors = priors_wide, aD = 0.95)
  
  
  # repeat for each model formula
  


