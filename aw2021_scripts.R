########################################################
##############   aw2021_sims_scripts   #################
# This document provides code for fitting trait divergence models using a grid-based approximate likelihood approach as in Anderson and Weir 2021 (submitted)
# The document is not meant to be run as a stand-alone script
# Instead, it provides code for each part of the analysis that can be worked into scripts by the user
# We present the code like this because running the analysis on the full suite of parameters is computationally intensive
# The use of a cluster or supercomputer will typically be required, and the correct scripting approach will vary based on the cluster
# Full scripts for running analyses on the Niagara or Cedar supercomputers in the Compute Canada consortium can be requested from seanas.anderson@mail.utoronto.ca

# ======================================================
###################  CONTENTS  #########################
# ======================================================

# Models
# ======================================================
# Attached is code for running the following tested in Anderson and Weir 2021:
# CD Model with Selection
# CD Model with no Selection
# CD Model with Selection in Sympatric Phase Only
# SS model with Selection
# SS model with no Selection
# SS model with no Externally-Imposed Wait Time to Secondary Contact
# AD model with Selection
# AD model with no Selection
# ======================================================

# Sections
# ======================================================
# The script is divided into the following sections, which represent steps in the analysis: 
# Defining Parameter Grids
# Model Fitting
# Summarizing Results 
# ======================================================

# Preliminary Setup
# ======================================================
# Load packages, data, and custom functions
require(diverge)
source("aw2021_functions.R")
billdata = read.csv("aw2021_billdata.csv")
# ======================================================


# ======================================================
######## SECTION 1: DEFINING PARAMETER GRIDS ###########
# ======================================================
# NOTE: ACTUAL PARAMETER VALUES TO BE SET BY USER
# Values here shown as an example only (most of these values were used in first exploratory round of model fits in the paper)
# For full list of values used, see aw2021_parameters_tested.R

## ========================
## CD Model With Selection ##
# In this model, paired lineages evolve by selection toward either a shared optimum or separate adaptive optima
# The distance between optima (i.e. psi, which can be zero) can change when pairs transition from allopatry to sympatry (if it doesn't change, no character displacement has occurred)
# Both the allopatric and sympatric psi values can vary with latitude
# Parameters are alpha, sig2, psi_allopatric_slope, psi_allopatric_intercept, psi_sympatric_slope, psi_sympatric_intercept, and ts_var (variance of the time-to-sympatry in units of my)
# note: ts_var can be replaced with lambda in the case of an exponential distribution for wait times to sympatry

# Define starting parameters
alpha=c(0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5))
sig2=seq(0.001, 0.009, 0.001)
ts_var=c(1,2)
rates=as.matrix(expand.grid(alpha, sig2, ts_var)) # note the 'as.matrix' part here is important - speeds it up dramatically for some reason
psi_allo_eq=c(0, c(0.001, 0.005), seq(0.01, 0.09, 0.02), seq(0.1, 0.5, 0.1)) 
psi_allo_polar=psi_allo_eq
psi_allos=model_generator(par_low=psi_allo_eq, par_high=psi_allo_polar, domain=c(0,60))
psi_sym_eq=c(0, c(0.001, 0.005), seq(0.01, 0.09, 0.02), seq(0.1, 0.5, 0.1)) 
psi_sym_polar=psi_sym_eq
psi_syms=model_generator(par_low=psi_sym_eq, par_high=psi_sym_polar, domain=c(0,60))

# Create grid of parameter combinations
# Note, we have to use nested loops here instead of expand.grid because not all psi slopes can be joined with all psi intercepts
# Slopes and intercepts are combined such that the psi parameter stays within a reasonable range across the domain of 0-60 degrees latitude
# (i.e. psi remains nonzero across the domain)
chardis_pars_selec <- matrix(nrow=nrow(rates)*nrow(psi_allos)*nrow(psi_syms), ncol=7)
colnames(chardis_pars_selec) <- c("psi_allo_slope", "psi_allo_int", "alpha", "sig2", "psi_sym_slope", "psi_sym_int", "ts_var")
counter=0
for(i in 1:nrow(psi_allos)) {
  psi_allo=psi_allos[i,]
    for(k in 1:nrow(rates)) {
      rate=rates[k,]
      for(j in 1:nrow(psi_syms)) {
        psi_sym=psi_syms[j,]
            counter=counter+1
            chardis_pars_selec[counter,1:2] <- psi_allo[1:2]
            chardis_pars_selec[counter,3:4] <- rate[1:2]
            chardis_pars_selec[counter,5:6] <- psi_sym[1:2]
            chardis_pars_selec[counter,7] <- rate[3]
    }
  }
}

## ========================
## CD Model With No Selection ##
# In this model, paired lineages evolve under a random walk process in both allopatric and sympatric phases
# The BM dispersion parameter can change between sympatry and allopatry (if it doesn't change, no character displacement has occurred)
# Both the allopatric and sympatric BM dispersion parameter can also vary with latitude
# Parameters are sig2_alpha_slope, sig2_alpha_int, sig2_sym_slope, sig2_sym_int, and ts_var (variance of the time-to-sympatry in units of my)
# note: ts_var can be replaced with lambda in the case of an exponential distribution for wait times to sympatry

# Define starting parameters
sig2_eq=c(0.001, seq(0.005, 0.009,0.001), 0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5)) 
sig2_polar = sig2_eq
sig2allos=model_generator(par_low=sig2_eq, par_high=sig2_polar, domain=c(0,60))
sig2syms = sig2allos
ts_var=seq(0.5, 3, 0.5)

# Create grid of parameter combinations
chardis_pars_rw <- matrix(nrow=nrow(sig2allos)*nrow(sig2syms)*length(ts_var), ncol=5)
colnames(chardis_pars_rw) <- c("sig2allo_slope", "sig2allo_int", "sig2sym_slope", "sig2sym_int", "ts_var")
counter=0
for(i in 1:nrow(sig2allos)) {
  sigallos=sig2allos[i,]
    for(j in 1:nrow(sig2syms)) {
      sigsyms=sig2syms[j,]
        for(m in 1:length(ts_var)) {
          ts=ts_var[m]
          counter=counter+1
          chardis_pars_rw[counter,1:2] <- sigallos
          chardis_pars_rw[counter,3:4] <- sigsyms
          chardis_pars_rw[counter,5] <- ts
        }
    }
}

## ========================
## CD Model With Selection in Sympatric Phase Only ##
# In this model, paired lineages evolve under a random walk process in allopatric phase but are driven apart by divergent selection (i.e. selection toward separate optima) in the sympatric phase
# The psi parameter in sympatry can vary with latitude
# Parameters are alpha, sig2, psi_slope, psi_int, and ts_var (variance of the time-to-sympatry in units of my)
# note: ts_var can be replaced with lambda in the case of an exponential distribution for wait times to sympatry

# Define starting parameters
alpha=c(0.001, 0.01, 0.05, seq(0.1, 0.9, 0.1), seq(1, 5, 0.5))
sig2=c(0.001, seq(0.005, 0.009,0.001), 0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5)) #july18
rates=as.matrix(expand.grid(alpha,sig2))
psi_eq= c(0, c(0.0025, 0.005, 0.0075), c(0.025, 0.05, 0.075), 0.1, seq(0.2,1,0.2), seq(1,10,1))
psi_polar=psi_eq
psis=model_generator(par_low=psi_eq, par_high=psi_polar, domain=c(0,60))
ts_var=seq(0.5, 3, 0.5)

# Create grid of parameter combinations
chardis_pars_rwallo <- matrix(nrow=nrow(rates)*nrow(psis)*length(ts_var), ncol=5)
colnames(chardis_pars_rwallo) <- c("alpha", "sig2", "psi_slope", "psi_int", "ts_var")
counter=0
for(i in 1:nrow(rates)) {
  rate=rates[i,]
    for(j in 1:nrow(psis)) {
      psipars=psis[j,]
        for(m in 1:length(ts_var)) {
          ts=ts_var[m]
          counter=counter+1
          chardis_pars_rwallo[counter,1:2] <- rate
          chardis_pars_rwallo[counter,3:4] <- psipars
          chardis_pars_rwallo[counter,5] <- ts
        }
    }
}

## ========================
## SS Model With Selection ##
# In this model, paired lineages evolve by selection toward either a shared optimum or separate adaptive optima
# The distance between optima (i.e. the psi parameter, which may be zero in the case of a single shared optimum) can vary with latitude
# After a period of geographic isolation (tc), pairs have the potential to become sympatric, and the probability of establishing sympatry is proportional to trait differences
# Parameters are psi_slope, psi_int, alpha, sig2, tc_var (variance of the time-to-contact in units of my), and inflection (the value of trait divergence at which the probability of becoming sympatric is 0.5)
# note: ts_var can be replaced with lambda in the case of an exponential distribution for wait times to sympatry

# Define starting parameters
alpha=c(0, 0.001, 0.01, 0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5), 10)
sig2=c(0.001, seq(0.005,0.009,0.001), 0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5))
tc_var = c(1,2)
inflection = seq(0.1,0.6,0.1)
psi_eq = c(0, 0.001, 0.005, seq(0.01, 0.09, 0.02), seq(.1,.5,.1),seq(1,9))
psi_polar = c(0, 0.001, 0.005, seq(0.01, 0.09, 0.02), seq(.1,.5,.1),seq(1,9))
pars = as.matrix(expand.grid(alpha, sig2, tc_var, inflection)) # note the 'as.matrix' part here is key - speeds it up dramatically
psi=model_generator(par_low=psi_eq, par_high=psi_polar, domain=c(0,60))

# Create grid of parameter combinations
# Again, we have to use nested loops here instead of expand.grid because not all psi slopes can be joined with all psi intercepts
sorting_pars_selec <- matrix(NA, nrow=nrow(pars)*nrow(psi), ncol=6)
colnames(sorting_pars_selec) <- c("psi_slope", "psi_int", "alpha", "sig2", "tc_var", "inflection")
counter=0
for(k in 1:nrow(pars)) {
  params=pars[k,]
  for(i in 1:nrow(psi)) {
    psis=psi[i,]
      counter=counter+1
      sorting_pars_selec[counter,1:2] <- psis
      sorting_pars_selec[counter,3:5] <- params[1:3]
      sorting_pars_selec[counter,6] <- params[4]
    }
}

## ========================
## SS Model With No Selection ##
# In this model, paired lineages evolve under a random walk process
# The BM dispersion parameter (sig2) can vary with latitude
# After a period of geographic isolation (tc), pairs have the potential to become sympatric, and the probability of establishing sympatry is proportional to trait differences
# Parameters are sig2_slope, sig2_int, ts_var (variance of the time-to-contact in units of my), and inflection (the value of trait divergence at which the probability of becoming sympatric is 0.5)
# note: ts_var can be replaced with lambda in the case of an exponential distribution for wait times to sympatry

# Define starting parameters
sig2_eq=c(0.001, seq(0.005, 0.009,0.001), 0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5)) 
sig2_polar = sig2_eq
sig2=model_generator(par_low=sig2_eq, par_high=sig2_polar, domain=c(0,60))
ts_var=c(seq(0.5, 7.5, 0.5),7.75)
inflection = seq(0.1,0.6,0.1)

# Create grid of parameter combinations
sorting_pars_rw <- matrix(nrow=nrow(sig2)*length(ts_var)*length(inflection), ncol=4)
colnames(sorting_pars_rw) <- c("sig2_slope", "sig2_int", "ts_var", "inflection")
counter=0
for(i in 1:nrow(sig2)) {
  sigs=sig2[i,]
    for(j in 1:length(ts_var)) {
      ts=ts_var[j]
        for(m in 1:length(inflection)) {
          infl=inflection[m]
          counter=counter+1
          sorting_pars_rw[counter,1:2] <- sigs
          sorting_pars_rw[counter,3] <- ts
          sorting_pars_rw[counter,4] <- infl
        }
    }
}

## ========================
## SS Model With no Externally-Imposed Wait Time to Secondary Contact ##
# In this model, paired lineages evolve by selection toward either a shared optimum or separate adaptive optima
# The distance between optima (i.e. the psi parameter, which may be zero in the case of a single shared optimum) can vary with latitude
# Unlike the other SS Models, this model does not impose a wait time after which pairs have the potential to become sympatric
# Pairs can become sympatric at any time in this model -- the probability of sympatry is dictated solely by trait differences
# Parameters are psi_slope, psi_int, alpha, sig2, and inflection (the value of trait divergence at which the probability of becoming sympatric is 0.5)

# Define starting parameters
alpha=c(0.001, 0.01, 0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5), 10)
sig2=c(0.001, seq(0.005, 0.009,0.001), 0.05, seq(0.1, 0.9, 0.2), seq(1, 2.5, 0.5))
inflection = seq(0.1, 2, 0.2)
pars = as.matrix(expand.grid(alpha, sig2, inflection))
psi_int = c(0, c(0.001, 0.005), seq(0.01, 0.09, 0.02), seq(0.1, 0.5, 0.1), 1, 5) 
psi_pol = psi_int
psi=model_generator(par_low=psi_int, par_high=psi_pol, domain=c(0,60))

# Create grid of parameter combinations
sorting_pars_notc <- matrix(NA, nrow=nrow(pars)*nrow(psi), ncol=5)
colnames(sorting_pars_notc) <- c("psi_slope", "psi_int", "alpha", "sig2", "inflection")
counter=0
for(k in 1:nrow(pars)) {
  params=pars[k,]
  for(i in 1:nrow(psi)) {
    psis=psi[i,]
      counter=counter+1
      sorting_pars_notc[counter,1:2] <- psis
      sorting_pars_notc[counter,3:5] <- params
    }
}

## ========================
## AD Model With Selection  ##
# In this model, paired lineages evolve by selection toward either a shared optimum or separate adaptive optima
# Note that CD models of selection in which psi_allo == psi_sym collapse into this model
# The distance between optima (i.e. the psi parameter, which may be zero in the case of a single shared optimum) can vary with latitude
# Unlike the SS models, all pairs have equal probability of establishing sympatry in the AD model
# Parameters are psi_slope, psi_int, alpha, and sig2

# Define starting parameters
psi_int = c(c(0, 0.001, 0.005), seq(0.01,0.09,0.02), seq(0.1, 0.5, 0.1), seq(1, 9, 1))
psi_pol = psi_int
psi=model_generator(par_low=psi_int, par_high=psi_pol, domain=c(0,60))
alpha = c(seq(0.002, 0.009, 0.001), seq(0.02, 0.09, 0.01))
sig2 = c(0.001, seq(0.005, 0.009, 0.001), 0.050, 0.100)
pars=as.matrix(expand.grid(alpha,sig2))

# Create grid of parameter combinations
agediff_pars_selec <- matrix(NA, nrow=nrow(pars)*nrow(psi), ncol=4)
colnames(agediff_pars_selec) <- c("psi_slope", "psi_int", "alpha", "sig2")
counter=0
for(k in 1:nrow(pars)) {
  params=pars[k,]
  for(i in 1:nrow(psi)) {
    psis=psi[i,]
    counter=counter+1
    agediff_pars_selec[counter,1:2] <- psis
    agediff_pars_selec[counter,3:4] <- params
  }
}

## ========================
## AD Model with No Selection ##
# In this model, paired lineages evolve under a random walk process
# The BM dispersion parameter (sig2) can vary with latitude
# Note that CD models with no selection in which psi_allo == psi_sym collapse into this model
# Parameters are sig2_slope, sig2_intercept

# Define starting parameters
sig2 = c(seq(0.0001, 0.0009, 0.0001), seq(0.001, 0.009, 0.001), seq(0.01, 1, 0.01), seq(1,10,1), 100)

# Create grid of parameter combinations
agediff_pars_rw = model_generator(domain=c(0,60), par_low = sig2, par_high = sig2)

# END OF SECTION 1
# ======================================================




# ======================================================
############  SECTION 2: MODEL FITTING  ################
# ======================================================
# NOTE 1: the code below uses the function mclapply, which uses a forking routine and is not available on windows
# It can be substitued by parLapply, which requires the creation of a socket cluster
# Note 2: most model-fitting functions in aw2021_functions.R have default arguments and won't work unless the 'billdata' dataset is in your environment
# While this is annoying, it speeds things up when you're running the function hundreds of thousands or millions of times.
# NOTE 3: obvious point: for testing and debugging, use small parameter grids!!
# NOTE 4: In the code below, we fit models using small subsets (the first ten rows) of the grids we defined in section 1

# Define no. simulations to run for each parameter combination (SET BY USER)
# NOTE: the function will not work if you set N to a value too low
# This is because the function fits density models to a simulated distribution, but we don't allow it to set density models if that distribution is comprised of too-few values, as this makes the resulting density function is less reliable.
# Thus, for testing purposes, user few parameter combinations instead of small N values.
N=5000

# Define number of cores available on the machine for parallel runs of the analysis (SET BY USER)
ncor = 2

# Separate allopatric and sympatric data vectors (req'd as default arguments when fitting CD and SS models)
allo_lat=abs(billdata$lat[billdata$patry==0]) 
allo_age=billdata$age[billdata$patry==0]
allo_ed=billdata$ed[billdata$patry==0]
sym_lat=abs(billdata$lat[billdata$patry==1])
sym_age=billdata$age[billdata$patry==1]
sym_ed=billdata$ed[billdata$patry==1]
# define the mean times to contact for allopatric pairs (based on Weir and Price (2011), req'd for SS models only)
meantc_allo=-0.0238*allo_lat + 3.18

## ========================
## Character Displacement Models ##
# 1. Fitting CD Model with Selection
# Convert parameter matrix to list (it is MUCH faster to run mclapply on a list than it is to run parallel operators on matrices)
cdpars=as.list(as.data.frame(t(chardis_pars_selec[1:10,])),all.names=TRUE) 

# Run the simlikely_cd function on each parameter combo
# likely_cdselec = simlikely_cd(pars=cdpars[[1]], N=N, exponential=FALSE) # single parameter combo, set exponential=TRUE in any of these codelines to model Tsym under that distribution
# likely_cdselec = lapply(cdpars, FUN=simlikely_cd, N=N) # multiple parameter combos in serial
likely_cdselec = mclapply(cdpars, FUN=simlikely_cd, N=N, mc.cores=ncor) # multiple parameter combos in parallel

# 2. Fitting CD Model with No Selection
# Convert parameter matrix to list
cdpars_rw=as.list(as.data.frame(t(chardis_pars_rw[1:10,]))) 

# Run the simlikely_cdbm function on each parameter combo
# likely_cdbm = simlikely_cdbm(pars=cdpars_rw[[1]], N=N, exponential=FALSE)
# likely_cdbm = lapply(cdpars_rw, FUN=simlikely_cdbm, N=N) 
likely_cdbm = mclapply(cdpars_rw, FUN=simlikely_cdbm, N=N, mc.cores=ncor) 

# 3. Fitting CD Model with Selection in the Sympatric Phase Only
# Convert parameter matrix to list
cdpars_rwallo=as.list(as.data.frame(t(chardis_pars_rwallo[1:10,]))) 

# Run the simlikely_cdbmallo function on each parameter combo
# likely_cdbmallo = simlikely_cdbmallo(pars=cdpars_rwallo[[1]], N=N, exponential=FALSE)
# likely_cdbmallo = lapply(cdpars_rwallo, FUN=simlikely_cdbmallo, N=N)
likely_cdbmallo = mclapply(cdpars_rwallo, FUN=simlikely_cdbmallo, N=N, mc.cores=ncor)


## ========================
## Species Sorting Models ##
# 1. Fitting SS Model with Selection
# Convert parameter matrix to list
ss_pars=as.list(as.data.frame(t(sorting_pars_selec[1:10,])))

# Run the simlikely_ss function on each parameter combo
#likely_ss_selec = simlikely_ss(pars=ss_pars[[1]], N=N, exponential=FALSE) # single parameter combo; set exponential = TRUE in the simlikely_ss function to model times to secondary contact from that distribution
#likely_ss_selec = lapply(ss_pars, FUN=simlikely_ss, N=N) # multiple parameter combos in serial
likely_ss_selec = mclapply(ss_pars, FUN=simlikely_ss, N=N, mc.cores=ncor) # multiple parameter combos in parallel

# 2. Fitting SS Model with no Selection
# Convert parameter matrix to list
ss_pars_rw=as.list(as.data.frame(t(sorting_pars_rw[1:10,])))

# Run the simlikely_ss function on each parameter combo
#likely_ssbm = simlikely_ss(pars=ss_pars_rw[[1]], N=N)
#likely_ssbm = lapply(ss_pars_rw, FUN=simlikely_ss, N=N)
likely_ssbm = mclapply(ss_pars_rw, FUN=simlikely_ss, N=N, mc.cores=ncor)

# 3. Fitting the SS Model with no Externally-Imposed Wait Time to Secondary Contact
# Convert parameter matrix to list
ss_pars_notc=as.list(as.data.frame(t(sorting_pars_notc[1:10,])))

# Run the simlikely_ss_notc function on each parameter combo
#likely_ss_notc = simlikely_ss_notc(pars=ss_pars_notc[[1]], N=N)
#likely_ss_notc = lapply(ss_pars_notc, FUN=simlikely_ss_notc, N=N)
likely_ss_notc = mclapply(ss_pars_notc, FUN=simlikely_ss_notc, N=N, mc.cores=ncor)


## ========================
## Age Difference Models ##
# 1. Fitting AD Model with Selection
# Convert parameter matrix to list
ad_pars=as.list(as.data.frame(t(agediff_pars_selec[1:10,])))

# Run the simlikely_ad function on each parameter combo
#likely_ad_selec = simlikely_ad(pars=ad_pars[[1]], N=N) # single parameter combo
#likely_ad_selec = lapply(ad_pars, FUN=simlikely_ad, N=N) # multiple parameter combos in serial
likely_ad_selec = mclapply(ad_pars, FUN=simlikely_ad, N=N, mc.cores=ncor) # multiple parameter combos in parallel

# 2. Fitting AD Model with no Selection
# Convert parameter matrix to list
ad_pars_rw=as.list(as.data.frame(t(agediff_pars_rw[1:10,])))

# Run the simlikely_ad function on each parameter combo
#likely_adbm = simlikely_ad(pars=ad_pars_rw[[1]], N=N)
#likely_adbm = lapply(ad_pars_rw, FUN=simlikely_ad, N=N)
likely_adbm = mclapply(ad_pars_rw, FUN=simlikely_ad, N=N, mc.cores=ncor)

# END OF SECTION 2
# ======================================================


# ======================================================
##########  SECTION 3: SUMMARIZING RESULTS  ############
# ======================================================
# To summarize results we use the 'summarizer' function, which calculates AICc for each paramter combination and returns a matrix in descending order from best to worst fit
# The first columns of each matrix are parameter values (number depends on the model)
# The last four columns are the likelihood and summary statistics corresponding to those parameter values
# note: the 'parmat' argument in the summarizer function refers to the original paramter matrix we created in section 1, not the parameter list we created in section 2

# Testing the effect of outliers
# note that the rm.outliers argument in the summarizer function allows you to summarize model fit based on a dataset with the two most extremely differentiated pairs removed
# These taxa are:
# Nasica_longirostris group; age: 6.7721018, bill shape differentiation: 1.77398886
# Tiaris_olivaceus group; age: 5.3425478, bill shape differentiation: 1.19927600

# To run the function, first rearrange age vector to correspond with sympatric and allopatric pairs
# This is required to make column names of the raw output line up properly when using the summarizer model
allos = billdata$age[billdata$patry==0]
syms = billdata$age[billdata$patry==1]
sp.pairs = c(syms, allos)

# CD MODELS
# CD Model with Selection
res_cd_selec = summarizer(parmat=chardis_pars_selec[1:10,], reslist=likely_cdselec, rm.outliers=FALSE) # set rm.outliers=TRUE to calculate model fit when the two most extremely differentiated pairs are removed
res_cd_selec
# CD Model with no Selection
res_cd_bm = summarizer(parmat=chardis_pars_rw[1:10,], reslist=likely_cdbm)
res_cd_bm
# CD Model with Selection Only in the Sympatric Phase
res_cd_bmallo = summarizer(parmat=chardis_pars_rwallo[1:10,], reslist=likely_cdbmallo)
res_cd_bmallo

# SS MODELS
# SS Model with Selection
res_ss_selec = summarizer(parmat=sorting_pars_selec[1:10,], reslist=likely_ss_selec, rm.outliers=FALSE) # set rm.outliers=TRUE to calculate model fit when the two most extremely differentiated pairs are removed
res_ss_selec

# SS Model with No Selection
res_ss_bm = summarizer(parmat=sorting_pars_rw[1:10,], reslist=likely_ssbm)
res_ss_bm

# SS Model with No Externally-Imposted Wait Time to Secondary Contact
res_ss_notc = summarizer(parmat=sorting_pars_notc[1:10,], reslist=likely_ss_notc)
res_ss_notc

# AD MODELS
# AD Model with Selection
res_ad_selec = summarizer(parmat=agediff_pars_selec[1:10,], reslist=likely_ad_selec, rm.outliers=FALSE) # set rm.outliers=TRUE to calculate model fit when the two most extremely differentiated pairs are removed
res_ad_selec

# AD Model with no Selection
res_ad_rw = summarizer(parmat=agediff_pars_rw[1:10,], reslist=likely_adbm)
res_ad_rw

# END OF SECTION 3
# ======================================================


