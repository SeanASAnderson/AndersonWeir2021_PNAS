# Dependencies
require(parallel)
require(locfit)
require(truncnorm)
require(VGAM)

#===========================================
# Character displacement function 1
# Fits a CD model in which trait differentiation results from divergent selection, where the magnitude of divergent selection
# is greater in sympatry than in allopatry (as expected from competition)
# Simulate trait divergence and calculate likelihood under character displacement model
#' @param sym_div: Vector of trait divergence values for sympatric pairs
#' @param allo_div: Vector of trait divergence values for allopatric pairs
#' @param sym_ages: Vector of species pair ages for sympatric pairs
#' @param allo_ages: Vector of species pair ages for allopatric pairs
#' @param sym_lats: Vector of species pair latitudes (see 'methods') for sympatric pairs
#' @param allo_lats: Vector of species pair latitudes (see 'methods') for allopatric pairs
#' @param pars: Vector of model parameter values (see second note below)
#' @param N: number of replicates to simulate for each age in 'ages'
#' @return A vector containing the log-likelihood of models parameters given each pair, plus the value of the sum of the log-likelihoods
# Vector returned of length == length(div) + 1
### NOTE: pars vector must be in the order: c(psi_allo_slope, psi_allo_intercept, alpha, sig2, psi_sym_slope, psi_sym_intercept, ts_var)
simlikely_cd = function(pars, N, sym_div=sym_ed, allo_div=allo_ed, allo_ages=allo_age, 
  allo_lats=allo_lat, sym_ages=sym_age, sym_lats=sym_lat, exponential=FALSE) {
  
  if(length(pars)!=7) {
    stop("Double check your parameter vector; you need to specify 7 parameters to run this function")
  }
  
  # Calculate mean time to sympatry
  meanTs=-0.0238*sym_lats + 3.18 # equation of line from Weir and Price (2011)
  
  # Define/calculate parameters
  psi1_allo=allo_lat*pars[1] + pars[2]
  psi1_sym=sym_lat*pars[1] + pars[2]
  alpha=pars[3]
  sig2=pars[4]
  psi2=sym_lat*pars[5] + pars[6]
  timepar=pars[7]
  
  # calculate standard deviation of trait differentiation
  sd_allo = sqrt((sig2/alpha)*(1-exp(-2*alpha*(allo_ages)))) 
  sd_sym = sqrt((sig2/alpha)*(1-exp(-2*alpha*(sym_ages))))
  
  # draw N divergences for sympatric taxa of each age/lat combo
  simdiv_sym = vector("list", length(sym_ages))
  for(i in 1:length(sym_ages)) {
    if(exponential==FALSE) {
      Tsym = rtruncnorm(n=N, a=0, b=as.numeric(sym_ages[i]), mean=meanTs[i], sd=sqrt(timepar))
    } else {
      Tsym = rexp(n=N, rate=timepar)
    }
    u2_sym=exp(-alpha*sym_ages[i])*(psi1_sym[i]*(exp(alpha*Tsym)-1)+psi2[i]*(exp(alpha*sym_ages[i])-exp(alpha*Tsym))) # vector of N u2s
    simdiv_sym[[i]] = rtruncnorm(n = N, a=0, mean= u2_sym[i], sd = sd_sym[i])
  }

  # draw N divergences for allopatirc taxa of each age/lat combo
  simdiv_allo = vector("list", length(allo_ages)) 
  u2_allo = psi1_allo*(1-exp(-alpha*allo_ages))
  for(i in 1:length(allo_ages)) {
    simdiv_allo[[i]] = rtruncnorm(n = N, a=0, mean= u2_allo[i], sd = sd_allo[i])
  }

  # generate density functions for allo and sympatric simulated datasets
  dens_sym = suppressWarnings(lapply(simdiv_sym, function(x) {
        j = max(x)
      locfit(~lp(x), xlim=c(0,j+0.5), maxk=200)})) 
  dens_allo = suppressWarnings(lapply(simdiv_allo, function(y) locfit(~lp(y), xlim=c(0,10), maxk=200)))

  # calculate likelihood of the different pairs given the model
  loglike_sym = rep(0, length(sym_ages))
  loglike_allo = rep(0, length(allo_ages))
  for(i in 1:length(sym_ages)) {
    loglike_sym[i] = log(predict(dens_sym[[i]], sym_div[i]))
  }
  for(i in 1:length(allo_ages)) {
    loglike_allo[i] = log(predict(dens_allo[[i]], allo_div[i]))
  }
  loglikely = c(loglike_sym, loglike_allo)
  loglikely = c(loglikely, sum(loglikely))

  return(loglikely)
}

## example function call
#load("aw2020_billdata.Rdata")
#allo_ed = abs(billdata$ed[billdata$patry==0])
#allo_lat=abs(billdata$lat[billdata$patry==0]) 
#allo_age=billdata$age[billdata$patry==0]
#sym_ed = abs(billdata$ed[billdata$patry==1])
#sym_lat=abs(billdata$lat[billdata$patry==1])
#sym_age=billdata$age[billdata$patry==1]
#simpars = c(0, 0.1, 0.8, 0.005, -0.006, 3, 2)
#res = simlikely_chardis(pars=simpars, N=5000)

#===========================================
# Character displacement function 2
# Here instead of divergent selection, trait divergence evolves under a Brownian motion process
# The intensity of random dispersion differs between allopatry and sympatry
# The intensity of random dispersion varies with latitude in both cases (allo and sympatric)
# params are the same as described in simlikely_cd
### NOTE: pars vector must be in the order: c(sig2_allo_slope, sig2_allo_intercept, sig2_sym_slope, sig2_sym_intercept, ts_var)
simlikely_cdbm = function(pars, N, sym_div=sym_ed, allo_div=allo_ed, allo_ages=allo_age, 
  allo_lats=allo_lat, sym_ages=sym_age, sym_lats=sym_lat, exponential=FALSE) {
  
  if(length(pars)!=5) {
    stop("Double check your parameter vector; you need 5 parameters to run this function")
  }
  
  # calc mean time to sympatry
  meanTs= -0.0238*sym_lat + 3.18 # equation of line from Weir and Price (2011)
  
  # Define/calculate parameters
  sig2allo=allo_lat*pars[1] + pars[2]
  sig2allosym = sym_lat*pars[1] + pars[2]
  sig2sym=sym_lat*pars[3] + pars[4]
  timepar=pars[5]

  # calculate standard deviations and simulate trait differentiation
  sd_allo = sqrt(2*sig2allo*allo_ages)
  simdiv_allo = lapply(sd_allo, rtruncnorm, mean=0, n=N, a=0)
  simdiv_sym = vector("list", length(sym_ages))
  for(i in 1:length(sym_ages)) {
    if(exponential==FALSE) {
      Tsym = rtruncnorm(n=N, a=0, b=as.numeric(sym_ages[i]), mean=meanTs[i], sd=sqrt(timepar))
    } else {
      Tsym = rexp(n=N, rate=timepar)
    }
    sd_sym = sapply(Tsym, function(x) sqrt(2*sig2allosym[i]*x + 2*sig2sym[i]*(sym_ages[i]-x)))
    simdiv_sym[[i]] = rtruncnorm(n = N, a=0, mean= 0, sd = sd_sym)
  }

  # generate density functions for allo and sympatric simulated datasets
  dens_allo = suppressWarnings(lapply(simdiv_allo, function(x) locfit(~lp(x), xlim=c(0,10), maxk=200)))
  dens_sym = suppressWarnings(lapply(simdiv_sym, function(x) locfit(~lp(x), xlim=c(0,10), maxk=200)))

  # calculate likelihood of the different pairs given the model
  loglike_sym = rep(0, length(sym_ages))
  for(i in 1:length(sym_ages)) loglike_sym[i] = log(predict(dens_sym[[i]], sym_div[i]))
  loglike_allo = rep(0, length(allo_ages))
  for(i in 1:length(allo_ages)) loglike_allo[i] = log(predict(dens_allo[[i]], allo_div[i]))
  loglikely = c(loglike_sym, loglike_allo, sum(loglike_sym, loglike_allo))
  return(loglikely)
}

## example function call
#load("aw2020_billdata.Rdata")
#allo_ed = abs(billdata$ed[billdata$patry==0])
#allo_lat=abs(billdata$lat[billdata$patry==0]) 
#allo_age=billdata$age[billdata$patry==0]
#sym_ed = abs(billdata$ed[billdata$patry==1])
#sym_lat=abs(billdata$lat[billdata$patry==1])
#sym_age=billdata$age[billdata$patry==1]
#simpars = c(0, 0.005, 0, 0.01, 2)
#res = simlikely_chardis_bm(pars=simpars, N=5000)


#===========================================
# Character displacement function 3
# Trait differentiation results from Brownian motion processes in allopatric stage only
# Upon sympatry, divergent selection kicks in
# params are the same as described in simlikely_cd
### NOTE: pars vector must be in the order: c(alpha, sig2, psi_sym_slope, psi_sym_intercept, ts_var)
simlikely_cdbmallo = function(pars, N, sym_div=sym_ed, allo_div=allo_ed, allo_ages=allo_age, 
  allo_lats=allo_lat, sym_ages=sym_age, sym_lats=sym_lat, exponential=FALSE) {
  
  if(length(pars)!=5) {
    stop("Double check your parameter vector; you need 5 parameters to run this function")
  }
  
  # calc mean time to sympatry
  meanTs= -0.0238*sym_lat + 3.18 # equation of line from Weir and Price (2011)
  
  # Define/calculate parameters
  alpha = pars[1]
  sig2 = pars[2]
  psi = sym_lat*pars[3] + pars[4]
  timepar=pars[5]

  # calculate standard deviations and simulate trait differentiation
  sd_allo = sqrt(2*sig2*allo_ages)
  simdiv_allo = lapply(sd_allo, rtruncnorm, mean=0, n=N, a=0)
  simdiv_sym = vector("list", length(sym_ages))
  for(i in 1:length(sym_ages)) {
    if(exponential==FALSE) {
      Tsym = rtruncnorm(n=N, a=0, b=as.numeric(sym_ages[i]), mean=meanTs[i], sd=sqrt(timepar))
    } else {
      Tsym = rexp(n=N, rate=timepar)
    }
    sd_sym = sapply(Tsym, function(x) sqrt(2*sig2*x + (sig2/alpha)*(1-exp(-2*alpha*(sym_ages[i]-x)))))
    u = sapply(Tsym, function(x) (1 - exp(-alpha*(sym_ages[i]-x)))*psi[i])
    simdiv_sym[[i]] = rtruncnorm(n = N, a=0, mean=u, sd = sd_sym)
  }

  # generate density functions for allo and sympatric simulated datasets
  #dens_sym = suppressWarnings(lapply(simdiv_sym, function(x) {
        #j = max(x)
      #locfit(~lp(x), xlim=c(0,j+0.5), maxk=200)})) 
  dens_allo = suppressWarnings(lapply(simdiv_allo, function(x) locfit(~lp(x), xlim=c(0,10), maxk=200)))
  dens_sym = suppressWarnings(lapply(simdiv_sym, function(x) locfit(~lp(x), xlim=c(0,10), maxk=200)))

  # calculate likelihood of the different pairs given the model
  loglike_sym = rep(0, length(sym_ages))
  for(i in 1:length(sym_ages)) loglike_sym[i] = log(predict(dens_sym[[i]], sym_div[i]))
  loglike_allo = rep(0, length(allo_ages))
  for(i in 1:length(allo_ages)) loglike_allo[i] = log(predict(dens_allo[[i]], allo_div[i]))
  loglikely = c(loglike_sym, loglike_allo, sum(loglike_sym, loglike_allo))
  return(loglikely)
}

## example function call
#load("aw2020_billdata.Rdata")
#allo_ed = abs(billdata$ed[billdata$patry==0])
#allo_lat=abs(billdata$lat[billdata$patry==0]) 
#allo_age=billdata$age[billdata$patry==0]
#sym_ed = abs(billdata$ed[billdata$patry==1])
#sym_lat=abs(billdata$lat[billdata$patry==1])
#sym_age=billdata$age[billdata$patry==1]
#simpars = c(0.8, 0.005, 0, 0.5, 2)
#res = simlikely_chardis_bmallo(pars=simpars, N=5000)

#===========================================
# Species Sorting function 1
# Fits an SS model in which trait differentiation results from divergent selection or a random walk process
#' @param sym_div: Vector of trait divergence values for sympatric pairs
#' @param allo_div: Vector of trait divergence values for allopatric pairs
#' @param sym_ages: Vector of species pair ages for sympatric pairs
#' @param allo_ages: Vector of species pair ages for allopatric pairs
#' @param sym_lats: Vector of species pair latitudes (see 'methods') for sympatric pairs
#' @param allo_lats: Vector of species pair latitudes (see 'methods') for allopatric pairs
#' @param meantc_allo: mean time to contact for allopatric pairs (based on their latitude)
#' @param meantc_sym: mean time to contact for sympatric pairs (based on their latitude)
#' @param pars: Vector of model parameter values (see second note below)
#' @param N: number of replicates to simulate for each age in 'ages'
#' @return A vector containing the log-likelihood of models parameters given each pair, plus the value of the sum of the log-likelihoods
# Vector returned of length == length(div) + 1
### NOTE: pars vector must be in the order: c(psi_slope, psi_intercept, alpha, sig2, tc_var, inflection)
simlikely_ss = function(pars, N, allo_ages=allo_age, allo_lats=allo_lat, sym_ages=sym_age, 
  sym_lats=sym_lat, sym_div=sym_ed, allo_div=allo_ed, mean_tc_allo=meantc_allo, exponential=FALSE) {

  # Define/calculate parameters
  if(length(pars)!=6 & length(pars)!=4) {
    stop("Double check your parameter vector; you need either 6 or 4 parameters to run this function")
  }

  if(length(pars) == 6) {
    psi_allo = allo_lats*pars[1] + pars[2]
    psi_sym = sym_lats*pars[1] + pars[2]
    alpha=pars[3]
    sig2=pars[4]
    timepar=pars[5]
    inflection=pars[6]
    # calculate the expected value and variance of trait differentiation for allopatric and sympatric pairs
    u_sym = psi_sym * (1 - exp(-alpha * sym_ages))
    u_allo = psi_allo * (1 - exp(-alpha * allo_ages))
    if(alpha > 0) {
      v_sym = (sig2/(alpha)) * (1 - exp(-2 * alpha * sym_ages))
      v_allo = (sig2/(alpha)) * (1 - exp(-2 * alpha * allo_ages))
    } else {
      v_sym = 2*sig2*sym_ages
      v_allo = 2*sig2*allo_ages
    }
  } 
  if(length(pars) == 4) {
    sig2_allo = allo_lats*pars[1] + pars[2]
    sig2_sym = sym_lats*pars[1] + pars[2]
    timepar=pars[3]
    inflection=pars[4]
    # calculate the expected value and variance of trait differentiation for allopatric and sympatric pairs
    u_sym=rep(0, length(sym_ages))
    u_allo=rep(0, length(allo_ages))
    v_sym=2*sig2_sym*sym_ages
    v_allo=2*sig2_allo*allo_ages
    }

  # Draw N times-to-contact for each allopatric pair
  if(exponential == FALSE) {
    tc_allo = sapply(mean_tc_allo, FUN=truncnorm::rtruncnorm, n=N, a=0, b=Inf, sd=sqrt(timepar))
  } else {
    tc_allo = rexp(n=N, rate=timepar)
  }
  
  ## Draw absolute trait divergences for allopatric pairs
  # And calculate likelihood of the model given each pair
  loglike_allo = rep(NA, length(allo_ages))
  for(i in 1:length(allo_ages)) {
    if(exponential==FALSE) tc=tc_allo[,i]
    if(exponential==TRUE) tc=tc_allo
    nolder = length(tc[tc<allo_ages[i]])
    p_older = nolder/N
    p_younger = 1 - p_older
    divs_allo = rtruncnorm(n=N, a=0, mean=u_allo[i], sd=sqrt(v_allo[i]))
    if(nolder > 0) {
      # find probability of each randomly-drawn trait divergence becoming sympatric
      # p(sympatric) here is function of absolute trait differentiation based on the curve from the cdf of triangle distribution
      # this curve is a piecewise function, so first we break up the randomly drawn divergences into sectors ...
      # of the domain corresponding to p=0, the convex portion, the concave portion, and p=1
      # the inflection point of the curve is at p=0.5 and occurs at the mean absolute value of divergence (u)
      null = divs_allo[divs_allo==0]
      conv = divs_allo[divs_allo > 0 & divs_allo <= inflection]
      conc = divs_allo[divs_allo > inflection & divs_allo < 2*inflection]
      full = divs_allo[divs_allo >= 2*inflection]
      divs_allo = c(null, conv, conc, full)
      psym_nul = rep(0,length(null))
      psym_conv = conv^2/(2*inflection^2)
      psym_conc = 1 - ((2*inflection - conc)^2)/(2*inflection^2)
      psym_ful = rep(1, length(full))
      psyms_allo = c(psym_nul, psym_conv, psym_conc, psym_ful)

      # run Bernoulli trials based on psyms, where a "success" means becoming sympatric
      # grab all values of trait divergence that failed to become sympatric
      bern_trials = rbinom(n=length(divs_allo), 1, prob=psyms_allo)

      # if enough pairs are allopatric after the trial, fit density function to their distribution and calc likelihood of func given real pair
      # if too few pairs are allopatric, density function will be of poor quality (poor fit)
      # in this case, set the likelihood for the model to zero
      # this effectively penalizes models for which p(allo) is too low (<0.01) to account for actual allopatric pairs of a given age
      if(length(which(bern_trials==0))>100) {
        div_older = divs_allo[which(bern_trials==0)]
        dens_older = suppressWarnings(locfit(~lp(div_older), xlim=c(0,10), maxk=200))
        like_older = predict(dens_older, allo_div[i])
        } else {
        like_older = 0 
        }
    }

    if(nolder < N) {
      # fit density function to simulated joint distribution & calc likelihood of func given real pair
      dens_younger = suppressWarnings(locfit(~lp(divs_allo), xlim=c(0,10), maxk=200))
      like_younger = predict(dens_younger, allo_div[i])
    }

    if(nolder > 0 & nolder < N) {
      # weight the likelihoods by the prob of being in either category and sum the marginal densities
      loglike_allo[i] = log(like_older*p_older + like_younger*p_younger)
    }
    if(nolder == 0) {
      loglike_allo[i] = log(like_younger)
    }
    if(nolder == N) {
      loglike_allo[i] = log(like_older)
    }
  }

  ## Draw absolute trait divergences for sympatric pairs
  # And calculate likelihood of the model given each pair
  loglike_sym = rep(NA, length(sym_ages))
  for(i in 1:length(sym_ages)) {
    divs_sym = rtruncnorm(n=N, a=0, mean=u_sym[i], sd=sqrt(v_sym[i]))

    # Find probability of each randomly-drawn trait divergence becoming sympatric
    null = divs_sym[divs_sym==0]
    conv = divs_sym[divs_sym > 0 & divs_sym <= inflection]
    conc = divs_sym[divs_sym > inflection & divs_sym < 2*inflection]
    full = divs_sym[divs_sym >= 2*inflection]
    divs_sym = c(null, conv, conc, full)
    psym_nul = rep(0,length(null))
    psym_conv = conv^2/(2*inflection^2)
    psym_conc = 1 - ((2*inflection - conc)^2)/(2*inflection^2)
    psym_ful = rep(1, length(full))
    psyms_sym = c(psym_nul, psym_conv, psym_conc, psym_ful)

    # run Bernoulli trials based on psyms and grab the divs that are sympatric
    bern_trials = rbinom(n=length(divs_sym), 1, prob=psyms_sym)

    # fit density function to the simulated distribution & calc likelihood of func given real pair
    if(length(which(bern_trials==1))>100) {
      syms = divs_sym[which(bern_trials==1)]
      dens_sym = suppressWarnings(locfit(~lp(syms), xlim=c(0,10), maxk=200))
      loglike_sym[i] = log(predict(dens_sym, sym_div[i]))
      } else {
      loglike_sym[i] = -Inf
      }
  }

  # combine results and calculate likelihood for the given parameter set
  loglikely = c(loglike_sym, loglike_allo)
  loglikely = c(loglikely, sum(loglikely))
  return(loglikely)
}
## example function call
#load("aw2020_billdata.Rdata")
#allo_ed = abs(billdata$ed[billdata$patry==0])
#allo_lat=abs(billdata$lat[billdata$patry==0]) 
#allo_age=billdata$age[billdata$patry==0]
#sym_ed = abs(billdata$ed[billdata$patry==1])
#sym_lat=abs(billdata$lat[billdata$patry==1])
#sym_age=billdata$age[billdata$patry==1]
#meantc_allo=-0.0238*allo_lat + 3.18
#simpars = c(0, 0.8, 0.5, 0.006, 2, 0.2)
#res = simlikely_es_inflection(pars=simpars, N=5000) 
#res = simlikely_es_inflection(pars=sorting_pars[1000,], N=5000) 
# note - I tested this with the old function (from submission 1) and it comes up with the same values, as expected


#===========================================
# Species Sorting function 2
# This is the species sorting function for when no external time to contact is imposed
simlikely_ss_notc = function(pars, N, allo_ages=allo_age, allo_lats=allo_lat, 
  sym_ages=sym_age, sym_lats=sym_lat, sym_div=sym_ed, allo_div=allo_ed) {

  # Define/calculate parameters
  if(length(pars)!=5 & length(pars)!=3) {
    stop("Double check your parameter vector; you need either 5 or 3 parameters to run this function")
  }
  if(length(pars) == 5) {
    psi_allo = allo_lat*pars[1] + pars[2]
    psi_sym = sym_lat*pars[1] + pars[2]
    alpha=pars[3]
    sig2=pars[4]
    inflection=pars[5]
    # calculate the expected value and variance of trait differentiation for allopatric and sympatric pairs
    u_sym = psi_sym * (1 - exp(-alpha * sym_ages))
    u_allo = psi_allo * (1 - exp(-alpha * allo_ages))
    if(alpha > 0) {
      v_sym = (sig2/(alpha)) * (1 - exp(-2 * alpha * sym_ages))
      v_allo = (sig2/(alpha)) * (1 - exp(-2 * alpha * allo_ages))
    } else {
      v_sym = 2*sig2*sym_ages
      v_allo = 2*sig2*allo_ages
    }
  } 
  if(length(pars) == 3) {
    sig2_allo = allo_lat*pars[1] + pars[2]
    sig2_sym = sym_lat*pars[1] + pars[2]
    inflection=pars[3]
    u_sym=rep(0, length(sym_ages))
    u_allo=rep(0, length(allo_ages))
    v_sym=2*sig2_sym*sym_ages
    v_allo=2*sig2_allo*allo_ages
    }

  ## Draw absolute trait divergences for allopatric pairs
  # And calculate likelihood of the model given each pair
  loglike_allo = rep(NA, length(allo_ages))
  for(i in 1:length(allo_ages)) {
    divs_allo = rtruncnorm(n=N, a=0, mean=u_allo[i], sd=sqrt(v_allo[i]))
    null = divs_allo[divs_allo==0]
    conv = divs_allo[divs_allo > 0 & divs_allo <= inflection]
    conc = divs_allo[divs_allo > inflection & divs_allo < 2*inflection]
    full = divs_allo[divs_allo >= 2*inflection]
    divs_allo = c(null, conv, conc, full)
    psym_nul = rep(0,length(null))
    psym_conv = conv^2/(2*inflection^2)
    psym_conc = 1 - ((2*inflection - conc)^2)/(2*inflection^2)
    psym_ful = rep(1, length(full))
    psyms_allo = c(psym_nul, psym_conv, psym_conc, psym_ful)

    # run Bernoulli trials based on psyms, where a "success" means becoming sympatric
    # grab all values of trait divergence that failed to become sympatric
    bern_trials = rbinom(n=length(divs_allo), size=1, prob=psyms_allo)

    # if enough pairs are allopatric after the trial, fit density function to their distribution and calc likelihood of func given real pair
    # if too few pairs are allopatric, density function will be of poor quality (poor fit)
    # in this case, set the likelihood for the model to zero
    # this effectively penalizes models for which p(allo) is too low (<0.01) to account for actual allopatric pairs of a given age
    if(length(which(bern_trials==0))>100) {
      allos = divs_allo[which(bern_trials==0)]
      dens_allo = suppressWarnings(locfit(~lp(allos), xlim=c(0,10), maxk=200))
      loglike_allo[i] = log(predict(dens_allo, allo_div[i]))
    } else {
      loglike_allo[i] = log(0) 
    }
  }

  ## Draw absolute trait divergences for sympatric pairs
  # And calculate likelihood of the model given each pair
  loglike_sym = rep(NA, length(sym_ages))
  for(i in 1:length(sym_ages)) {
    divs_sym = rtruncnorm(n=N, a=0, mean=u_sym[i], sd=sqrt(v_sym[i]))

    # Find probability of each randomly-drawn trait divergence becoming sympatric
    null = divs_sym[divs_sym==0]
    conv = divs_sym[divs_sym > 0 & divs_sym <= inflection]
    conc = divs_sym[divs_sym > inflection & divs_sym < 2*inflection]
    full = divs_sym[divs_sym >= 2*inflection]
    divs_sym = c(null, conv, conc, full)
    psym_nul = rep(0,length(null))
    psym_conv = conv^2/(2*inflection^2)
    psym_conc = 1 - ((2*inflection - conc)^2)/(2*inflection^2)
    psym_ful = rep(1, length(full))
    psyms_sym = c(psym_nul, psym_conv, psym_conc, psym_ful)

    # run Bernoulli trials based on psyms and grab the divs that are sympatric
    bern_trials = rbinom(n=length(divs_sym), size=1, prob=psyms_sym)

    # fit density function to the simulated distribution & calc likelihood of func given real pair
    if(length(which(bern_trials==1))>100) {
      syms = divs_sym[which(bern_trials==1)]
      dens_sym = suppressWarnings(locfit(~lp(syms), xlim=c(0,10), maxk=200))
      loglike_sym[i] = log(predict(dens_sym, sym_div[i]))
      } else {
      loglike_sym[i] = log(0)
      }
  }
  # combine results and calculate likelihood for the given parameter set
  loglikely = c(loglike_sym, loglike_allo)
  loglikely = c(loglikely, sum(loglikely))
  return(loglikely)
}

#===========================================
# Age-difference Function
# Fits an 'age-difference' model
#' @param sym_div: Vector of trait divergence values for sympatric pairs
#' @param allo_div: Vector of trait divergence values for allopatric pairs
#' @param sym_ages: Vector of species pair ages for sympatric pairs
#' @param allo_ages: Vector of species pair ages for allopatric pairs
#' @param sym_lats: Vector of species pair latitudes (see 'methods') for sympatric pairs
#' @param allo_lats: Vector of species pair latitudes (see 'methods') for allopatric pairs
#' @param pars: Vector of model parameter values (see second note below)
#' @param N: number of replicates to simulate for each age in 'ages'
#' @return A vector containing the log-likelihood of models parameters given each pair, plus the value of the sum of the log-likelihoods
# Vector returned of length == length(div) + 1
simlikely_ad = function(pars, N, allo_ages=allo_age, allo_lats=allo_lat, 
  sym_ages=sym_age, sym_lats=sym_lat, sym_div=sym_ed, allo_div=allo_ed) {
  lats=c(sym_lats, allo_lats)
  ages=c(sym_ages, allo_ages)
  div=c(sym_div, allo_div)
  if(length(pars)==2) {
    sig2 = lats*pars[1] + pars[2]
    sd = sqrt(2*ages*sig2)
    u2=0
  } else {
    psi = lats*pars[1] + pars[2]
    alpha = pars[3]
    sig2 = pars[4]
    u2 = psi*(1-exp(-alpha*ages))
    if(alpha > 0) sd = sqrt((sig2/alpha)*(1-exp(-2*alpha*(ages))))
    if(alpha == 0) sd = sqrt(2*ages*sig2)
  }
  simdiv = lapply(sd, FUN=truncnorm::rtruncnorm, n=N, a=0, mean=u2)
  dens_div = suppressWarnings(lapply(simdiv, function(x) locfit(~lp(x), xlim=c(0,10), maxk=200)))
  loglike = rep(0, length(simdiv))
  for(i in 1:length(simdiv)) {
    loglike[i] = log(predict(dens_div[[i]], div[i]))
  }
  loglikely = c(loglike, sum(loglike))
  return(loglikely)
}
## example function call
# load("aw2020_billdata.Rdata")
# simpars = c(0, 0.005) # if testing AD with no selection
# simpars = c(-0.004, 0.5, 0.8, 0.2) # if testing AD with selection
# res = simlikely_ad(pars=simpars, N=5000)
# loglike = res[nrow(billdata)+1]
# loglike


#===========================================
# Model-fitting summarizer function
#' @param parmat: parameter matrix
#' @param reslist: list of results from the model fitting step, where the i'th element in the list corresponds to the likelihoods based on the parameters in the i'th row of the parameter matrix
#' @param sp.pairs: dataset the model was fit to
#' @param rm.outliers: logical indicating whether to calculate fit after removing the most extremely differentiated pairs
#' @return A matrix in which  the sum of log likelihoods matched up to their respective parameter combinations
summarizer = function(parmat, reslist, sps=sp.pairs, rm.outliers=FALSE) {
  res_cd = cbind(parmat, do.call(rbind,reslist))
  colnames(res_cd)[(ncol(parmat)+1):ncol(res_cd)] = c(round(sps, 3), "sumloglike")
  if(rm.outliers==TRUE) {
    x=which(colnames(res_cd)==6.772 | colnames(res_cd) == 5.343)
    cd_full = res_cd[,-c(x, ncol(res_cd))]
    sumloglike = apply(cd_full, MARGIN=1, function(x) sum(x[(ncol(parmat)+1):length(x)]))
    cd_full = cbind(cd_full[,1:ncol(parmat)], sumloglike)
  } else {
    cd_full = res_cd[,c(1:ncol(parmat),ncol(res_cd))]
  }
  k = as.numeric(apply(cd_full, MARGIN=1, function(x) ncol(parmat)-sum(x==0)))
  AIC = 2*k - 2*cd_full[,"sumloglike"]
  if(rm.outliers==TRUE) AICc = AIC + 2*k*(k+1)/(length(sp.pairs)-2-k-1)
  if(rm.outliers==FALSE) AICc = AIC + 2*k*(k+1)/(length(sp.pairs)-k-1)
  cd_full = cbind(cd_full, k, AIC, AICc)
  cd_full = cd_full[order(cd_full[,"AICc"], decreasing=F),]
  return(cd_full)
}





