## Below are parameters used in likelihood searches for CD, SS, and AD models in Anderson and Weir 2021

# ========================================================
## CHARACTER DISPLACEMENT ##

# A: Assuming a truncated normal for the time to sympatry (Ts)
# 1. Under regimes of adaptation

# param definitions:
# alpha: constraint parameter in OU process (roughly analagous to strength of selection)
# sig2: Brownian motion dispersion parameter
# ts_var: variance in the time to sympatry
# psi_allo_eq: psi for allopatric pairs at the equator
# psi_allo_polar: psi for allopatric pairs at the poles
# psi_sym_eq: psi for sympatric pairs at the equator
# psi_sym_polar: psi for sympatric pairs at the poles
# NOTE: psi_polar is not an actual parameter, but is used to calculate the slope of psi, which is a parameter

alpha = 0.000 0.001 0.010 0.020 0.030 0.040 0.050 0.060 0.070 0.080 0.090 
	0.100 0.100 0.150 0.200 0.250 0.300 0.350 0.500 0.700 0.900 1.000 1.500 2.000 2.500
ts_var = 1 2 3 4
sig2 = 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.009 0.010 0.050
    0.100 0.300 0.500 0.700 0.900 1.000 1.500 2.000 2.500
psisym_int = 0.000 0.001 0.005 0.010 0.050 0.100 0.200 0.300 0.400 0.500 0.600 0.700
    0.800 0.900 1.000 2.000 3.000 4.000 5.000 6.000 7.000
psiallo_int = 0.000 0.001 0.002 0.004 0.005 0.006 0.008 0.010 0.025 0.050 0.075 0.100
    0.500 1.000 5.000
psiallo_pol = 0.000 0.001 0.002 0.004 0.005 0.006 0.008 0.010 0.025 0.050 0.075 0.100
    0.500 1.000 5.000
psisym_pol  = 0.000 0.001 0.005 0.010 0.050 0.100 0.200 0.300 0.400 0.500 0.600 0.700
   0.800 0.900 1.000 2.000 3.000 4.000 5.000

# 2. Under regimes of unconstrained random walks in the allopatric phase only 

alpha = 0.00010 0.00025 0.00050 0.00075 0.00100 0.01000 0.05000 0.10000 0.20000
   0.30000 0.40000 0.50000 0.60000 0.70000 0.80000 0.90000 1.00000 1.20000
   1.30000 1.40000 1.50000 2.00000 2.50000 3.00000 3.50000 4.00000 4.50000
   5.00000 6.00000 7.00000 8.00000
sig2 = 0.001 0.005 0.006 0.007 0.008 0.009 0.013 0.014 0.015 0.016 0.017 0.018
   0.019 0.020 0.050
ts_var = 4.00 5.00 6.00 7.00 7.75
psi_int = 0.0000 0.0025 0.0050 0.0075 0.0250 0.0400 0.0500 0.0600 0.0750
   0.1000 0.1250 0.1500 0.1750 0.2000 0.4000 0.6000 0.8000 1.0000
   2.0000 3.0000 4.0000 5.0000 6.0000 7.0000 8.0000 9.0000 10.0000
psi_pol = 0.000 0.002 0.003 0.005 0.007 0.008 0.025 0.040 0.050 0.060
   0.075 0.100 0.125 0.150 0.175 0.200 0.400 0.600 0.800 1.000
  2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000 10.000

# 3. Under regimes of unconstrained random walks (i.e. BM) in both the allopatric and sympatric phase

sig2allo_int = 0.001 0.005 0.006 0.007 0.008 0.009 0.010 0.020 0.030 0.040 0.050 0.100
  0.500 0.900
sig2sym_int = 5e-04 6e-04 7e-04 8e-04 9e-04 1e-03 3e-03 4e-03 5e-03 6e-03 7e-03 8e-03
  9e-03 1e-02 2e-02 3e-02 4e-02 5e-02 6e-02 7e-02 8e-02 9e-02 1e-01 5e-01 9e-01
sig2ts_var = 4.00 5.00 6.00 7.00 7.75
sig2allo_pol = 0.001 0.005 0.006 0.007 0.008 0.009 0.050 0.100 0.300 0.500 0.700 0.900
 1.000 1.500 2.000 2.500
sig2sym_pol = 0.001 0.005 0.006 0.007 0.008 0.009 0.050 0.100 0.300 0.500 0.700 0.900
1.000 1.500 2.000 2.500

# B: Assuming selection and an exponential distribution for the time to sympatry (Ts)
psiallo_int = 0.000 0.001 0.005 0.010 0.025 0.100 0.500 1.000 2.500 5.000
psiallo_pol = 0.000 0.001 0.005 0.010 0.025 0.100 0.500 1.000 2.500 5.000
psisym_int = 0.000 0.001 0.005 0.010 0.025 0.100 0.500 1.000 2.500 5.000
psisym_pol = 0.000 0.001 0.005 0.010 0.025 0.100 0.500 1.000 2.500 5.000
tc_rate = 0.10 0.25 0.50 0.75 1.00 3.00
alpha = 0.001 0.010 0.050 0.100 0.300 0.500 0.700 0.900 1.000 2.000 5.000
sig2 = 0.001 0.005 0.006 0.007 0.008 0.009 0.050 0.100 0.500 2.000


# ========================================================
## SPECIES SORTING ##

# param definitions:
# alpha: constraint parameter in OU process (roughly analagous to strength of selection)
# sig2: Brownian motion dispersion parameter
# ts_var: variance in the time to sympatry
# psi_int: psi at the equator
# psi_pol: psi at the poles
# sig2_int: the BM dispersion parameter at the equator (in models of SS with BM only)
# sig2_polar: the BM dispersion parameter at the poles (in models of SS with BM only)
# NOTE: psi_polar is not an actual parameter, but is used to calculate the slope of psi, which is a parameter

# A: assuming a truncated normal distribution for the time to contact (Tc)

# under a regime of selection 
alpha = 0.0000  0.0001  0.0005  0.0009  0.0010  0.0100  0.0500  0.1000  0.3000
  0.5000  0.7000  0.9000  1.0000  1.5000  2.0000  2.5000 10.0000
sig2 = 0.001 0.005 0.006 0.007 0.008 0.009 0.050 0.100 0.300 0.500 0.700 0.900
  1.000 2.500 5.000
psi_int = 0.000 0.001 0.005 0.010 0.030 0.050 0.070 0.090 0.100 0.200 0.300 0.400
  0.500 1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000
psi_pol = 0.000 0.001 0.005 0.010 0.030 0.050 0.070 0.090 0.100 0.200 0.300 0.400
  0.500 1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000
inflection = 0.0005 0.0010 0.0025 0.0040 0.0050 0.0100 0.0250 0.0400 0.0500 0.0600
  0.0700 0.0800 0.0900 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000
tc_var = 0.10 0.25 0.50 0.75 1.00 2.00 3.00 5.00 6.50 7.75

# under a regime of unconstrained random walks 
tc_var = 0.50 1.00 1.50 2.00 2.50 3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50 7.00 7.50 7.75
inflection = 0.01 0.02 0.03 0.04 0.05 0.06 0.10 0.20 0.30 0.40 0.50 0.60
sig2_int = 0.001 0.005 0.006 0.007 0.008 0.009 0.050 0.100 0.300 0.500 0.700 0.900
  1.000 1.500 2.000 2.500
sig2_pol = 0.001 0.005 0.006 0.007 0.008 0.009 0.050 0.100 0.300 0.500 0.700 0.900
  1.000 1.500 2.000 2.500

#ES_notc
psi_int = 0.00000 0.00025 0.00050 0.00075 0.00100 0.00500 0.01000 0.03000 0.05000
  0.07000 0.09000 0.10000 0.20000 0.30000 0.40000 0.50000 1.00000 5.00000
psi_pol = 0.000 0.001 0.005 0.010 0.030 0.050 0.070 0.090 0.100 0.200 0.300 0.400
  0.500 1.000 5.000
alpha = 0.000  0.001  0.003  0.005  0.007  0.009  0.010  0.030  0.050  0.070
  0.090  0.100  0.200  0.300  0.400  0.500  0.600  0.700  0.800  0.900
  1.000  1.100  1.200  1.300  1.400  1.500  1.600  1.700  1.800  1.900
  2.000  2.500 10.000
sig2 = 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.010 0.011 0.012
  0.013 0.014 0.015 0.016 0.017 0.018 0.019 0.020 0.030 0.040 0.050 0.060
  0.100 0.300 0.500 0.700 0.900 1.000 1.500 2.000 2.500
inflection = 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.7 0.9 1.1 1.3 1.5 1.7 1.9

# B: Assuming an exponential distribution for the time to contact (Tc)
alpha =  0.0007  0.0008  0.0009  0.0010  0.0020  0.0030  0.0040  0.0050  0.0060
  0.0070  0.0080  0.0090  0.0100  0.0200  0.0300  0.0400  0.0500  0.1000
  0.2000  0.3000  0.4000  0.5000  0.7000  0.9000  1.0000  1.5000  2.0000
  2.5000 10.0000
sig2 = 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.010 0.020 0.030
  0.040 0.050 0.100 0.300 0.500 0.700 0.900 1.000 1.500 2.000 2.500
tc_rate = 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.75 1.00 3.00
psi_int = 0.000 0.001 0.005 0.010 0.030 0.050 0.070 0.090 0.100 0.200 0.300 0.400
  0.500 1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000
psi_polar = 0.000 0.001 0.005 0.010 0.030 0.050 0.070 0.090 0.100 0.200 0.300 0.400
  0.500 1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000
inflection = c((0.02,0.09,0.02),(0.1,0.6,0.1))


# ========================================================
## AGE DIFFERENCE ##

# Note: AD models are often a special case of CD models (in which psi_sym == psi_allo at all latitudes)
# and so parameters are included in the above
# However, AD can also be a purely BM process, for which case we explored the following parameters:
sig2_int = c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(1, 10, 1), 100)
sig2_polar = sig2_int

# We also added additional fine-grained parameters for a regime of age difference under selection
psi_int = 0.000 0.001 0.005 0.010 0.030 0.050 0.070 0.090 0.100 0.200 0.300 0.400
  0.500 1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000
psi_pol = 0.000 0.001 0.005 0.010 0.030 0.050 0.070 0.090 0.100 0.200 0.300 0.400
  0.500 1.000 2.000 3.000 4.000 5.000 6.000 7.000 8.000 9.000
alpha = 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009 0.020 0.030 0.040 0.060
  0.070 0.080 0.090 0.150 0.200 0.250 0.350 0.400 0.450 0.550 0.600 0.650
  0.750 0.800 0.850
sig2 = 0.001 0.005 0.006 0.007 0.008 0.009 0.050 0.100


# ========================================================
## CALCULATING SLOPES OF PARAMETERS FROM ENDPOINTS ##

# e.g. with psi_allo from the AD model:
require(diverge)
sig2_int = c(seq(0.001, 0.009, 0.001), seq(0.01, 0.99, 0.01), seq(1, 10, 1), 100)
sig2_polar = sig2_int
AD_pars = model_generator(domain=c(0, 60), par_low = sig2_int, par_high = sig2_polar)
head(AD_pars)





