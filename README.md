# structural_forecast_BVAR
User Guide for "Macroeconomic Effects of Economic Policy Uncertainty Shocks in US"
Last version: June 08, 2023

Contact: 
Renato Vassallo | renato.vassallo@bse.eu
Barcelona School of Economics

System Requirements 
This system is tested compatible MATLAB R2022a/b. 

0. "data" and "functions"
--------------------------

It is important at the beginning to add to the path of the working directory the folders "data" and "functions", which contain the database used for the estimation; as well as useful functions that serve as inputs into the generated code routines.


The package contains three main m files that are described below. 

1. bvar_estimation_irf
-----------------------

Contains the routine for estimation, computation of impulse responses and generation of graphs with the posterior distribution of the parameters.

Estimates a Bayesian VAR model with independent Normal-Wishart Priors. The Gibbs Sampling algorithm performs 40,000 simulations and discards the first 30,000 (burn-in).


2. bvar_uncond_forecast
-----------------------

Contains the code to perform an unconditional forecast 12 periods ahead. The routine was modified to estimate the parameters with information up to 2019, but the forecast is made using March 2023 as the last observed point.


3. bvar_cond_forecast
-----------------------

Make a structural forecast 5 periods ahead, where a fixed path is assumed for the first variable of the VAR (EPU). The routine includes the generation of Fan Charts for both the unconditional and structural forecast.
