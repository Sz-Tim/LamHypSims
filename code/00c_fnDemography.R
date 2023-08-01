# KELPER
# Demography functions
# Tim Szewczyk


# Demographic functions used by the KELPER population models





# growth & loss -----------------------------------------------------------

#' Calculate growth in frond area
#' 
#' Calculates Frond Area Index at the end of the time step
#'
#' @param FAI_orig Vector of initial frond areas
#' @param N_orig Vector of initial population abundances
#' @param A.mx Transition matrix
#' @param kappa FAI_orig/FAI_K
#' @param logAreaFrond.stage Vector of log frond areas calculated allometrically for each size class
#' @param growRateFrond Vector of frond growth rate for each size class
#'
#' @return
#' @export
#'
#' @examples
growFrondArea <- function(FAI_orig, N_orig, A.mx, kappa, logAreaFrond.stage, growRateFrond) {
  # Laminaria hyperborea exhibits a 'May cast', where the previous year's growth 
  # is shed essentially in its entirety. For other less deciduous species, 
  # FAI_new[j] would include FAI_orig[j] * A.mx[j,j] in the sum
  FAI_new <- rep(0, length(N_orig))
  # recruits: (surviving FAI) + (new growth within stage)
  FAI_new[1] <- N_orig[1]*A.mx[1,1]*growRateFrond[1]*(1-kappa)
  for(j in 2:length(N_orig)) {
    # others: (surviving FAI) + (new growth within stage) + (FAI of newcomers)
    FAI_new[j] <- N_orig[j]*A.mx[j,j]*growRateFrond[j] +
      N_orig[j-1]*A.mx[j,j-1]*exp(logAreaFrond.stage[j])
  }
  
  return(FAI_new)
}










# allometry ---------------------------------------------------------------


#' Calculate biomass across years, seasons, and stages
#'
#' @param N Array with densities
#' @param FAI Array with frond area indexes
#' @param lwtStipe Vector of log stipe weight for each stage
#' @param lmFit Fitted model for logWtFrond ~ logAreaFrond
#' @param ndraws Number of draws from the posterior distribution to use
#' @param env.df Dataframe with covariates for input
#' @param scale.df Dataframe with mean, sd to de-center and de-scale
#' @param stages Vector (integer) with stages to use in biomass calculation
#'
#' @return
#' @export
#'
#' @examples
calcBiomass <- function(N, FAI, lwtStipe, lmFit, ndraws, env.df, scale.df, stages=3) {
  stipeMass <- frondMass <- array(0, dim=dim(N))
  biomass <- matrix(0, nrow=dim(N)[2], ncol=dim(N)[3])
  for(year in 1:(dim(N)[2])) {
    stipeMass[,year,] <- N[,year,] * exp(lwtStipe[year,])
    
  }
  for(stage in stages) {
    for(season in 1:(dim(N)[3])) {
      frondMass[stage,,season] <- N[stage,,season] * 
        exp(getPrediction(lmFit, ndraws, 
                          bind_cols(logAreaFrond=log(FAI[stage,,season]/N[stage,,season]), 
                                    env.df),
                          scale.df, "logWtFrond"))
    }
  }
  if(any(is.na(frondMass))) {
    mass.na <- which(is.na(frondMass), arr.ind=T)
    for(i in 1:dim(mass.na)[1]) {
      frondMass[mass.na[i,1], mass.na[i,2], mass.na[i,3]] <- 0
    }
  }
  for(year in 1:(dim(N)[2])) {
    for(season in 1:(dim(N)[3])) {
      biomass[year,season] <- sum(stipeMass[stages,year,season]) + sum(frondMass[stages,year,season])
    }
  }
  
  return(biomass/1e3)
}









#' Convert between metrics of kelp size 
#'
#' @param orig Vector of original measurements
#' @param beta_mn Mean intercept + slope(s) for species
#' @param shape Shape of allometric relationship ("linear", "quadratic", "exponential", "log")
#' @param resid_err Residual standard deviation among individuals (default:
#'   NULL). If included, results are stochastic with Norm(0, resid_err) added to
#'   each individual
#' @param beta_sd Standard deviation in beta estimates (default: NULL). If
#'   included, results are stochastic with beta ~Norm(beta_mn, beta_sd)
#'
#' @return
#' 
convertAllometry <- function(orig, beta_mn, shape="linear", 
                             resid_err=NULL, beta_sd=NULL) {
  
  # Incorporate parameter uncertainty
  if(!is.null(beta_sd)) {
    beta <- rnorm(length(beta_mn), beta_mn, beta_sd)
  } else {
    beta <- beta_mn
  }
  
  # Incorporate environmental stochasticity
  if(!is.null(resid_err)) {
    err <- rnorm(length(orig), 0, resid_err)
  } else {
    err <- 0
  }
  
  switch(shape,
         linear=beta[1] + beta[2]*orig + err,
         quadratic=beta[1] + beta[2]*orig + beta[3]*orig^2 + err,
         exponential=exp(beta[1] + beta[2]*orig + err),
         log=log(beta[1] + beta[2]*orig + err))
}




