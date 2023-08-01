# KELPER
# Simulation functions
# Tim Szewczyk


# Population simulation function

#' Simulate one population
#'
#' @param pars List of parameter values
#' @param N0 Vector of initial density for each stage
#' @param ndraws Number of draws from the posterior distribution to use
#'
#' @return
#' @export
#'
#' @examples
simulatePopulation <- function(pars, N0=NULL, ndraws=4e3) {
  
  library(tidyverse); library(brms)
  #---- setup landscape
  env.df <- pars$env %>% 
    mutate(PAR_atDepth=PAR * exp(-KD * pars$depth),
           lPAR_atDepth=log(PAR),
           location=NA) %>%
    select(PAR_atDepth, lPAR_atDepth, SST, fetch, location)
  if(nrow(pars$env) == 1) {
    env.df <- env.df %>% uncount(pars$tmax)
  }
  #---- setup parameters
  if(pars$stochParams) {
    par.yr <- list(loss=rbeta(pars$tmax, prod(pars$lossRate), (1-pars$lossRate[,1])*pars$lossRate[,2]),
                   settlement=pmax(0, rnorm(pars$tmax, pars$settlementRate[,1], pars$settlementRate[,2])),
                   surv=apply(pars$survRate, 1, function(x) rbeta(pars$tmax, prod(x), (1-x[1])*x[2])),
                   growStipeMax=apply(pars$growthRateStipeMax, 1, function(x) rnorm(pars$tmax, x[1], x[2])),
                   growFrond=apply(pars$growthRateFrond, 1, function(x) rnorm(pars$tmax, x[1], x[2])))
  } else {
    par.yr <- list(loss=rep(pars$lossRate[1], pars$tmax),
                   settlement=rep(pars$settlementRate$mn[1], pars$tmax),
                   surv=apply(pars$survRate, 1, function(x) rep(x[1], pars$tmax)),
                   growStipeMax=apply(pars$growthRateStipeMax, 1, function(x) rep(x[1], pars$tmax)),
                   growFrond=apply(pars$growthRateFrond, 1, function(x) rep(x[1], pars$tmax)))
  }
  #---- setup storm effects
  if(is.null(pars$stormIntensity)) {
    par.yr$surv_strm <- par.yr$surv
  } else {
    # vector of storm intensities affects winter survival, loss
    par.yr$loss <- qbeta(pnorm(pars.sim$storms, 0, 1), 
                         prod(pars$lossRate), 
                         (1-pars$lossRate[,1])*pars$lossRate[,2])
    par.yr$surv_strm <- apply(pars$survRate, 1, 
                              function(x) qbeta(pnorm(-pars.sim$storms, 0, 1), 
                                                prod(x), (1-x[1])*x[2]))
  }
  par.yr$surv <- sqrt(par.yr$surv) # annual rates to 1/2 year rates
  par.yr$surv_strm <- sqrt(par.yr$surv_strm)
  
  #---- global parameters
  ## maxStipeLen = maximum expected canopy height
  ## sizeClassLimits = boundaries between stages
  ## sizeClassMdPts = midpoint per stage
  ## K_N = carrying capacity on abundance / m2
  ## K_FAI = carrying capacity on frond area / m2
  maxStipeLen <- getPrediction(pars$canopyHeight, ndraws, 
                               env.df %>% summarise(across(.fns=mean)) %>% mutate(location=NA_character_),
                               pars$sc.df$canopyHeight.lm, "maxStipeLen")
  sizeClassLimits <- maxStipeLen * c(0, 0.333, 0.75, 1.25)#(0:3)/3
  sizeClassMdpts <- (sizeClassLimits+lag(sizeClassLimits))[-1]/2
  
  K_N <- pmax(1e-2, getPrediction(pars$N_canopy, ndraws, env.df,
                                  pars$sc.df$N_canopy.lm, "N"))
  K_FAI <- pmax(1e-2, getPrediction(pars$FAI, ndraws, env.df,
                                    pars$sc.df$FAI.lm, "FAI"))
  if(is.null(N0)) N0 <- K_N[1] * c(5,2,1)/2
  
  
  #---- per capita mass, area by stage
  logWtStipe.stage <- log(sizeClassMdpts) %>%
    map(~getPrediction(pars$lenSt_to_wtSt, ndraws, bind_cols(logLenStipe=.x, env.df),
                       pars$sc.df$lenSt_to_wtSt.lm, "logWtStipe")) %>%
    do.call('cbind', .)
  logWtFrond.stage <- log(sizeClassMdpts) %>%
    map(~getPrediction(pars$lenSt_to_wtFr, ndraws, bind_cols(logLenStipe=.x, env.df),
                       pars$sc.df$lenSt_to_wtFr.lm, "logWtFrond")) %>%
    do.call('cbind', .)
  logAreaFrond.stage <- log(sizeClassMdpts) %>%
    map(~getPrediction(pars$wtFr_to_arFr, ndraws, bind_cols(logWtFrond=.x, env.df),
                       pars$sc.df$wtFr_to_arFr, "logAreaFrond")) %>%
    do.call('cbind', .)
  
  
  #---- storage & initialization
  ## N[stage,year,season] = density/m2
  ## FAI[stage,year,season] = frond area / m2
  ## harvest[year] = log(grams harvested) / m2
  ## kappa[year,season,FAI|N] = proportion of K
  ## biomass[year,season] = kg / m2
  N <- FAI <- array(dim=c(3, pars$tmax, 3))
  N[,1,1] <- N0
  FAI[,1,1] <- N[,1,1] * exp(logAreaFrond.stage[1,])
  harvest <- rep(0, pars$tmax) 
  kappa <- array(dim=c(pars$tmax, 3, 2))
  
  
  #---- transition matrices
  ## A[[growing|non-growing]][stageTo,stageFrom]
  A <- map(1:2, ~matrix(0, 3, 3))
  
  
  #---- simulation loop
  for(year in 1:pars$tmax) {
    harvestYear <- year %% pars$freqHarvest == 0
    
    #---- growing season: 
    season <- 1
    kappa[year,season,] <- pmin(1, c(FAI[3,year,season]/K_FAI[year], N[3,year,season]/K_N[year]))
    
    # growth
    growRate_i <- par.yr$growStipeMax[year,] * (1-kappa[year,season,1]^pars$growthRateDensityShape)
    prGrowToNext <- pmin(1, pmax(0, growRate_i/(sizeClassLimits-lag(sizeClassLimits))[-1]))
    A[[1]][2,1] <- par.yr$surv[year,1]*prGrowToNext[1]
    A[[1]][3,2] <- par.yr$surv[year,2]*prGrowToNext[2]
    # survival
    A[[1]][1,1] <- par.yr$surv[year,1] - A[[1]][2,1]
    A[[1]][2,2] <- par.yr$surv[year,2] - A[[1]][3,2]
    A[[1]][3,3] <- par.yr$surv[year,3]
    # update population
    N[,year,season+1] <- A[[1]] %*% N[,year,season]
    FAI[,year,season+1] <- growFrondArea(FAI[,year,season], N[,year,season], 
                                         A[[1]], kappa[year,season,1],
                                         logAreaFrond.stage[year,], par.yr$growFrond[year,])
    
    
    #---- harvest
    season <- 2
    kappa[year,season,] <- pmin(1, c(FAI[3,year,season]/K_FAI[year], N[3,year,season]/K_N[year]))
    
    # harvest
    if(harvestYear) {
      N[,year,season+1] <- (1-pars$prFullHarvest) * N[,year,season]
      FAI[,year,season+1] <- (1-pars$prFullHarvest) * FAI[,year,season]
    } else {
      N[,year,season+1] <- N[,year,season]
      FAI[,year,season+1] <- FAI[,year,season]
    }
    
    
    #---- non-growing season
    season <- 3
    kappa[year,season,] <- pmin(1, c(FAI[3,year,season]/K_FAI[year], N[3,year,season]/K_N[year]))
    
    if(year < pars$tmax) {
      # survival
      diag(A[[2]]) <- pmax(0, par.yr$surv_strm[year,])
      # update population
      N[,year+1,1] <- A[[2]] %*% N[,year,season]
      # reproduction
      N[1,year+1,1] <- par.yr$settlement[year]*(1-max(kappa[year,season,]))
      FAI[,year+1,1] <- FAI[,year,season] * pmax(0, diag(A[[2]]) - par.yr$loss[year])
    }
    
  }
  
  
  # biomass calculation
  biomass <- calcBiomass(N, FAI, logWtStipe.stage, pars$arFr_to_wtFr,
                         ndraws, env.df, pars$sc.df$arFr_to_wtFr.lm, stages=3)
  
  return(list(N=N, FAI=FAI, harvest=harvest, kappa=kappa, K_FAI=K_FAI, K_N=K_N,
              biomass=biomass, PAR=env.df$PAR_atDepth))
}







# simulation wrappers -----------------------------------------------------


#' Simulate multiple depths within one cell
#'
#' @param x Index for parLapply
#' @param grid.i Dataframe with row for each grid cell
#' @param grid.sim Simulated grid
#' @param grid.id Simulated grid id
#' @param gridRes Grid resolution
#' @param pars.sim List of parameters
#' @param surv.df Survival rate dataframe
#' @param fecund.df Settlement dataframe
#' @param lm.fit List of model fits
#' @param lm.mnsd List of mn, sd to de-center and de-scale lm predictions
#'
#' @return
#' @export
#'
#' @examples
simDepthsWithinCell <- function(x, grid.i, grid.sim=NULL, grid.id=NA, gridRes, 
                                pars.sim, surv.df, fecund.df, lm.fit, lm.mnsd, par.dir) {
  library(glue); library(lubridate); library(sf); library(brms)
  library(Matrix); library(tidyverse)
  sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")
  
  if(pars.sim$landscape=="dynamic") {
    cell.env <- grid.sim %>% filter(id==grid.i$id[x])
  } else {
    cell.env <- grid.i[x,]
  }
  
  params <- dir(par.dir, "par_") %>%
    map(~suppressMessages(read_csv(glue(par.dir, .x))) %>%
          mutate(exposure=case_when(exposure=="low"~1, exposure=="high"~2))) %>%
    setNames(str_sub(dir(par.dir, "par_"), 5, -5))
  
  pop.ls <- mass.ls <- vector("list", length(pars.sim$depths))
  for(j in 1:length(pars.sim$depths)) {
    par_i <- setParameters(
      tmax=pars.sim$tmax, 
      survRate=params$survival %>% 
        filter(exposure==cell.env$fetchCat[1]) %>%
        select(mn, prec) %>% cbind,
      settlementRate=params$recruitment %>% 
        filter(exposure==cell.env$fetchCat[1]) %>%
        select(mn, sd) %>% cbind,
      lossRate=params$erosion %>% select(mn, prec) %>% cbind,
      stochParams=pars.sim$stochParams,
      stormIntensity=pars.sim$storms,
      extraPars=list(
        depth=pars.sim$depths[j],
        env=cell.env,
        lenSt_to_wtSt=lm.fit$lenSt_to_wtSt.lm,
        lenSt_to_wtFr=lm.fit$lenSt_to_wtFr.lm,
        wtFr_to_arFr=lm.fit$wtFr_to_arFr.lm,
        arFr_to_wtFr=lm.fit$arFr_to_wtFr.lm,
        N_canopy=lm.fit$N_canopy.lm,
        FAI=lm.fit$FAI.lm,
        canopyHeight=lm.fit$canopyHeight.lm,
        sc.df=lm.mnsd)
    )
    out <- map(1:pars.sim$nSim, 
               ~simulatePopulation(par_i, ndraws=ifelse(par_i$stochParams, 50, 4e3)))
    pop.ls[[j]] <- imap_dfr(out, 
                            ~tibble(sim=.y, 
                                    year=rep(1:par_i$tmax, 3),
                                    month=rep(c(1,6,7), each=par_i$tmax),
                                    date=ymd(glue("{year+1942}-{month}-01")),
                                    PAR_atDepth=rep(.x$PAR, 3),
                                    K_N=rep(.x$K_N, 3),
                                    K_FAI=rep(.x$K_FAI, 3),
                                    FAI=c(.x$FAI[3,,]),
                                    N.recruits=c(.x$N[1,,]),
                                    N.subcanopy=c(.x$N[2,,]),
                                    N.canopy=c(.x$N[3,,]),
                                    kappa_FAI=c(.x$kappa[,,1]),
                                    kappa_N=c(.x$kappa[,,2])) %>%
                              pivot_longer(contains("N."), names_to="stage", values_to="N") %>%
                              mutate(stage=factor(str_sub(stage, 3, -1), 
                                                  levels=c("canopy", "subcanopy", "recruits")))) %>%
      mutate(id=x) %>%
      left_join(., par_i$env) %>%
      mutate(depth=par_i$depth,
             landscape=pars.sim$landscape,
             stochParams=pars.sim$stochParams,
             grid.id=grid.id)
    mass.ls[[j]] <- imap_dfr(out,
                             ~tibble(sim=.y,
                                     year=rep(1:par_i$tmax, 3),
                                     month=rep(c(1,6,7), each=par_i$tmax),
                                     date=ymd(glue("{year+1942}-{month}-01")),
                                     PAR_atDepth=rep(.x$PAR, 3),
                                     FAI=c(.x$FAI[3,,]),
                                     biomass=c(.x$biomass),
                                     kappa_FAI=c(.x$kappa[,,1]),
                                     kappa_N=c(.x$kappa[,,2]))) %>%
      mutate(id=x) %>%
      left_join(., par_i$env) %>%
      mutate(depth=par_i$depth,
             landscape=pars.sim$landscape,
             stochParams=pars.sim$stochParams,
             grid.id=grid.id)
  }
  sim.info <- glue("{str_pad(x,4,'left','0')}_{gridRes}_{pars.sim$landscape}",
                   "_{ifelse(pars.sim$stochParams, 'stochPar', 'optPar')}",
                   "_g{str_pad(grid.id,3,'left','0')}")
  saveRDS(do.call(rbind, pop.ls), glue("out{sep}storms{sep}pop_{sim.info}.rds"))
  saveRDS(do.call(rbind, mass.ls), glue("out{sep}storms{sep}mass_{sim.info}.rds"))
  return(x)
}









#' Simulate a range of parameters across depths within a cell
#'
#' @param x Index for parLapply
#' @param grid.i Dataframe with row for each grid cell
#' @param gridRes Grid resolution
#' @param pars.sens List of parameters
#' @param lm.fit List of model fits
#' @param lm.mnsd List of mn, sd to de-center and de-scale lm predictions
#' @param parSets List of parameter values
#'
#' @return
#' @export
#'
#' @examples
simSensitivityDepthsWithinCell <- function(x, grid.i, gridRes, pars.sens, 
                                           lm.fit, lm.mnsd, parSets) {
  library(glue); library(lubridate); library(sf); library(brms)
  library(Matrix); library(tidyverse)
  sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")
  
  # setup landscape and parameters
  cell.env <- grid.i[x,]
  grStipe.draws <- parSets[[3]] %>% filter(param=="growStipe") %>% group_split(parDraw)
  grFrond.draws <- parSets[[3]] %>% filter(param=="growFrond") %>% group_split(parDraw)
  loss.draws <- parSets[[3]] %>% filter(param=="loss")
  surv.draws <- parSets[[cell.env$fetchCat]] %>% filter(param=="surv") %>% group_split(parDraw)
  settle.draws <- parSets[[cell.env$fetchCat]] %>% filter(param=="settlement")
  densShape.draws <- parSets[[3]] %>% filter(param=="densityEffShape")
  
  pop.depth <- mass.depth <- vector("list", length(pars.sens$depths))
  for(j in 1:length(pars.sens$depths)) {
    pop.ls <- mass.ls <- vector("list", pars.sens$nParDraws)
    for(k in 1:pars.sens$nParDraws) {
      par_i <- setParameters(
        tmax=pars.sens$tmax, 
        growthRateStipeMax=cbind(grStipe.draws[[k]]$val),
        growthRateFrond=cbind(grFrond.draws[[k]]$val),
        lossRate=loss.draws$val[k],
        survRate=cbind(surv.draws[[k]]$val),
        settlementRate=settle.draws$val[k],
        growthRateDensityShape=densShape.draws$val[k],
        stochParams=pars.sens$stochParams,
        stormIntensity=pars.sens$storms,
        extraPars=list(
          depth=pars.sens$depths[j],
          env=cell.env,
          lenSt_to_wtSt=lm.fit$lenSt_to_wtSt.lm,
          lenSt_to_wtFr=lm.fit$lenSt_to_wtFr.lm,
          wtFr_to_arFr=lm.fit$wtFr_to_arFr.lm,
          arFr_to_wtFr=lm.fit$arFr_to_wtFr.lm,
          N_canopy=lm.fit$N_canopy.lm,
          FAI=lm.fit$FAI.lm,
          canopyHeight=lm.fit$canopyHeight.lm,
          sc.df=lm.mnsd)
      )
      out <- map(1:pars.sens$nSim, ~simulatePopulation(par_i, ndraws=4e3))
      pop.ls[[k]] <- imap_dfr(out, 
                              ~tibble(sim=.y, 
                                      year=rep(1:par_i$tmax, 3),
                                      month=rep(c(1,6,7), each=par_i$tmax),
                                      date=ymd(glue("{year}-{month}-01")),
                                      N.recruits=c(.x$N[1,,]),
                                      N.subcanopy=c(.x$N[2,,]),
                                      N.canopy=c(.x$N[3,,])) %>%
                                pivot_longer(contains("N."), names_to="stage", values_to="N") %>%
                                mutate(stage=factor(str_sub(stage, 3, -1), 
                                                    levels=c("canopy", "subcanopy", "recruits")))) %>%
        filter(month != 6 & year > pars.sens$tskip) %>%
        mutate(N=replace(N, N<=0, 1e-4)) %>%
        group_by(sim, month, stage) %>%
        summarise(N_mn=mean(N), N_md=median(N), N_sd=sd(N),
                  lN_mn=mean(log(N+1)), lN_md=median(log(N+1)), lN_sd=sd(log(N+1))) %>%
        mutate(parDraw=k)
      mass.ls[[k]] <- imap_dfr(out,
                               ~tibble(sim=.y,
                                       year=rep(1:par_i$tmax, 3),
                                       month=rep(c(1,6,7), each=par_i$tmax),
                                       date=ymd(glue("{year}-{month}-01")),
                                       PAR_atDepth=rep(.x$PAR, 3),
                                       K_N=rep(.x$K_N, 3),
                                       K_FAI=rep(.x$K_FAI, 3),
                                       FAI=c(.x$FAI[3,,]),
                                       biomass=c(.x$biomass),
                                       kappa_FAI=c(.x$kappa[,,1]),
                                       kappa_N=c(.x$kappa[,,2]))) %>%
        filter(month != 6 & year > pars.sens$tskip) %>%
        mutate(kappa=pmax(kappa_FAI, kappa_N)) %>%
        group_by(sim, month) %>%
        summarise(PAR_atDepth=first(PAR_atDepth), K_N=first(K_N), K_FAI=first(K_FAI),
                  FAI_mn=mean(FAI), FAI_md=median(FAI), FAI_sd=sd(FAI),
                  biomass_mn=mean(biomass), biomass_md=median(biomass), biomass_sd=sd(biomass),
                  p_kFAI=mean(kappa_FAI==kappa), p_kN=mean(kappa_N==kappa),
                  kappa_mn=mean(kappa), kappa_md=median(kappa), kappa_sd=sd(kappa),
                  kFAI_mn=mean(kappa_FAI), kFAI_md=median(kappa_FAI), kFAI_sd=sd(kappa_FAI),
                  kN_mn=mean(kappa_N), kN_md=median(kappa_N), kN_sd=sd(kappa_N)) %>%
        mutate(parDraw=k) 
    }
    pop.depth[[j]] <- do.call(rbind, pop.ls) %>% mutate(depth=pars.sens$depths[j])
    mass.depth[[j]] <- do.call(rbind, mass.ls) %>% mutate(depth=pars.sens$depths[j])
  }
  sim.info <- glue("{str_pad(x,4,'left','0')}_{gridRes}")
  saveRDS(do.call(rbind, pop.depth) %>% mutate(id=x), 
          glue("out{sep}sensitivity{sep}pop_{sim.info}.rds"))
  saveRDS(do.call(rbind, mass.depth) %>% mutate(id=x), 
          glue("out{sep}sensitivity{sep}mass_{sim.info}.rds"))
  return(x)
}








