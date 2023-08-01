# KELPER
# Helper functions
# Tim Szewczyk


# Miscellaneous helper functions associated with the KELPER project



# data preparation --------------------------------------------------------


#' Compile datasets from digitized csvs
#'
#' @param data.dir Directory with digitized csv files
#' @param supp.f Directory with collaborator files
#'
#' @return
#' @export
#'
#' @examples
compileDatasets <- function(data.dir, supp.f=NULL) {
  
  library(tidyverse)
  
  if(!is.null(supp.f)) {
    # read in field data -- requires different formatting
    library(readxl)
    supp.sheets <- excel_sheets(supp.f)[-1]
    supp.sites <- read_xlsx(supp.f, "siteLocations", skip=1) %>%
      select(Location, lat, lon) %>%
      rename(location=Location)
    supp.raw <- map(supp.sheets, ~read_xlsx(supp.f, .x, skip=1)) %>% 
      setNames(supp.sheets)
    supp.ls <- list()
    supp.raw$allometry <- supp.raw$allometry %>%
      mutate(across(contains("length"), ~.x*10),
             across(contains("width"), ~.x*10)) %>%
      rename(depth=Depth, location=Location)
    supp.ls$lengthStipe_weightStipe <- supp.raw$allometry %>% 
      select(lengthStipe, weightStipe, depth, location)
    supp.ls$lengthStipe_weightFrond <- supp.raw$allometry %>% 
      select(lengthStipe, weightFrond, depth, location)
    supp.ls$weightStipe_weightFrond <- supp.raw$allometry %>% 
      select(weightStipe, weightFrond, depth, location)
    supp.ls$lengthStipe_weightTotal <- supp.raw$allometry %>% 
      mutate(weightTotal=weightHoldfast + weightStipe + weightFrond) %>%
      select(lengthStipe, weightTotal, depth, location)
    supp.ls$depth_N <- supp.raw$abundStage %>%
      rename(depth=Depth, location=Location, N_canopy=CF_N, N_subcanopy=SC_N) %>%
      mutate(N_recruits=DJ_N+UJ_N) %>%
      select(depth, starts_with("N_"), location) %>%
      bind_rows(supp.raw$densityCanopy %>% 
                  rename(depth=Depth, location=Location, N_canopy=CanopyDensity) %>%
                  select(location, N_canopy, depth)) %>%
      rename(NperSqM=N_canopy) %>%
      mutate(data_id=row_number())
    supp.ls$depth_prCover <- supp.raw$abundStage %>%
      rename(depth=Depth, location=Location, pr_canopy=CF_pct, pr_subcanopy=SC_pct) %>%
      mutate(pr_recruits=DJ_N+UJ_N) %>%
      select(depth, "pr_canopy", location, Rep) %>%
      bind_rows(supp.raw$densityCanopy %>% 
                  rename(depth=Depth, location=Location, 
                         pr_canopy=CanopyPctCover,
                         Rep=QuadratRep) %>%
                  select(location, Rep, pr_canopy, depth)) %>%
      mutate(across(starts_with("pr_"), ~.x/100),
             data_id=row_number())
    supp.ls <- map(supp.ls, ~left_join(.x, supp.sites, by="location")) %>%
      setNames(., paste0("Smale_Moore_", 1:length(.))) %>%
      imap(~.x %>% mutate(reference=.y))
    
  }
  
  site.i <- read_csv(paste0(data.dir, "sitesDigitized.csv"), col_select=1:3, show_col_types=F)
  digitizedMetadata <- read_csv(paste0(data.dir, "metadata.csv"), show_col_types=F) %>%
    left_join(., site.i, by="location")
  fig.f <- unique(str_remove(dir(data.dir, "_fig.*csv"), "[a-z]?.csv"))
  tab.f <- unique(str_remove(dir(data.dir, "_table.*csv"), "[a-z]?.csv"))
  
  digitized <- c(fig.f, tab.f) %>% 
    map(~dir(data.dir, paste0(.x, "[a-z]?.csv")) %>%
          map_dfr(~read_csv(paste0(data.dir, .x), show_col_types=F) %>%
                    mutate(reference=str_remove(.x, ".csv"))) %>%
          left_join(., digitizedMetadata, by="reference", suffix=c("", ".meta")) %>%
          select(-contains(".meta"))) %>%
    setNames(c(fig.f, tab.f))
  
  
  
  fig_var <- imap_dfr(c(digitized, supp.ls), 
                      ~tibble(V1=names(.x)[1], V2=names(.x)[2], reference=.y)) %>%
    arrange(V1, V2) %>%
    mutate(vars=paste0(V1, "_", V2))
  
  unique_comparisons <- unique(fig_var$vars)
  varCombined.ls <- vector("list", length(unique_comparisons)) %>% 
    setNames(unique_comparisons)
  for(i in seq_along(varCombined.ls)) {
    refs <- fig_var$reference[fig_var$vars==unique_comparisons[i]]
    varCombined.ls[[i]] <- do.call('bind_rows', c(digitized, supp.ls)[refs])
  }
  return(varCombined.ls)
}








#' Load covariates summarised to mn, sd across years
#'
#' @param gis.dir GIS directory
#' @param bbox Bounding box
#' @param loadFile .rds file to load; takes priority if !NULL 
#' @param saveFile .rds file to save
#'
#' @return
#' @export
#'
#' @examples
loadCovariates <- function(gis.dir=NULL, bbox=NULL, loadFile=NULL, saveFile=NULL) {
  if(is.null(loadFile)) {
    rast.proj <- dir(paste0(gis.dir, "climate"), "MODISA_L3m_SST.*11W", full.names=T)[1] %>%
      raster() %>% crs
    covars.ls <- list(
      fetch=dir(paste0(gis.dir, "fetch"), "log10_UK200m_depth40_0na.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      logSlope=dir(paste0(gis.dir, "bathymetry"), "UK_log10slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      slope=dir(paste0(gis.dir, "bathymetry"), "UK_slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      sstDayGrow_mn=dir(paste0(gis.dir, "climate"), 
                        "MODISA_L3m_SST.*11W.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      sstDayGrow_sd=dir(paste0(gis.dir, "climate"), 
                        "MODISA_L3m_SST.*11W_49N.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      KD_mn=dir(paste0(gis.dir, "attenuation"), 
                "MODISA_L3m_KD.*11W_49N.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      KD_sd=dir(paste0(gis.dir, "attenuation"), 
                "MODISA_L3m_KD.*11W_49N.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      PAR_mn=dir(paste0(gis.dir, "light"), 
                 "MODISA_L3m_PAR.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(mean, na.rm=T),
      PAR_sd=dir(paste0(gis.dir, "light"), 
                 "MODISA_L3m_PAR.*tif$", full.names=T) %>%
        map(raster) %>% stack %>% crop(bbox) %>% calc(sd, na.rm=T),
      irrad_monthly=dir(paste0(gis.dir, "light"), 
                        "POWER_Regional_monthly", full.names=T) %>%
        map_dfr(~read_csv(.x, skip=9, show_col_types=F)) %>% 
        filter(JAN != -999) %>% select(-ANN) %>%
        group_by(PARAMETER, YEAR, LAT, LON) %>%
        summarise(across(everything(), mean)) %>%
        group_by(PARAMETER, YEAR) %>% mutate(id=row_number()) %>% ungroup %>%
        pivot_longer(5:16, names_to="MONTH", values_to="PAR") %>%
        st_as_sf(coords=c("LON", "LAT"), crs=4326)
    )
    
    irrad.grid <- (st_bbox(covars.ls$irrad_monthly) + .25*c(-1,-1,1,1)) %>%
      st_make_grid(cellsize=c(0.5, 0.5)) %>% st_sf(id=1:length(.))
    
    covars.ls$irrad_growing.month <- covars.ls$irrad_monthly %>%
      filter(MONTH %in% c("JAN", "FEB", "MAR", "APR", "MAY", "JUN")) %>%
      group_by(YEAR, id) %>%
      summarise(PAR=mean(PAR)) %>%
      group_by(id) %>%
      summarise(mnPAR=mean(PAR), sdPAR=sd(PAR)) %>%
      ungroup() %>% select(-id) %>%
      st_join(irrad.grid %>% select(-id), .)
    covars.ls <- covars.ls[-which(names(covars.ls)=="irrad_monthly")]
    
    if(!is.null(saveFile)) {
      saveRDS(covars.ls, saveFile)
    }
  } else {
    covars.ls <- readRDS(loadFile)
  }
  return(covars.ls)
}



#' Load covariates unsummarised across years
#'
#' @param gis.dir GIS directory
#' @param bbox Bounding box
#' @param loadFile .rds file to load; takes priority if !NULL 
#' @param saveFile .rds file to save
#'
#' @return
#' @export
#'
#' @examples
loadCovariates_full <- function(gis.dir=NULL, bbox=NULL, loadFile=NULL, saveFile=NULL) {
  library(stars)
  if(is.null(loadFile)) {
    rast.proj <- CRS("+init=epsg:4326")
    covars.ls <- list(
      fetch=dir(paste0(gis.dir, "fetch"), "log10_UK200m_depth40_0na.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      logSlope=dir(paste0(gis.dir, "bathymetry"), "UK_log10slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      slope=dir(paste0(gis.dir, "bathymetry"), "UK_slope40.tif$", full.names=T) %>%
        raster %>% projectRaster(., crs=rast.proj) %>% crop(bbox),
      sstDayGrow=dir(paste0(gis.dir, "climate"), 
                     "MODISA_L3m_SST.*11W_49N.*tif$", full.names=T) %>%
        setNames(., str_sub(., 84, 87)) %>%
        map(raster) %>% stack() %>% projectRaster(., crs=rast.proj) %>% crop(bbox) %>% 
        st_as_stars() %>% st_as_sf() %>% 
        pivot_longer(starts_with("X"), names_to="YEAR", values_to="SST") %>%
        mutate(YEAR=str_sub(YEAR, 2, -1)),
      KD_mn=dir(paste0(gis.dir, "attenuation"), 
                "MODISA_L3m_KD.*11W_49N.*tif$", full.names=T) %>%
        setNames(., str_sub(., 84, 87)) %>%
        map(raster) %>% stack() %>% projectRaster(., crs=rast.proj) %>% crop(bbox) %>% 
        st_as_stars() %>% st_as_sf() %>% 
        pivot_longer(starts_with("X"), names_to="YEAR", values_to="KD") %>%
        mutate(YEAR=str_sub(YEAR, 2, -1)),
      PAR_mn=dir(paste0(gis.dir, "light"), 
                 "MODISA_L3m_PAR.*tif$", full.names=T)  %>%
        setNames(., str_sub(., 76, 79)) %>%
        map(raster) %>% stack() %>% projectRaster(., crs=rast.proj) %>% crop(bbox) %>% 
        st_as_stars() %>% st_as_sf() %>% 
        pivot_longer(starts_with("X"), names_to="YEAR", values_to="PAR") %>%
        mutate(YEAR=str_sub(YEAR, 2, -1)),
      irrad_monthly=dir(paste0(gis.dir, "light"), 
                        "POWER_Regional_monthly", full.names=T) %>%
        map_dfr(~read_csv(.x, skip=9, show_col_types=F)) %>% 
        filter(JAN != -999) %>% select(-ANN) %>%
        group_by(PARAMETER, YEAR, LAT, LON) %>%
        summarise(across(everything(), mean)) %>%
        group_by(PARAMETER, YEAR) %>% mutate(id=row_number()) %>% ungroup %>%
        pivot_longer(5:16, names_to="MONTH", values_to="PAR") %>%
        st_as_sf(coords=c("LON", "LAT"), crs=rast.proj)
    )
    
    irrad.grid <- (st_bbox(covars.ls$irrad_monthly) + .25*c(-1,-1,1,1)) %>%
      st_make_grid(cellsize=c(0.5, 0.5)) %>% st_sf(id=1:length(.), crs=rast.proj)
    
    covars.ls$irrad_growing.month <- covars.ls$irrad_monthly %>%
      filter(MONTH %in% c("JAN", "FEB", "MAR", "APR", "MAY", "JUN")) %>%
      group_by(YEAR, id) %>%
      summarise(PAR=mean(PAR)) %>%
      ungroup() %>% select(-id) %>%
      st_join(irrad.grid %>% select(-id), .) %>%
      filter(!is.na(YEAR))
    covars.ls <- covars.ls[-which(names(covars.ls)=="irrad_monthly")]
    
    if(!is.null(saveFile)) {
      saveRDS(covars.ls, saveFile)
    }
  } else {
    covars.ls <- readRDS(loadFile)
  }
  return(covars.ls)
}







#' Extract covariate values to dataset locations
#'
#' @param data.ls List loaded with compileDatasets()
#' @param covars.ls Covariates loaded with loadCovariates()
#' @param PAR_datasource PAR source (MODIS vs POWER)
#' @param grid.sf Grid sf object
#'
#' @return
#' @export
#'
#' @examples
extractCovarsToDatasets <- function(data.ls, covars.ls=NULL, PAR_datasource, grid.sf=NULL) {
  if(is.null(grid.sf)) {
    fetch_thresh <- 4.02 # bad practice, but based on Pedersen sites (low, med+high)
    for(i in seq_along(data.ls)) {
      if(any(!is.na(data.ls[[i]]$lat))) {
        i_sf <- data.ls[[i]] %>% 
          filter(!is.na(lat) & !is.na(lon)) %>%
          st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
          mutate(sstDay_mn=raster::extract(covars.ls$sstDayGrow_mn, ., buffer=8e3,
                                           fun=mean, small=T),
                 sstDay_sd=raster::extract(covars.ls$sstDayGrow_sd, ., buffer=8e3, 
                                           fun=mean, small=T),
                 KD_mn=raster::extract(covars.ls$KD_mn, ., buffer=8e3, 
                                       fun=mean, small=T),
                 KD_sd=raster::extract(covars.ls$KD_sd, ., buffer=8e3, 
                                       fun=mean, small=T),
                 PAR_mn=raster::extract(covars.ls$PAR_mn, ., buffer=8e3, 
                                        fun=mean, small=T),
                 PAR_sd=raster::extract(covars.ls$PAR_sd, ., buffer=8e3, 
                                        fun=mean, small=T),
                 PAR_POWER=st_join(., 
                                   covars.ls$irrad_growing.month, 
                                   left=T, join=st_nearest_feature)$mnPAR,
                 logSlope=raster::extract(covars.ls$logSlope, ., buffer=8e3, 
                                          fun=mean, small=T),
                 slope=raster::extract(covars.ls$slope, ., buffer=8e3, 
                                       fun=mean, small=T),
                 fetch=raster::extract(covars.ls$fetch, ., buffer=8e3, 
                                       fun=mean, small=T),
                 fetchCat=(fetch < fetch_tresh)+1)
        data.ls[[i]] <- left_join(data.ls[[i]], st_drop_geometry(i_sf))
        if(any(!is.na(data.ls[[i]]$depth) & 
               !is.na(data.ls[[i]]$KD_mn) & 
               !is.na(data.ls[[i]]$PAR_mn))) {
          data.ls[[i]] <- data.ls[[i]] %>%
            mutate(PAR_atDepth_POWER=PAR_POWER * exp(-KD_mn * depth),
                   PAR_atDepth_MODIS=PAR_mn * exp(-KD_mn * depth))
          data.ls[[i]]$PAR_atDepth <- data.ls[[i]][[paste0("PAR_atDepth_", PAR_datasource)]]
        }
      }
    }
  } else {
    for(i in seq_along(data.ls)) {
      if(any(!is.na(data.ls[[i]]$lat))) {
        i_sf <- data.ls[[i]] %>% 
          filter(!is.na(lat) & !is.na(lon)) %>%
          st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
          st_join(., grid.sf, left=T)
        data.ls[[i]] <- left_join(data.ls[[i]], st_drop_geometry(i_sf))
        if(any(!is.na(data.ls[[i]]$depth) & 
               !is.na(data.ls[[i]]$KD) & 
               !is.na(data.ls[[i]]$PAR))) {
          data.ls[[i]] <- data.ls[[i]] %>%
            mutate(PAR_atDepth=PAR * exp(-KD * depth))
        }
      }
    }
  }
  
  return(data.ls)
}




#' Extract covariate values to point locations
#'
#' @param site.i Dataframe with site lat, lon
#' @param covars.ls Covariates loaded with loadCovariates()
#' @param PAR_datasource PAR source (MODIS vs POWER)
#'
#' @return
#' @export
#'
#' @examples
extractCovarsToPts <- function(site.i, covars.ls, PAR_datasource) {
  fetch_thresh <- 4.02 # bad practice, but based on Pedersen sites (low, med+high)
  site.sf <- site.i %>% 
    filter(!is.na(lat) & !is.na(lon)) %>%
    st_as_sf(coords=c("lon", "lat"), crs=4326) %>%
    mutate(sstDay_mn=raster::extract(covars.ls$sstDayGrow_mn, ., 
                                     fun=mean, small=T, method="bilinear"),
           sstDay_sd=raster::extract(covars.ls$sstDayGrow_sd, ., 
                                     fun=mean, small=T, method="bilinear"),
           KD_mn=raster::extract(covars.ls$KD_mn, ., 
                                 fun=mean, small=T, method="bilinear"),
           KD_sd=raster::extract(covars.ls$KD_sd, ., 
                                 fun=mean, small=T, method="bilinear"),
           PAR_mn=raster::extract(covars.ls$PAR_mn, ., 
                                  fun=mean, small=T, method="bilinear"),
           PAR_sd=raster::extract(covars.ls$PAR_sd, ., 
                                  fun=mean, small=T, method="bilinear"),
           PAR_POWER=st_join(., 
                             covars.ls$irrad_growing.month, 
                             left=T, join=st_nearest_feature)$mnPAR,
           logSlope=raster::extract(covars.ls$logSlope, ., 
                                    fun=mean, small=T, method="bilinear"),
           slope=raster::extract(covars.ls$slope, ., 
                                 fun=mean, small=T, method="bilinear"),
           fetch=raster::extract(covars.ls$fetch, ., 
                                 fun=mean, small=T, method="bilinear"),
           fetchCat=(fetch < fetch_tresh)+1)
  site.sf$PAR_surface <- site.sf[[ifelse(PAR_datasource=="MODIS", "PAR_mn", "PAR_POWER")]]
  return(site.sf)
}





#' Extract covariate values to grid cells
#'
#' @param grid.domain Grid sf file
#' @param covars.ls Covariates loaded with loadCovariates()
#' @param PAR_datasource PAR source (MODIS vs POWER)
#'
#' @return
#' @export
#'
#' @examples
extractCovarsToGrid <- function(grid.domain, covars.ls, PAR_datasource) {
  fetch_thresh <- 4.02 # bad practice, but based on Pedersen sites (low, med+high)
  grid.domain <- grid.domain %>% 
    mutate(sstDay_mn=raster::extract(covars.ls$sstDayGrow_mn, grid.domain, fun=mean, na.rm=T),
           sstDay_sd=raster::extract(covars.ls$sstDayGrow_sd, grid.domain, fun=mean, na.rm=T),
           KD_mn=raster::extract(covars.ls$KD_mn, grid.domain, fun=mean, na.rm=T),
           KD_sd=raster::extract(covars.ls$KD_sd, grid.domain, fun=mean, na.rm=T),
           PAR_mn=raster::extract(covars.ls$PAR_mn, grid.domain, fun=mean, na.rm=T),
           PAR_sd=raster::extract(covars.ls$PAR_sd, grid.domain, fun=mean, na.rm=T),
           logSlope=raster::extract(covars.ls$logSlope, grid.domain, fun=mean, na.rm=T),
           slope=raster::extract(covars.ls$slope, grid.domain, fun=mean, na.rm=T),
           fetch=raster::extract(covars.ls$fetch, grid.domain, fun=mean, na.rm=T)
    ) %>%
    st_join(., 
            covars.ls$irrad_growing.month %>% select(mnPAR) %>% rename(PAR_POWER=mnPAR), 
            left=T) %>%
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    mutate(fetchCat=(fetch < fetch_tresh)+1)
  
  grid.domain$PAR_surface <- grid.domain[[ifelse(PAR_datasource=="MODIS", "PAR_mn", "PAR_POWER")]]
  return(grid.domain)
}









# misc. processing --------------------------------------------------------



#' Lag multiple variables at once
#'
#' From https://stackoverflow.com/questions/55814028/multiple-lags-with-dplyr
#'
#' @param data Dataframe
#' @param ... Unquoted variable names to lag
#' @param n Number of lags
#'
#' @return
#' @export
#'
#' @examples
multijetlag <- function(data, ..., n=10){
  library(rlang)
  variable <- enquos(...)
  
  indices <- seq_len(n)
  combos <- crossing(indices, var =as.list(variable))
  
  quosures <- map2(combos$indices, combos$var,
                   ~quo(lag(!!.y, !!.x)) ) %>% 
    set_names(paste("lag", combos$indices, map_chr(combos$var, quo_text), sep = "_"))
  mutate( data, !!!quosures )
  
}





#' Generate natural scale prediction from fitted model
#'
#' @param mod Fitted model
#' @param ndraws Number of draws from the posterior distribution to use
#' @param new.df Dataframe with covariates for input
#' @param scale.df Dataframe with mean, sd to de-center and de-scale
#' @param y_var Response variable
#'
#' @return
#' @export
#'
#' @examples
getPrediction <- function(mod, ndraws, new.df, scale.df, y_var) {
  # center and scale predictors to match regression inputs
  for(i in 1:ncol(new.df)) {
    if(is.numeric(new.df[[i]])) {
      scale.row <- which(scale.df$Par==names(new.df)[i])
      new.df[[i]] <- (new.df[[i]] - scale.df$ctr[scale.row])/scale.df$scl[scale.row]
    }
  }
  
  pred <- colMeans(posterior_epred(mod, newdata=new.df, ndraws=ndraws, re.form=NA))
  
  # de-center and de-scale predictions to natural scale
  if(y_var != "N") {
    scale.row <- which(scale.df$Par==y_var)
    pred <- pred * scale.df$scl[scale.row] + scale.df$ctr[scale.row] 
  }
  return(pred)
}







getHDIs <- function(data.df, response) {
  library(HDInterval)
  data.df %>% 
    summarise(mn=mean({{response}}),
              md=median({{response}}),
              lo1=hdi({{response}}, 0.2)[1],
              hi1=hdi({{response}}, 0.2)[2],
              lo2=hdi({{response}}, 0.3)[1],
              hi2=hdi({{response}}, 0.3)[2],
              lo3=hdi({{response}}, 0.4)[1],
              hi3=hdi({{response}}, 0.4)[2],
              lo4=hdi({{response}}, 0.5)[1],
              hi4=hdi({{response}}, 0.5)[2],
              lo5=hdi({{response}}, 0.6)[1],
              hi5=hdi({{response}}, 0.6)[2],
              lo6=hdi({{response}}, 0.7)[1],
              hi6=hdi({{response}}, 0.7)[2],
              lo7=hdi({{response}}, 0.8)[1],
              hi7=hdi({{response}}, 0.8)[2],
              lo8=hdi({{response}}, 0.9)[1],
              hi8=hdi({{response}}, 0.9)[2],
              lo9=hdi({{response}}, 0.95)[1],
              hi9=hdi({{response}}, 0.95)[2])
} 

getQuantiles <- function(data.df, response) {
  data.df %>% 
    summarise(mn=mean({{response}}),
              md=median({{response}}),lo1=quantile({{response}}, probs=0.45),
              hi1=quantile({{response}}, probs=0.55),
              lo2=quantile({{response}}, probs=0.4),
              hi2=quantile({{response}}, probs=0.6),
              lo3=quantile({{response}}, probs=0.35),
              hi3=quantile({{response}}, probs=0.65),
              lo4=quantile({{response}}, probs=0.3),
              hi4=quantile({{response}}, probs=0.7),
              lo5=quantile({{response}}, probs=0.25),
              hi5=quantile({{response}}, probs=0.75),
              lo6=quantile({{response}}, probs=0.2),
              hi6=quantile({{response}}, probs=0.8),
              lo7=quantile({{response}}, probs=0.15),
              hi7=quantile({{response}}, probs=0.85),
              lo8=quantile({{response}}, probs=0.1),
              hi8=quantile({{response}}, probs=0.9),
              lo9=quantile({{response}}, probs=0.05),
              hi9=quantile({{response}}, probs=0.95))
} 









# simulation preparation -----------------------------------------------


#' Generate list of parameters for simulations
#' 
#' All rates must already be scaled to the timeframe of the transition matrix (e.g., growth rate for the season rather than daily).
#'
#' @param path Path to .rds file with stored parameters; if \code{NULL}, creates list instead
#' @param override If `path` is not NULL, named list of parameters to override (e.g., \code{list(N_stages=2)})
#' @param N_stages Number of stages
#' @param N_seasons Number of seasons per year
#' @param tmax Number of years
#' @param growthRateStipeMax Vector of maximum stipe growth rates during spring
#' @param growthRateFrond Vector of maximum frond area growth rates during spring 
#' @param frondAreaMax Carrying capacity for Frond Area Index
#' @param growthRateStipeDensityShape Shape parameter for stipe growth density dependence
#' @param sizeClassLimits $\delta$ Size limits between stages
#' @param survRate $\s$ Survival rate
#' @param settlementRateBg Background recruit settlement rate
#' @param lossRate Vector of frond loss rates during winter
#' @param prFullHarvest Proportion of population to be harvested in full
#' @param freqHarvest Harvest frequency (years)
#' @param harvestTarget Kelp part for biomass calculation (all, stipe, frond)
#' @param extraPars Any additional parameters to include (e.g., regressions)
#'
#' @return
#' @export
#'
#' @examples
setParameters <- function(path=NULL,
                          override=NULL,
                          tmax=20,
                          growthRateStipeMax=cbind(c(194, 195, 58),
                                                   c(23, 21, 10)), # Kain 1976
                          growthRateFrond=cbind(c(1787, 2299, 3979)/1e4,
                                                c(179, 397, 417)/1e4), # m^2/plant/growing season
                          frondAreaMax=5500/1e4,
                          growthRateDensityShape=2,
                          sizeClassLimits=(1000 * (0:6)/6)[c(1,2,5,7)],
                          survRate=c(0.4, 0.7, 0.9),
                          settlementRate=400,
                          lossRate=0.2345,
                          prFullHarvest=0,
                          freqHarvest=1e5,
                          harvestTarget=c("all", "stipe", "frond")[1],
                          extraPars=NULL,
                          stochParams=FALSE,
                          stormIntensity=NULL
) {
  
  if(!is.null(path)) {
    pars <- readRDS(path)
    if(!is.null(override)) {
      for(i in seq_along(override)) {
        pars[names(override)[i]] <- override[i]
      }
    }
  } else {
    pars <- list(
      
      # simulation 
      tmax=tmax,
      stochParams=stochParams,
      stormIntensity=stormIntensity,
      
      # growth
      growthRateStipeMax=growthRateStipeMax,
      growthRateFrond=growthRateFrond,
      frondAreaMax=frondAreaMax,
      growthRateDensityShape=growthRateDensityShape,
      sizeClassLimits=sizeClassLimits,
      sizeClassMdpts=(sizeClassLimits+lag(sizeClassLimits))[-1]/2,
      
      # survival
      survRate=survRate,
      
      # settlement
      settlementRate=settlementRate,
      
      # loss
      lossRate=lossRate,
      
      # harvest
      prFullHarvest=prFullHarvest,
      freqHarvest=freqHarvest,
      harvestTarget=harvestTarget
    ) 
  }
  
  if(!is.null(extraPars)) {
    for(i in seq_along(extraPars)) {
      pars[names(extraPars)[i]] <- extraPars[i]
    }
  }
  
  return(pars)
}






#' Simulate landscape for a covariate
#'
#' @param grid.sf Grid sf object
#' @param var.sf Covariate sf object (element of loadCovariates_full() output)
#' @param nYr Number of years to simulate
#' @param colName Name of column within var.sf with focal covariate
#'
#' @return
#' @export
#'
#' @examples
simulateLandscape <- function(grid.sf, var.sf, nYr, colName) {
  library(mvtnorm)
  # pair var with grid, then extract to matrix
  var.mx <- grid.sf %>% select(id) %>%
    st_join(., 
            var.sf %>%
              pivot_wider(names_from="YEAR", values_from=colName, names_prefix="y_") %>%
              st_as_sf(),
            left=T) %>%
    group_by(id) %>% summarise(across(everything(), mean, na.rm=T)) %>% ungroup %>%
    st_drop_geometry() %>% select(-id) %>% as.matrix()
  # generate simulated values
  na.ids <- which(is.na(var.mx), arr.ind=T)
  if(nrow(na.ids)>0) {
    for(i in 1:nrow(na.ids)) {
      var.mx[na.ids[i,1], na.ids[i,2]] <- mean(var.mx[na.ids[i,1],], na.rm=T) +
        rnorm(1, 0, 0.01)
    }
  }
  if(colName=="KD") var.mx <- log(var.mx)
  var.sim <- t(rmvnorm(nYr, rowMeans(var.mx), cov(t(var.mx))))
  colnames(var.sim) <- paste0(colName, "_", 1:ncol(var.sim))
  if(colName=="KD") var.sim <- exp(var.sim)
  var.sim[var.sim < 0] <- min(var.sf[[colName]], na.rm=T)
  # bind to grid
  grid.updated <- grid.sf %>% bind_cols(as_tibble(var.sim)) %>% 
    pivot_longer(starts_with(paste0(colName, "_")), 
                 names_to="year", values_to=colName) %>%
    mutate(year=as.numeric(str_sub(year, nchar(colName)+2, -1)))
  return(grid.updated)
}













# sensitivity analysis ----------------------------------------------------


#' Emulate global sensitivity analysis output
#'
#' Emulate the output from a global sensitivity analysis using boosted
#' regression trees with different interaction depths. Based on function
#' described in Prowse et al 2016. Writes files to brt.dir
#' @param sens.out Dataframe of the parameter sets and simulation summaries;
#'   \code{.$results} from \link{global_sensitivity}
#' @param params Vector of parameter names to include
#' @param n.cores \code{1} Number of cores for fitting subsample BRTs
#' @param n.sub \code{10} Number of subsamples for each emulation
#' @param td \code{c(1,3,5)} Vector of regression tree interaction depths to
#'   test
#' @param resp Which response summary to use (column name from \code{sens.out})
#' @param brt.dir \code{"out/brt/"} Directory to store boosted regression tree
#'   output
#' @return Success message
#' @keywords parameters, sensitivity, save, output
#' @export

emulate_sensitivity <- function(sens.out, params, resp, brt.dir, siminfo, 
                                n.sub=10, td=c(1,3,5)) {
  library(tidyverse)
  sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")
  if(!dir.exists(brt.dir)) dir.create(brt.dir)
  sub.prop <- seq(0.75, 1, length.out=n.sub)
  
  for(i in 1:n.sub) {
    # subset sensitivity results
    sub.samp <- sample_frac(sens.out, sub.prop[i]) %>% as.data.frame
    n <- nrow(sub.samp)
    
    # fit BRTs of different tree complexities for given response variable
    for(j in 1:length(td)) {
      td_j <- td[j]
      brt.fit <- dismo::gbm.step(sub.samp, gbm.x=params, gbm.y=resp, 
                                 max.trees=200000, n.folds=5, 
                                 family="gaussian", tree.complexity=td_j,
                                 bag.fraction=0.8, silent=T, plot.main=F)
      saveRDS(brt.fit, glue::glue("{brt.dir}{sep}{resp}_{siminfo}_td-{td_j}-{n}.rds"))
    }
  }
  return(paste("Finished emulations of", resp))
}









#' Summarize BRT emulations
#'
#' Summarize the output from global sensitivity analysis boosted regression
#' trees emulations. Based on function described in Prowse et al 2016. Reads BRT
#' output saved via emulate_sensitivity
#' @param resp Which response summary to use (column name from \code{sens.out})
#' @return List of two dataframes: relative influences, and cross-validation
#'   deviances, each across different subsample sizes and tree complexities.
#' @param brt.dir \code{"out/brt/"} Directory with stored boosted regression
#'   tree output
#' @keywords parameters, sensitivity, save, output
#' @export

emulation_summary <- function(resp, brt.dir, siminfo) {
  library(gbm); library(tidyverse)
  sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")
  f <- dir(brt.dir, paste0(resp, "_", siminfo))
  f.i <- str_split_fixed(f, "-", 3)
  id <- as.numeric(str_split_fixed(siminfo, "_", 4)[,1])
  month <- str_split_fixed(siminfo, "_", 4)[,3]
  depth <- str_split_fixed(siminfo, "_", 4)[,4]
  cvDev.df <- tibble(response=resp,
                     td=f.i[,2],
                     smp=as.numeric(str_remove(f.i[,3], ".rds")),
                     id=id,
                     month=month,
                     depth=depth,
                     Dev=NA)
  ri.ls <- vector("list", length(f))
  for(i in seq_along(f)) {
    brt <- readRDS(paste0(brt.dir, f[i]))
    # cross validation deviance
    cvDev.df$Dev[i] <- brt$cv.statistics$deviance.mean
    # relative influence
    ri.ls[[i]] <- as_tibble(brt$contributions) %>%
      mutate(td=cvDev.df$td[i], 
             smp=cvDev.df$smp[i], 
             response=cvDev.df$response[i])
  }
  ri.df <- bind_rows(ri.ls) %>% 
    group_by(response, td, smp, var) %>%
    summarise(rel.inf=sum(rel.inf)) %>% 
    ungroup %>% group_by(response, td, smp) %>%
    mutate(rel.inf=rel.inf/sum(rel.inf)) %>%
    ungroup %>% arrange(desc(td), desc(smp), desc(rel.inf)) %>%
    mutate(id=id,
           month=month,
           depth=depth)
  
  write_csv(cvDev.df, glue::glue("{brt.dir}{sep}summaries{sep}{resp}_cv_{siminfo}.csv"))
  write_csv(ri.df, glue::glue("{brt.dir}{sep}summaries{sep}{resp}_ri_{siminfo}.csv"))
  return()
}








#' Run Boosted Regression Trees
#' 
#' Wrapper for emulate_sensitivity() > emulation_summary()
#'
#' @param x Index for parLapply
#' @param pop.f Vector of files with population output
#' @param mass.f Vector of files with biomass output
#' @param parSets List of parameter values
#' @param grid.i Dataframe with grid information
#' @param meta.cols Deprecated
#' @param params Vector of parameters to evaluate
#' @param brt.dir BRT output directory
#' @param resp Vector of response variables to evaluate
#'
#' @return
#' @export
#'
#' @examples
runBRTs <- function(x, pop.f=NULL, mass.f, parSets, grid.i, meta.cols, params, 
                    brt.dir, resp=NULL) {
  library(glue); library(tidyverse)
  sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")
  siminfo <- str_remove(str_split_fixed(mass.f[x], "mass_", 2)[,2], ".rds")
  mass.i <- readRDS(mass.f[x]) %>%
    left_join(., parSets[[3]]) %>%
    left_join(., parSets[[filter(grid.i, id==.$id[1])$fetchCat]]) %>%
    group_by(depth) %>%
    group_split()
  depths <- map_dbl(mass.i, ~.x$depth[1])
  if(is.null(resp)) {
    resp <- c("biomass_mn", "biomass_sd")
  }
  for(j in 1:length(resp)) {
    for(k in 1:length(depths)) {
      emulate_sensitivity(mass.i[[k]] %>% filter(month==1), params,
                          resp[j], brt.dir, glue("{siminfo}_jan_{depths[k]}m"))
      emulate_sensitivity(mass.i[[k]] %>% filter(month==7), params,
                          resp[j], brt.dir, glue("{siminfo}_july_{depths[k]}m"))
      try(emulation_summary(resp[j], brt.dir, glue("{siminfo}_jan_{depths[k]}m")))
      try(emulation_summary(resp[j], brt.dir, glue("{siminfo}_july_{depths[k]}m")))
    }
  }
}




