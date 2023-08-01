# KELPER
# Sensitivity analysis
# Tim Szewczyk

# This script runs simulations in each grid cell across specified depths with a
# survival and/or loss rate dependent on annual storm intensity




# set up ------------------------------------------------------------------


# libraries and local functions
pkgs <- c("glue", "tidyverse", "sf", "brms", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# directories & files
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
out.dir <- glue("out{sep}storms{sep}")
par.dir <- glue("data{sep}")

# switches & settings
gridRemake <- F
rerun <- F
nCores <- 4
nSimGrid <- 500
set.seed(789)
options(mc.cores=nCores)
pars.sim <- list(depths=c(2, 5, 10, 15, 20), 
                 tskip=0,
                 nSim=1, 
                 stochParams=T,
                 landscape="static")

# load files
grid.sf <- st_read(glue("data{sep}grid_0.1_MODIS.gpkg")) %>% 
  mutate(fetchCat=pmin(fetchCat, 2))
grid.i <- grid.sf %>% st_drop_geometry() %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn) %>%
  select(id, SST, KD, PAR, fetch, fetchCat)
data.ls <- compileDatasets(data.dir, supp.f)
lm.fit <- readRDS(glue("data{sep}fits_0.1.rds"))
lm.mnsd <- readRDS(glue("data{sep}dfs_mn_sd_0.1.rds"))

# set sim length, storm intensity
pars.sim$tmax <- nrow(data.ls$year_stormIndex)
pars.sim$storms <- data.ls$year_stormIndex$stormIndex


# simulate landscapes if needed
if(pars.sim$landscape=="dynamic") {
  if(gridRemake) {
    covars.full <- loadCovariates_full(loadFile=glue("data{sep}covarFull_ls_reproj.rds"))
    for(i in 1:nSimGrid) {
      grid.sf %>% select(id) %>%
        simulateLandscape(covars.full$sstDayGrow, pars.sim$tmax, "SST") %>%
        full_join(., 
                  simulateLandscape(grid.sf %>% select(id), 
                                    covars.full$KD_mn, pars.sim$tmax, "KD") %>%
                    st_drop_geometry()) %>%
        full_join(., 
                  simulateLandscape(grid.sf %>% select(id), 
                                    covars.full$PAR_mn, pars.sim$tmax, "PAR") %>%
                    st_drop_geometry()) %>%
        left_join(., 
                  grid.sf %>% st_drop_geometry() %>% select(id, fetch, fetchCat)) %>%
        st_drop_geometry() %>%
        saveRDS(glue("data{sep}gridSim{sep}gridSim_0.1_{str_pad(i,3,'left','0')}.rds"))
    }
    rm(covars.full)
  }
}






# run simulations ---------------------------------------------------------

if(rerun) {
  cl <- makeCluster(nCores, outfile=glue("temp{sep}sim_out.txt"))
  
  for(i in 1:nSimGrid) {
    if(pars.sim$landscape=="dynamic") {
      grid.sim <- readRDS(glue("data{sep}gridSim{sep}gridSim_0.1_{str_pad(i,3,'left','0')}.rds")) %>% 
        mutate(fetchCat=pmin(fetchCat, 2))
    } else {
      grid.sim <- NULL
    }
    obj.exclude <- c("data.ls", "grid.sf")
    obj.include <- ls()
    
    clusterExport(cl, obj.include[-match(obj.exclude, obj.include)])
    Sys.sleep(1)
    cat("Starting runs: grid", i, "of", nSimGrid, "\n")
    out.ls <- parLapply(cl, X=1:nrow(grid.i), fun=simDepthsWithinCell,
                        grid.i=grid.i, grid.sim=grid.sim, grid.id=i, gridRes=gridRes,
                        pars.sim=pars.sim, par.dir=par.dir,
                        lm.fit=lm.fit, lm.mnsd=lm.mnsd)
  }
  
  stopCluster(cl)
}

