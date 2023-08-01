# KELPER
# Sensitivity analysis
# Tim Szewczyk

# This script runs simulations in each grid cell across specified depths for
# combinations of parameter values




# set up ------------------------------------------------------------------


# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms", "parallel")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")


# directories & files
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
sens.dir <- glue("out{sep}sensitivity{sep}")
parSets.f <- "000_paramSets_exposure_"

# switches & settings
rerun <- T
reanalyse <- T
nCores <- 4
options(mc.cores=nCores)
pars.sens <- list(nParDraws=1e3,
                  depths=c(2, 5, 10, 15),
                  tmax=30,
                  tskip=10,
                  nSim=1,
                  stochParams=F, 
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
surv.df <- read_csv(glue("data{sep}par_survival.csv"))
recruitment.df <- read_csv(glue("data{sep}par_recruitment.csv"))
erosion.df <- read_csv(glue("data{sep}par_erosion.csv"))
stipeGrow.df <- read_csv(glue("data{sep}par_growStipe.csv"))
frondGrow.df <- read_csv(glue("data{sep}par_growFrond.csv"))








# assign parameter values -------------------------------------------------


# parameter ranges: 90% CIs
p <- c(0.05, 0.95)
par.rng <- list(surv=surv.df %>%
                  mutate(valMean=mn,
                         valPrec=prec,
                         valMin=qbeta(p[1], shp1, shp2),
                         valMax=qbeta(p[2], shp1, shp2)) %>%
                  select(stage, exposure, valMean, valPrec, valMin, valMax, shp1, shp2),
                settlement=recruitment.df %>%
                  mutate(valMean=mn,
                         valSD=sd,
                         valMin=truncnorm::qtruncnorm(p[1], 0, Inf, valMean, valSD),
                         valMax=truncnorm::qtruncnorm(p[2], 0, Inf, valMean, valSD)) %>%
                  select(exposure, valMean, valSD, valMin, valMax),
                growStipe=stipeGrow.df %>%
                  mutate(valMean=mn, 
                         valSD=sd,
                         valMin=qnorm(p[1], valMean, valSD),
                         valMax=qnorm(p[2], valMean, valSD)) %>%
                  select(stage, exposure, valMean, valSD, valMin, valMax),
                growFrond=frondGrow.df %>%
                  mutate(valMean=mn, 
                         valSD=sd,
                         valMin=qnorm(p[1], valMean, valSD),
                         valMax=qnorm(p[2], valMean, valSD)) %>%
                  select(stage, exposure, valMean, valSD, valMin, valMax),
                loss=erosion.df %>%
                  mutate(valMean=mn,
                         valPrec=prec,
                         valMin=qbeta(p[1], shp1, shp2),
                         valMax=qbeta(p[2], shp1, shp2)) %>%
                  select(stage, exposure, valMean, valPrec, valMin, valMax, shp1, shp2),
                densityEffShape=tibble(valMean=1, valMin=0.5, valMax=1.5)) %>%
  imap_dfr(., ~.x %>% mutate(param=.y)) %>%
  select(param, stage, exposure, valMean, valSD, valPrec, valMin, valMax, shp1, shp2)
saveRDS(par.rng, glue("{sens.dir}000_parameter_ranges.rds"))

if(rerun) {
  parSets <- map_dfr(1:pars.sens$nParDraws, 
                     ~par.rng %>% mutate(parDraw=.x, 
                                         val=runif(n(), valMin, valMax))) %>%
    select(parDraw, param, stage, exposure, val) %>%
    mutate(exposureNum=case_when(exposure=="low" ~ 1, exposure=="high" ~ 2)) %>%
    group_by(exposureNum) %>% group_split()
  walk2(parSets, c(1:2, "NA"), ~write_csv(.x, glue("{sens.dir}{parSets.f}{.y}.csv")))
}










# run simulations ---------------------------------------------------------


if(rerun) {
  cl <- makeCluster(nCores, outfile=glue("temp{sep}sensitivity_out.txt"))
  obj.exclude <- c("data.ls", "grid.sf", "surv.df", "fecund.df")
  obj.include <- ls()
  
  cat("Exporting cluster.\n")
  clusterExport(cl, obj.include[-match(obj.exclude, obj.include)])
  Sys.sleep(5)
  cat("Exported. Starting parallel runs.\n")
  out.ls <- parLapply(cl, X=1:nrow(grid.i), fun=simSensitivityDepthsWithinCell,
                      grid.i=grid.i, gridRes=0.1, pars.sens=pars.sens, 
                      lm.fit=lm.fit, lm.mnsd=lm.mnsd, parSets=parSets)
  stopCluster(cl)
  cat("Finished simulations.\n")
}










# run BRTs ----------------------------------------------------------------


if(reanalyse) {
  meta.cols <- c("sim", "month", "stage", "id", "parDraw", "depth",
                 "PAR_atDepth", "K_N", "K_FAI")
  pop.f <- dir(sens.dir, glue("pop_...._0.1.rds"), full.names=T)
  mass.f <- dir(sens.dir, glue("mass_...._0.1.rds"), full.names=T)
  parSets <- map(dir(sens.dir, "000_.*csv", full.names=T), read_csv) %>%
    map(~.x %>% 
          mutate(param=str_remove(paste(param, stage, sep="."), ".NA")) %>%
          select(parDraw, param, val) %>% 
          pivot_wider(values_from="val", names_from="param"))
  params <- map(parSets, names) %>% unlist %>% unique
  params <- params[params!="parDraw"]
  
  cl <- makeCluster(nCores, outfile=glue("temp{sep}brt_out.txt"))
  obj.exclude <- c("data.ls", "grid.sf", "surv.df", "fecund.df", "lm.fit", "lm.mnsd")
  obj.include <- ls()
  
  cat("Exporting cluster.\n")
  clusterExport(cl, obj.include[-match(obj.exclude, obj.include, nomatch=0)])
  Sys.sleep(5)
  cat("Exported. Starting BRTs.\n")
  out.ls <- parLapply(cl, X=1:nrow(grid.i), fun=runBRTs,
                      pop.f=pop.f, mass.f=mass.f, parSets=parSets,
                      grid.i=grid.i, meta.cols=meta.cols, params=params,
                      brt.dir=glue("{sens.dir}{sep}BRTs{sep}"))
  stopCluster(cl)
  cat("Finished BRTs.\n")
}








# aggregate RIs -----------------------------------------------------------


dir(glue("{sens.dir}{sep}BRTs{sep}summaries{sep}"), "_ri_", full.names=T) %>%
  map_dfr(~read_csv(.x, show_col_types=F)) %>% 
  write_csv(glue("{sens.dir}{sep}RelInf_0.1.csv"))













