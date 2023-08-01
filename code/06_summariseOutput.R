# KELPER
# Summarise output
# Tim Szewczyk

# This script summarises output for use in 07_results.R



# setup -------------------------------------------------------------------

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms", "ggpubr")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00[a-z]_.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# directories
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
out.dir <- glue("out{sep}storms{sep}")
sens.dir <- glue("out{sep}sensitivity{sep}")

# grid and data
grid.sf <- st_read(glue("data{sep}grid_0.1_MODIS.gpkg")) %>%
  filter(id > 5)
grid.i <- grid.sf %>% st_drop_geometry() %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn) %>%
  select(id, SST, KD, PAR, fetch, fetchCat)
obs.ls <- readRDS(glue("data{sep}dfs_0.1.rds"))
data.ls <- compileDatasets(data.dir, supp.f)
sites.sf <- bind_rows(
  read_csv(glue("data{sep}raw{sep}collab{sep}collab_sites.csv"), skip=1) %>%
    select(lat, lon) %>% mutate(source="collab"),
  read_csv(glue("data{sep}raw{sep}digitized{sep}sitesDigitized.csv")) %>%
    select(lat, lon) %>% mutate(source="digitized")
) %>% 
  filter(complete.cases(.)) %>%
  st_as_sf(coords=c("lon", "lat"), crs=4326)






# output
mass.df <- readRDS(glue("summaries{sep}mass_df_0.1.rds")) %>% 
  filter(id > 5) %>%
  filter(!stochParams)


# lagged effects ----------------------------------------------------------

nLags <- 15

stormIndexLag <- data.ls$year_stormIndex %>% select(year) %>%
  bind_cols(imap_dfc(setNames(1:nLags, glue("lag_{1:nLags}_strm")), 
                     ~lag(data.ls$year_stormIndex$stormIndex, .x)))
gridSim.df <- readRDS(glue("data{sep}gridSim_{gridRes}.rds")) %>%
  select(id, grid.id, year, PAR, SST, KD)

b.df <- mass.df %>% 
  filter(month==7) %>% 
  select(id, grid.id, depth, year, date, biomass) %>%
  left_join(., gridSim.df) %>%
  mutate(year=year(date),
         PAR=PAR * exp(-depth * KD)) %>%
  group_by(id, grid.id, depth) %>%
  mutate(biomass=c(scale(biomass)),
         PAR=c(scale(PAR)),
         SST=c(scale(SST))) %>%
  multijetlag(PAR, SST, n=nLags) %>%
  ungroup %>%
  full_join(stormIndexLag, by="year") %>% 
  group_by(id, grid.id, depth) %>%
  summarise(across(starts_with("lag_"), ~coef(lm(biomass~.x))[2])) %>%
  pivot_longer(starts_with("lag_"), names_to="predictor", values_to="b") %>%
  mutate(lag=as.numeric(str_split_fixed(predictor, "_", 3)[,2]),
         covar=str_split_fixed(predictor, "_", 3)[,3])
b.df %>%
  group_by(lag, depth, covar) %>%
  summarise(b_mn=mean(b), b_md=median(b),
            b_q10=quantile(b, probs=0.1),
            b_q90=quantile(b, probs=0.9)) %>%
  saveRDS(glue("out/processed/b_lagEffects_biomass.rds"))
b.df %>%
  group_by(lag, depth, covar, id) %>%
  summarise(b_mn=mean(b), b_md=median(b),
            b_q10=quantile(b, probs=0.1),
            b_q90=quantile(b, probs=0.9)) %>%
  saveRDS(glue("out/processed/b_lagEffects_id_biomass.rds"))

for(i in c("recruits", "subcanopy", "canopy")) {
  b.df <- readRDS(glue("summaries{sep}pop_df_0.1_{i}.rds")) %>%
    filter(id > 5, month==7) %>%
    select(id, grid.id, depth, year, date, N) %>%
    left_join(., gridSim.df) %>%
    mutate(year=year(date),
           PAR=PAR * exp(-depth * KD)) %>%
    group_by(id, grid.id, depth) %>%
    mutate(N=c(scale(log(N+1))),
           PAR=c(scale(PAR)),
           SST=c(scale(SST))) %>%
    multijetlag(PAR, SST, n=nLags) %>%
    ungroup %>%
    full_join(stormIndexLag, by="year") %>%
    group_by(id, grid.id, depth) %>%
    summarise(across(starts_with("lag_"), ~coef(lm(N~.x))[2])) %>%
    pivot_longer(starts_with("lag_"), names_to="predictor", values_to="b") %>%
    mutate(lag=as.numeric(str_split_fixed(predictor, "_", 3)[,2]),
           covar=str_split_fixed(predictor, "_", 3)[,3])
  b.df %>%
    group_by(lag, depth, covar) %>%
    summarise(b_mn=mean(b), b_md=median(b),
              b_q10=quantile(b, probs=0.1),
              b_q90=quantile(b, probs=0.9)) %>%
    saveRDS(glue("out/processed/b_lagEffects_{i}.rds"))
  b.df %>%
    group_by(lag, depth, covar, id) %>%
    summarise(b_mn=mean(b), b_md=median(b),
              b_q10=quantile(b, probs=0.1),
              b_q90=quantile(b, probs=0.9)) %>%
    saveRDS(glue("out/processed/b_lagEffects_id_{i}.rds"))
}




# detritus ----------------------------------------------------------------

loss.df <- mass.df %>%
  filter(landscape=="dynamic" & !stochParams) %>%
  select(-landscape, -stochParams) %>%
  group_by(id, grid.id, depth) %>%
  arrange(date) %>%
  mutate(delta_mass=biomass-lag(biomass),
         delta_pct=delta_mass/lag(biomass)) %>%
  ungroup
saveRDS(loss.df, "out/processed/loss_df.rds")
loss.df %>% 
  filter(delta_mass < 0) %>%
  filter(month==1 & !is.na(delta_mass)) %>%
  mutate(delta_mass=-delta_mass) %>%
  group_by(id, depth) %>%
  arrange(delta_mass) %>%
  mutate(decile=paste0(floor(row_number()/(n()+1)*10)*10, "th")) %>%
  group_by(id, depth, decile) %>%
  summarise(delta_mass=sum(delta_mass)) %>%
  group_by(id, depth) %>%
  mutate(delta_prop=delta_mass/sum(delta_mass)) %>%
  saveRDS("out/processed/loss_deciles.rds")


# parameter ranges --------------------------------------------------------

par.rng <- readRDS(glue("{sens.dir}000_parameter_ranges.rds")) %>%
  mutate(exposure=case_when(exposure=="low"~1, exposure=="high"~2)) 
surv.rng <- par.rng %>% filter(param=="surv") %>%
  group_by(exposure) %>% group_split()
loss.rng <- par.rng %>% filter(param=="loss")

data.ls$year_stormIndex %>%
  select(year, stormIndex) %>%
  mutate(surv_canopy_1=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[1]], stage=="canopy")$shp1,
                             filter(surv.rng[[1]], stage=="canopy")$shp2),
         surv_canopy_2=qbeta(pnorm(-stormIndex, 0, 1), 
                             filter(surv.rng[[2]], stage=="canopy")$shp1,
                             filter(surv.rng[[2]], stage=="canopy")$shp2),
         surv_subcanopy_1=qbeta(pnorm(-stormIndex, 0, 1), 
                                filter(surv.rng[[1]], stage=="subcanopy")$shp1,
                                filter(surv.rng[[1]], stage=="subcanopy")$shp2),
         surv_subcanopy_2=qbeta(pnorm(-stormIndex, 0, 1), 
                                filter(surv.rng[[2]], stage=="subcanopy")$shp1,
                                filter(surv.rng[[2]], stage=="subcanopy")$shp2),
         surv_recruit_1=qbeta(pnorm(-stormIndex, 0, 1), 
                              filter(surv.rng[[1]], stage=="recruits")$shp1,
                              filter(surv.rng[[1]], stage=="recruits")$shp2),
         surv_recruit_2=qbeta(pnorm(-stormIndex, 0, 1), 
                              filter(surv.rng[[2]], stage=="recruits")$shp1,
                              filter(surv.rng[[2]], stage=="recruits")$shp2)) %>%
  pivot_longer(starts_with("surv_"), names_to="param", values_to="value") %>%
  mutate(stage=str_split_fixed(param, "_", 3)[,2] %>%
           factor(levels=c("recruit", "subcanopy", "canopy")),
         exposure=str_split_fixed(param, "_", 3)[,3] %>% 
           factor(labels=c("Low", "High"))) %>%
  mutate(param=paste0("italic(s[", stage, "])"),
         value=sqrt(value)) %>%
  bind_rows(data.ls$year_stormIndex %>%
              select(year, stormIndex) %>%
              mutate(param="italic(epsilon)",
                     exposure="Low",
                     value=qbeta(pnorm(stormIndex, 0, 1),
                                 loss.rng$shp1, loss.rng$shp2))) %>%
  bind_rows(data.ls$year_stormIndex %>%
              select(year, stormIndex) %>%
              mutate(param="italic(epsilon)",
                     exposure="High",
                     value=qbeta(pnorm(stormIndex, 0, 1),
                                 loss.rng$shp1, loss.rng$shp2))) %>%
  mutate(param=factor(param, levels=rev(unique(param))),
         exposure=factor(exposure, levels=c("Low", "High"))) %>%
  saveRDS("out/processed/storm_effect_df.rds")

st_read(glue("data{sep}grid_0.1_MODIS.gpkg")) %>%
  filter(id > 5) %>%
  rename(PAR=PAR_mn, SST=sstDay_mn, KD=KD_mn) %>%
  mutate(PAR2=PAR*exp(-KD*2), PAR5=PAR*exp(-KD*5), PAR10=PAR*exp(-KD*10),
         PAR15=PAR*exp(-KD*15), PAR20=PAR*exp(-KD*20)) %>%
  select(id, SST, fetch, fetchCat, PAR, PAR2, PAR5, PAR10, PAR15, PAR20) %>%
  pivot_longer(2:10, names_to="predictor", values_to="value") %>%
  mutate(predictor=factor(predictor, 
                          levels=c("SST", "fetch", "fetchCat", "PAR", "PAR2",
                                   "PAR5", "PAR10", "PAR15", "PAR20"),
                          labels=c("SST (growing season)", "ln(fetch)", 
                                   "Exposure category",
                                   "PAR (surface)", "PAR (2m)", "PAR (5m)",
                                   "PAR (10m)", "PAR (15m)", "PAR (20m)"))) %>%
  st_write("out/processed/grid_long.gpkg")

