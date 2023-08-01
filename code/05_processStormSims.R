# KELPER
# Process output
# Tim Szewczyk

# This script processes the output produced in 04_stormSimulations.R





# set up ------------------------------------------------------------------

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
theme_set(theme_bw())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# directories
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
out.dir <- glue("out{sep}storms{sep}")

# load files
grid.sf <- st_read(glue("data{sep}grid_0.1_MODIS.gpkg")) %>%
  select(id, geom, PAR_surface)
sim.info <- glue("...._0.1")
pop.f <- dir(out.dir, glue("pop_{sim.info}"), full.names=T) %>%
  grep("optPar", ., value=T) %>%
  grep("_g0", ., value=T)
pop.df <- map_dfr(pop.f, ~readRDS(.x) %>% 
                    select(id, grid.id, depth, year, month, date, stage, N) %>%
                    filter(month!=6 & stage != "recruits"))
mass.f <- dir(out.dir, glue("mass_{sim.info}"), full.names=T)
mass.df <- map_dfr(mass.f, ~readRDS(.x) %>% filter(month!=6))




# summarise output --------------------------------------------------------

y_vars <- c("FAI", "N", "biomass", "logN", "logBiomass", "kappa_N", "kappa_FAI")
x_vars <- c("K_N", "K_FAI", "SST", "KD", "PAR", "PAR_atDepth", "fetch", "fetchCat")
pop.sum <- pop.df %>% 
  select(-date, -year) %>%
  mutate(logN=log(N+1)) %>%
  group_by(id, sim, month, stage, depth, landscape, stochParams) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max, 
                            pOcc=~sum(.x>1)/n(), q90=~quantile(.x, probs=0.9))))
mass.sum <- mass.df %>%
  select(-date, -year) %>%
  mutate(logBiomass=log(biomass+1)) %>%
  group_by(id, sim, month, depth, landscape, stochParams) %>%
  summarise(across(any_of(c(x_vars, y_vars)), .names="{.col}_{.fn}", 
                   .fn=list(md=median, mn=mean, sd=sd, min=min, max=max, q90=~quantile(.x, probs=0.9))))

mass.df %>% group_by(id, date, depth, stochParams) %>%
  summarise(biomass_mn=mean(biomass), 
            biomass_md=median(biomass),
            biomass_sd=sd(biomass),
            biomass_025=quantile(biomass, probs=0.025),
            biomass_10=quantile(biomass, probs=0.1),
            biomass_25=quantile(biomass, probs=0.25),
            biomass_75=quantile(biomass, probs=0.75),
            biomass_90=quantile(biomass, probs=0.9),
            biomass_975=quantile(biomass, probs=0.975)) %>%
  write_csv(glue("summaries{sep}mass_simGrid_0.1.csv"))

pop.df %>% group_by(id, date, depth, stage) %>%
  summarise(N_mn=mean(N), 
            N_md=median(N),
            N_sd=sd(N),
            N_025=quantile(N, probs=0.025),
            N_10=quantile(N, probs=0.1),
            N_25=quantile(N, probs=0.25),
            N_75=quantile(N, probs=0.75),
            N_90=quantile(N, probs=0.9),
            N_975=quantile(N, probs=0.975)) %>%
  write_csv(glue("summaries{sep}pop_simGrid_0.1.csv"))




# save output -------------------------------------------------------------

saveRDS(pop.df, glue("summaries{sep}pop_df_0.1_optPar.rds"))
saveRDS(mass.df, glue("summaries{sep}mass_df_0.1.rds"))
saveRDS(pop.sum, glue("summaries{sep}pop_sum_0.1.rds"))
saveRDS(mass.sum, glue("summaries{sep}mass_sum_0.1.rds"))



