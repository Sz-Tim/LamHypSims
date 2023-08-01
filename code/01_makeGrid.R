# KELPER
# Make grids
# Tim Szewczyk

# This script constructs the grids to use in the simulations and extracts the
# mean environmental conditions within each cell.






# set up ------------------------------------------------------------------

# libraries and local functions
pkgs <- c("raster", "lubridate", "tidyverse", "sf", "glue")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00.*R", full.names=T), source)
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# switches
PAR_datasource <- c("MODIS", "POWER")[1]
gridRes <- 0.1 # arc-minutes

# directories
gis.dir <- glue("..{sep}..{sep}00_gis{sep}")

# maximum extent for covariates
UK_bbox <- st_bbox(c(xmin=-11, xmax=2, ymin=49, ymax=61), crs=st_crs(4326))

# datasets
# covarsFull.ls <- loadCovariates_full(gis.dir, UK_bbox, saveFile=glue("data{sep}covarFull_ls.rds"))
covars.ls <- loadCovariates(gis.dir, UK_bbox, loadFile=glue("data{sep}covar_ls.rds"))




# generate grid -----------------------------------------------------------

grid.mn <- UK_bbox %>%
  st_make_grid(cellsize=c(gridRes, gridRes)) %>%
  st_sf(id=1:length(.)) %>%
  extractCovarsToGrid(., covars.ls, PAR_datasource) %>%
  filter(!is.na(KD_sd) & !is.na(PAR_surface) & !is.na(fetch)) %>%
  mutate(id=row_number())




# save output -------------------------------------------------------------

# save as gpkg
st_write(grid.mn, glue("data{sep}grid_{gridRes}_{PAR_datasource}.gpkg"), append=F)

