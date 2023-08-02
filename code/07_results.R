# KELPER
# Summarise results and make figures
# Tim Szewczyk

# This script makes the figures for publication and summarises results




# setup -------------------------------------------------------------------

# libraries and local functions
pkgs <- c("raster", "lubridate", "glue", "tidyverse", "sf", "brms", "ggpubr")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
walk(dir("code", "^00[a-z]_.*R", full.names=T), source)
theme_set(theme_classic())
sep <- ifelse(.Platform$OS.type=="unix", "/", "\\")

# directories
data.dir <- glue("data{sep}raw{sep}digitized{sep}")
supp.f <- glue("data{sep}raw{sep}collab{sep}collab_all.xlsx")
out.dir <- glue("out{sep}storms{sep}")
sens.dir <- glue("out{sep}sensitivity{sep}")
fig.dir <- glue("../drafts/LamHyp_ms/EcologicalModelling/figs/")

# grid and data
grid.sf <- st_read(glue("data{sep}grid_0.1_MODIS.gpkg")) %>%
  filter(id > 5)
grid.i <- grid.sf %>% st_drop_geometry() %>%
  rename(SST=sstDay_mn, PAR=PAR_surface, KD=KD_mn) %>%
  select(id, SST, KD, PAR, fetch, fetchCat)
gridSim.sum <- readRDS(glue("data{sep}gridSim_0.1.rds")) %>%
  select(id, grid.id, year, PAR, SST, KD) %>%
  group_by(id) %>%
  summarise(SST_mn=mean(SST),
            PAR_0=mean(PAR),
            PAR_2=mean(PAR*exp(-KD*2)),
            PAR_5=mean(PAR*exp(-KD*5)),
            PAR_10=mean(PAR*exp(-KD*10)),
            PAR_15=mean(PAR*exp(-KD*15)),
            PAR_20=mean(PAR*exp(-KD*20))) %>%
  ungroup %>%
  pivot_longer(starts_with("PAR_"), names_to="depth", values_to="PAR_atDepth") %>%
  mutate(depth=as.numeric(str_sub(depth, 5, -1))) %>%
  full_join(., grid.i %>% select(id, fetch, fetchCat))
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

# plot themes, etc
pub_theme <- theme(axis.title=element_text(size=9),
                   axis.text=element_text(size=8),
                   legend.title=element_text(size=9),
                   legend.text=element_text(size=8),
                   strip.text=element_text(size=9),
                   strip.background=element_rect(linewidth=0.4),
                   axis.line=element_line(linewidth=0.3),
                   axis.ticks=element_line(linewidth=0.25))
map_base.gg <- ggplot() +
  theme(legend.position="bottom", 
        panel.grid=element_blank(),
        axis.title=element_blank()) +
  scale_x_continuous(breaks=c(-8,0)) +
  scale_y_continuous(breaks=c(52, 58)) 

# output
mass.sum <- readRDS(glue("summaries{sep}mass_sum_0.1.rds")) %>% 
  filter(id > 5) %>%
  left_join(., gridSim.sum)
mass.sum.sf_canopy <- full_join(grid.sf, mass.sum, by="id") %>% 
  filter(month==7 & stochParams==F) %>%
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         month=month.name[month])
pop.sum_canopy <- readRDS(glue("summaries{sep}pop_sum_0.1.rds")) %>% 
  filter(stage=="canopy") %>%
  filter(id > 5) %>% 
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         month=month.name[month])
b.df_biomass <- readRDS("out/processed/b_lagEffects_biomass.rds")
ri.df <- read_csv(glue("{sens.dir}{sep}RelInf_0.1.csv")) %>%
  group_by(response, id, month, depth) %>%
  filter(smp==max(smp), td==max(td))
stages <- c("can.", "recr.", "subc.")
stages_ord <- c("can.", "subc.", "recr.")
par.rng <- readRDS(glue("{sens.dir}000_parameter_ranges.rds")) %>%
  mutate(exposure=case_when(exposure=="low"~1, exposure=="high"~2)) 
sens.param <- tibble(messy=sort(unique(ri.df$var)),
                     full=c("Density shape", 
                            paste("Frond growth:", stages),
                            paste("Stipe growth:", stages),
                            "Erosion", "Settlement", 
                            paste("Survival:", stages)),
                     eq=c("theta", 
                          paste0("omega[", stages, "]"),
                          paste0("gamma[", stages, "]"),
                          "epsilon", "z",
                          paste0("s[", stages, "]")),
                     grp=c("NDD",
                           rep("Growth: Frond", 3),
                           rep("Growth: Stipe", 3),
                           "Erosion", "Repro.",
                           rep("Survival", 3))) %>%
  mutate(eq=paste0("italic(", eq, ")")) %>%
  mutate(eq=factor(eq, 
                   levels=c("italic(epsilon)",
                            "italic(theta)", 
                            "italic(z)",
                            paste0("italic(omega[", rev(stages_ord), "])"),
                            paste0("italic(gamma[", rev(stages_ord), "])"),
                            paste0("italic(s[", rev(stages_ord), "])")
                   )))
loss.deciles <- readRDS("out/processed/loss_deciles.rds")
storm_effect.df <- readRDS("out/processed/storm_effect_df.rds")
grid.long <- st_read("out/processed/grid_long.gpkg")



# Fig 2 -------------------------------------------------------------------

fig2a <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy, 
          colour=NA, aes(fill=logBiomass_mn)) +
  scale_fill_gradient("Canopy\nbiomass\nmean\n(ln kg/m2)", 
                      low="grey90", high="green4", limits=c(0,NA)) +
  facet_grid(.~depth) + 
  pub_theme +
  theme(legend.position="right",
        legend.key.width=unit(0.2, "cm"))
fig2b <- map_base.gg +
  geom_sf(data=mass.sum.sf_canopy, 
          colour=NA, aes(fill=logBiomass_sd)) +
  scale_fill_viridis_c("Canopy\nbiomass\nstd. dev.\n(ln kg/m2)", 
                       limits=c(0, NA)) +
  facet_grid(.~depth) + 
  pub_theme +
  theme(legend.position="right",
        legend.key.width=unit(0.2, "cm"))
ggsave(glue("{fig.dir}Fig_2.png"), 
       ggarrange(fig2a, fig2b, ncol=1, labels="auto"),
       width=190, height=120, dpi=300, units="mm")
ggsave(glue("{fig.dir}Fig_2.pdf"), 
       ggarrange(fig2a, fig2b, ncol=1, labels="auto"),
       width=190, height=120, dpi=300, units="mm")


# Fig 3 -------------------------------------------------------------------

fig3a <- loss.deciles %>%
  filter(depth %in% c(2, 5, 10, 15)) %>%
  mutate(depth_F=factor(depth, levels=c(2, 5, 10, 15)),
         grp=factor(paste(depth_F, decile))) %>%
  ggplot(aes(decile, delta_mass/1e3, fill=depth, group=grp)) + 
  geom_boxplot(size=0.2, outlier.size=0.1, colour="grey30") +
  scale_fill_viridis_c("Depth (m)", direction=-1) +
  guides(fill=guide_colorbar(reverse=T)) +
  labs(x="Decile among years", y=expression(paste("Winter detritus (kg/", m^2, ")"))) + 
  pub_theme +
  theme(panel.grid.major.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position=c(0.15, 0.75),
        legend.key.width=unit(0.2, "cm"))

fig3b <- loss.deciles %>%
  filter(depth %in% c(2, 5, 10, 15)) %>%
  mutate(depth_F=factor(depth, levels=c(2, 5, 10, 15)),
         grp=factor(paste(depth_F, decile))) %>%
  ggplot(aes(decile, delta_prop, fill=depth, group=grp)) + 
  geom_boxplot(size=0.2, outlier.size=0.1, colour="grey30") +
  scale_fill_viridis_c("Depth (m)", direction=-1) +
  guides(fill=guide_colorbar(reverse=T)) +
  labs(x="Decile among years", y="Proportion of total detritus within cell") + 
  pub_theme +
  theme(panel.grid.major.x=element_blank(),
        legend.position="none",
        legend.key.width=unit(0.2, "cm"))

ggsave(glue(fig.dir, "Fig_3.png"), 
       ggarrange(fig3a, fig3b, ncol=1, labels="auto", heights=c(0.95, 1)), 
       width=90, height=150, dpi=300, units="mm")
ggsave(glue(fig.dir, "Fig_3.pdf"), 
       ggarrange(fig3a, fig3b, ncol=1, labels="auto", heights=c(0.95, 1)), 
       width=90, height=150, dpi=300, units="mm")



# Fig 4 -------------------------------------------------------------------

fig4 <- b.df_biomass %>% filter(depth %in% c(2,5,10,15)) %>%
  mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
         covar=factor(covar, levels=c("strm", "PAR", "SST"),
                      labels=c("Storms", "PAR", "SST"))) %>%
  ggplot(aes(lag, b_md, ymin=b_q10, ymax=b_q90, colour=covar, fill=covar)) +
  geom_hline(yintercept=0, colour="grey80") +
  geom_ribbon(alpha=0.3, colour=NA) +
  geom_line() +
  scale_colour_manual("Variable", values=c("grey40", "goldenrod", "red4")) +
  scale_fill_manual("Variable", values=c("grey40", "goldenrod", "red4")) +
  facet_grid(depth~.) + 
  scale_x_continuous(breaks=c(5, 10, 15)) +
  scale_y_continuous(breaks=c(-0.5, 0, 0.5)) +
  pub_theme +
  theme(legend.position=c(0.86, 0.05),
        legend.key.height=unit(0.25, "cm"),
        legend.key.width=unit(0.25, "cm"),
        legend.background=element_blank()) +
  labs(x="Lag (years)", 
       y=glue("Effect on biomass\n(standardized slope: median + middle 80%)"))
ggsave(glue("{fig.dir}Fig_4.png"), fig4,
       width=90, height=180, dpi=300, units="mm")
ggsave(glue("{fig.dir}Fig_4.pdf"), fig4, 
       width=90, height=180, dpi=300, units="mm")





# Fig 5 -------------------------------------------------------------------

fig5 <- ri.df %>%
  filter(var != "growStipe.canopy") %>%
  filter(month=="july") %>%
  group_by(var, month, depth, response) %>%
  summarise(md=median(rel.inf*100),
            q10=quantile(rel.inf*100, probs=0.1),
            q25=quantile(rel.inf*100, probs=0.25),
            q75=quantile(rel.inf*100, probs=0.75),
            q90=quantile(rel.inf*100, probs=0.9)) %>%
  ungroup %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=as.numeric(str_remove(depth, "m")),
         param=sens.param$eq[match(var, sens.param$messy)],
         grp=sens.param$grp[match(var, sens.param$messy)]) %>%
  mutate(grp=factor(grp, levels=c("Survival", "Repro.", "Erosion", "NDD",
                                  "Growth: Frond", "Growth: Stipe"))) %>%
  ggplot(aes(param, md, colour=depth, group=as.factor(depth))) +
  geom_hline(yintercept=0, colour="grey90") +
  geom_errorbar(aes(ymin=q10, ymax=q90), position=position_dodge(width=-0.75), width=0.5) +
  geom_linerange(aes(ymin=q25, ymax=q75), position=position_dodge(width=-0.75), linewidth=1) +
  geom_point(position=position_dodge(width=-0.75)) +
  facet_grid(grp~response, scales="free_y", space="free_y", switch="y") + 
  coord_flip() +
  scale_colour_viridis_c("Depth (m)", direction=-1, end=0.95) +
  guides(colour=guide_colorbar(reverse=T)) +
  scale_x_discrete(labels=scales::parse_format()) +
  scale_y_continuous(breaks=c(0, 20, 40)) + 
  pub_theme +
  theme(legend.position=c(0.86, 0.15),
        legend.key.width=unit(0.2, "cm"),
        legend.key.height=unit(0.4, "cm"),
        legend.background=element_blank()) +
  labs(y="Relative influence (%)", x="")
ggsave(glue("{fig.dir}Fig_5.png"), fig5, 
       width=90, height=160, dpi=300, units="mm")
ggsave(glue("{fig.dir}Fig_5.pdf"), fig5,
       width=90, height=160, dpi=300, units="mm")


# Fig 6 -------------------------------------------------------------------

fig6_df <- ri.df %>% filter(month=='july') %>% 
  filter(depth %in% paste0(c(2, 10), "m")) %>%
  filter(response=="biomass_mn") %>%
  group_by(var) %>%
  mutate(max_ri=max(rel.inf)) %>%
  ungroup %>%
  filter(max_ri > 0.1) %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20),"m")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  left_join(grid.sf, .) 
fig6 <- map_base.gg + 
  geom_sf(data=fig6_df, aes(fill=log10(rel.inf*100)), colour=NA) +
  scale_fill_viridis_c(expression(log[10](Rel.~Inf.)), option="cividis", 
                       limits=c(NA,2), breaks=c(0,1,2)) +
  pub_theme +
  theme(legend.key.height=unit(0.2, "cm")) +
  facet_grid(depth~param, 
             labeller=labeller(response=label_value, depth=label_value, param=label_parsed))
ggsave(glue("{fig.dir}/Fig_6.png"), fig6, 
       width=140, height=87, dpi=300, units="mm")
ggsave(glue("{fig.dir}/Fig_6.pdf"), fig6, 
       width=140, height=87, dpi=300, units="mm")


ri.df %>% filter(month=='july') %>% 
  group_by(var) %>%
  mutate(max_ri=max(rel.inf)) %>%
  ungroup %>%
  arrange(response, id, depth, desc(rel.inf)) %>%
  group_by(response, id, depth) %>%
  slice_head(n=1) %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20),"m")),
         param=sens.param$eq[match(var, sens.param$messy)]) %>%
  left_join(grid.sf, .) %>%
  ggplot(aes(fill=var)) + 
  geom_sf(colour=NA) + 
  scale_fill_manual(values=c("#1b7837", "#a6dba0", "#762a83", "#c2a5cf", "#9970ab")) +
  pub_theme + 
  facet_grid(response~depth)




# Supp A ------------------------------------------------------------------


# * Fig A.1 ---------------------------------------------------------------

map_base.gg +
  geom_sf(data=grid.sf, colour=NA, fill="#80cdc1") + 
  geom_sf(data=sites.sf %>% st_crop(grid.sf), colour="#8c510a") +
  pub_theme +
  ggtitle("UK data source locations")
ggsave(glue("{fig.dir}Fig_A1.png"), width=3, height=4.5, dpi=300)



# * Fig A.2 ---------------------------------------------------------------

lm.fit <- readRDS(glue("data{sep}fits_0.1.rds"))
out.eff <- map(lm.fit, ~conditional_effects(.x, surface=T))

ggarrange(plotlist=plot(out.eff[[1]], ask=F, points=T, plot=F, stype="raster",
                        point_args=list(shape=1, size=0.75, alpha=0.8), 
                        theme=pub_theme), 
          nrow=1, widths=c(1,1,1.2))
ggsave(glue("{fig.dir}Fig_A2a.png"), width=12, height=4)

ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[2]], 
                                    effects=c("logLenStipe", "lPAR_atDepth", "SST")), 
                ask=F, points=T, plot=F, stype="raster",
                point_args=list(shape=1, size=0.75, alpha=0.8), 
                theme=pub_theme), 
  nrow=1, widths=c(1,1,1))
ggsave(glue("{fig.dir}Fig_A2b-1.png"), width=12, height=4)
ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[2]], surface=T,
                                    effects=c("logLenStipe:lPAR_atDepth", "lPAR_atDepth:SST")), 
                ask=F, points=T, plot=F, stype="raster",
                point_args=list(shape=1, size=0.75, alpha=0.8), 
                theme=pub_theme), 
  nrow=1, widths=c(1,1))
ggsave(glue("{fig.dir}Fig_A2b-2.png"), width=11, height=4)

ggarrange(plotlist=plot(out.eff[[4]], ask=F, points=T, plot=F, stype="raster",
                        point_args=list(shape=1, size=0.75, alpha=0.8), 
                        theme=pub_theme), 
          nrow=1, widths=c(1,1,1.2))
ggsave(glue("{fig.dir}Fig_A2c.png"), width=12, height=4)

ggarrange(plotlist=plot(out.eff[[3]], ask=F, points=T, plot=F, stype="raster",
                        point_args=list(shape=1, size=0.75, alpha=0.8), 
                        theme=pub_theme), 
          nrow=1, widths=c(1,1,1.2))
ggsave(glue("{fig.dir}Fig_A2d.png"), width=12, height=4)

ggarrange(plotlist=plot(out.eff[[5]], ask=F, points=T, plot=F, 
                        point_args=list(shape=1, size=0.75, alpha=0.8), 
                        theme=pub_theme))
ggsave(glue("{fig.dir}Fig_A2e.png"), width=4, height=4)

ggarrange(plotlist=plot(out.eff[[6]], ask=F, points=T, plot=F, 
                        point_args=list(shape=1, size=0.75, alpha=0.8), 
                        theme=pub_theme), 
          nrow=1)
ggsave(glue("{fig.dir}Fig_A2f.png"), width=8, height=4)

ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[7]], 
                                    effects=c("PAR_atDepth", "SST", "fetch")), ask=F,
                points=T, plot=F, 
                point_args=list(shape=1, size=0.75, alpha=0.8), 
                theme=pub_theme),
  nrow=1)
ggsave(glue("{fig.dir}Fig_A2g-1.png"), width=12, height=4)
ggarrange(
  plotlist=plot(conditional_effects(lm.fit[[7]], surface=T,
                                    effects=c("PAR_atDepth:SST", "PAR_atDepth:fetch")), ask=F, 
                points=T, plot=F, stype="raster",
                point_args=list(shape=1, size=0.75, alpha=0.8), 
                theme=pub_theme),
  nrow=1)
ggsave(glue("{fig.dir}Fig_A2g-2.png"), width=10, height=4)





# * Fig A.3 ---------------------------------------------------------------

ggplot(data.ls$year_stormIndex, aes(year, stormIndex)) + 
  stat_smooth(method="loess", colour="grey", se=F, span=0.5, size=0.5) +
  geom_point() + geom_line(size=0.25) +
  pub_theme +
  labs(x="Year", y="Storm index")
ggsave(glue("{fig.dir}Fig_A3.png"), height=4, width=6, dpi=300)



# * Fig A.4 ---------------------------------------------------------------

ggplot(storm_effect.df, aes(year, value, colour=param)) + 
  geom_line() + 
  pub_theme +
  scale_x_continuous(breaks=c(1950, 1975, 2000)) +
  scale_y_continuous("Rate (non-growing season)", breaks=c(0, 0.5, 1)) +
  scale_colour_manual("Parameter", 
                      values=c("#bf812d", "#80cdc1", "#35978f", "#003c30"), 
                      labels=scales::parse_format()) +
  theme(axis.title.x=element_blank(),
        legend.text.align=0) +
  facet_grid(.~exposure) 
ggsave(glue("{fig.dir}Fig_A4.png"), width=7, height=2.5, dpi=300)




# * Fig A.5 ---------------------------------------------------------------

fetch_cols <- c("grey70", "black")
stage_cols <- c(recruits="#80cdc1", subcanopy="#35978f", canopy="#003c30")

# s_recruits
s_r <- list(xlim=c(0,1), main="Recruit survival rate",
            xlab=expression(italic(s[recruits])), ylab="density")
s_r$x.full <- seq(s_r$xlim[1], 0.9, length.out=1e3)
s_r$c.full <- map(1:2, 
                  ~with(filter(par.rng, param=="surv" & stage=="recruits" & exposure==.x),
                        dbeta(s_r$x.full, shp1, shp2)))
s_r$x.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="recruits" & exposure==.x),
                           c(valMin, valMax)))
s_r$c.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="recruits" & exposure==.x),
                           dbeta(s_r$x.95[[.x]], shp1, shp2)))
s_r$ylim <- c(0, max(unlist(s_r$c.full)[is.finite(unlist(s_r$c.full))]))

# s_subcanopy
s_s <- list(xlim=c(0,1),main="Subcanopy survival rate",
            xlab=expression(italic(s[subcanopy])), ylab="density")
s_s$x.full <- seq(s_s$xlim[1], s_s$xlim[2], length.out=1e3)
s_s$c.full <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="subcanopy" & exposure==.x),
                             dbeta(s_s$x.full, shp1, shp2)))
s_s$x.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="subcanopy" & exposure==.x),
                           c(valMin, valMax)))
s_s$c.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="subcanopy" & exposure==.x),
                           dbeta(s_s$x.95[[.x]], shp1, shp2)))
s_s$ylim <- c(0, max(unlist(s_s$c.full)[is.finite(unlist(s_s$c.full))]))

# s_canopy
s_c <- list(xlim=c(0,1), main="Canopy survival rate",
            xlab=expression(italic(s[canopy])), ylab="density")
s_c$x.full <- seq(s_c$xlim[1], s_c$xlim[2], length.out=5e2)
s_c$c.full <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="canopy" & exposure==.x),
                             dbeta(s_c$x.full, shp1, shp2)))
s_c$x.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="canopy" & exposure==.x),
                           c(valMin, valMax)))
s_c$c.95 <- map(1:2, ~with(filter(par.rng, param=="surv" & stage=="canopy" & exposure==.x),
                           dbeta(s_c$x.95[[.x]], shp1, shp2)))
s_c$ylim <- c(0, max(unlist(s_c$c.full)[is.finite(unlist(s_c$c.full))]))

# z
z <- list(xlim=c(0,5e3), main="Settlement rate",
          xlab=expression(italic(z)~(N/m^2)), ylab="density")
z$x.full <- seq(z$xlim[1], z$xlim[2], length.out=5e2)
z$c.full <- map(1:2, ~with(filter(par.rng, param=="settlement" & exposure==.x),
                           dnorm(z$x.full, valMean, valSD)))
z$x.95 <- map(1:2, ~with(filter(par.rng, param=="settlement" & exposure==.x),
                         c(valMin, valMax)))
z$c.95 <- map(1:2, ~with(filter(par.rng, param=="settlement" & exposure==.x),
                         dnorm(z$x.95[[.x]], valMean, valSD)))
z$ylim <- c(0, max(unlist(z$c.full)))

# growStipe
gamma <- list(xlim=c(0,300), main="Stipe growth rate",
              xlab=expression(italic(gamma[x])~(mm/year)), ylab="density")
gamma$x.full <- seq(gamma$xlim[1], gamma$xlim[2], length.out=5e2)
gamma$c.full <- map(1:3, ~with(filter(par.rng, param=="growStipe" & stage==names(stage_cols)[.x]),
                               dnorm(gamma$x.full, valMean, valSD)))
gamma$x.95 <- map(1:3, ~with(filter(par.rng, param=="growStipe" & stage==names(stage_cols)[.x]),
                             c(valMin, valMax)))
gamma$c.95 <- map(1:3, ~with(filter(par.rng, param=="growStipe" & stage==names(stage_cols)[.x]),
                             dnorm(gamma$x.95[[.x]], valMean, valSD)))
gamma$ylim <- c(0, max(unlist(gamma$c.full)))

# growFrond
omega <- list(xlim=c(0,0.6), main="Frond growth rate",
              xlab=expression(italic(omega[x])~(m^2/year)), ylab="density")
omega$x.full <- seq(omega$xlim[1], omega$xlim[2], length.out=5e2)
omega$c.full <- map(1:3, ~with(filter(par.rng, param=="growFrond" & stage==names(stage_cols)[.x]),
                               dnorm(omega$x.full, valMean, valSD)))
omega$x.95 <- map(1:3, ~with(filter(par.rng, param=="growFrond" & stage==names(stage_cols)[.x]),
                             c(valMin, valMax)))
omega$c.95 <- map(1:3, ~with(filter(par.rng, param=="growFrond" & stage==names(stage_cols)[.x]),
                             dnorm(omega$x.95[[.x]], valMean, valSD)))
omega$ylim <- c(0, max(unlist(omega$c.full)))

# erosion
epsilon <- list(xlim=c(0,1), main="Frond erosion rate",
                xlab=expression(italic(epsilon)~(prop.~frond/non-growing~season)), ylab="density")
epsilon$x.full <- seq(epsilon$xlim[1], epsilon$xlim[2], length.out=5e2)
epsilon$c.full <- with(filter(par.rng, param=="loss"), 
                       dbeta(epsilon$x.full, shp1, shp2))
epsilon$x.95 <- with(filter(par.rng, param=="loss"), c(valMin, valMax))
epsilon$c.95 <- with(filter(par.rng, param=="loss"),
                     dbeta(epsilon$x.95, shp1, shp2))
epsilon$ylim <- c(0, max(unlist(epsilon$c.full)))

# density effect shape
theta <- list(xlim=c(0,2), main="Density effect shape",
              xlab=expression(italic(theta)), ylab="density")
theta$x.full <- seq(theta$xlim[1], theta$xlim[2], length.out=5e2)
theta$c.full <- dnorm(theta$x.full, 1, 0.25)
theta$x.95 <- with(filter(par.rng, param=="densityEffShape"), c(valMin, 1.5))
theta$c.95 <- with(filter(par.rng, param=="densityEffShape"), dnorm(theta$x.95, 1, 0.25))
theta$ylim <- c(0, max(unlist(theta$c.full)))




png(glue("{fig.dir}Fig_A5.png"), width=9, height=4.5, res=300, units="in")
{
  par(mfrow=c(2,4), mar=c(5,4,2,1)+0.1)
  
  # s_recruits
  plot(NA, NA, xlim=s_r$xlim, ylim=s_r$ylim, main=s_r$main, xlab=s_r$xlab, ylab=s_r$ylab)
  for(i in 1:2) {
    segments(x0=s_r$x.95[[i]], y0=0, y1=s_r$c.95[[i]], col=fetch_cols[i])
    lines(s_r$x.full, s_r$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low",  "High"), bty="n", title="Wave fetch")
  
  # s_subcanopy
  plot(NA, NA, xlim=s_s$xlim, ylim=s_s$ylim, main=s_s$main, xlab=s_s$xlab, ylab=s_s$ylab)
  for(i in 1:2) {
    segments(x0=s_s$x.95[[i]], y0=0, y1=s_s$c.95[[i]], col=fetch_cols[i])
    lines(s_s$x.full, s_s$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low", "High"), bty="n", title="Wave fetch")
  
  # s_canopy
  plot(NA, NA, xlim=s_c$xlim, ylim=s_c$ylim, main=s_c$main, xlab=s_c$xlab, ylab=s_c$ylab)
  for(i in 1:2) {
    segments(x0=s_c$x.95[[i]], y0=0, y1=s_c$c.95[[i]], col=fetch_cols[i])
    lines(s_c$x.full, s_c$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low", "High"), bty="n", title="Wave fetch")
  
  # z
  plot(NA, NA, xlim=z$xlim, ylim=z$ylim, main=z$main, xlab=z$xlab, ylab=z$ylab)
  for(i in 1:2) {
    segments(x0=z$x.95[[i]], y0=0, y1=z$c.95[[i]], col=fetch_cols[i])
    lines(z$x.full, z$c.full[[i]], col=fetch_cols[i])    
  }
  legend("topright", col=fetch_cols, lty=1, c("Low", "High"), bty="n", title="Wave fetch")
  
  # gamma
  plot(NA, NA, xlim=gamma$xlim, ylim=gamma$ylim, main=gamma$main, xlab=gamma$xlab, ylab=gamma$ylab)
  for(i in 1:3) {
    segments(x0=gamma$x.95[[i]], y0=0, y1=gamma$c.95[[i]], col=stage_cols[i])
    lines(gamma$x.full, gamma$c.full[[i]], col=stage_cols[i])    
  }
  legend("topright", col=stage_cols, lty=1, names(stage_cols), bty="n", title="Stage")
  
  # omega
  plot(NA, NA, xlim=omega$xlim, ylim=omega$ylim, main=omega$main, xlab=omega$xlab, ylab=omega$ylab)
  for(i in 1:3) {
    segments(x0=omega$x.95[[i]], y0=0, y1=omega$c.95[[i]], col=stage_cols[i])
    lines(omega$x.full, omega$c.full[[i]], col=stage_cols[i])    
  }
  legend("topright", col=stage_cols, lty=1, names(stage_cols), bty="n", title="Stage")
  
  # epsilon
  plot(NA, NA, xlim=epsilon$xlim, ylim=epsilon$ylim, main=epsilon$main, xlab=epsilon$xlab, ylab=epsilon$ylab)
  segments(x0=epsilon$x.95, y0=0, y1=epsilon$c.95)
  lines(epsilon$x.full, epsilon$c.full)    
  
  # theta
  plot(NA, NA, xlim=theta$xlim, ylim=theta$ylim, main=theta$main, xlab=theta$xlab, ylab=theta$ylab)
  segments(x0=theta$x.95, y0=0, y1=theta$c.95)
  lines(theta$x.full, theta$c.full)    
}
dev.off()



# * Fig A.6 ---------------------------------------------------------------

sst.p <- map_base.gg +
  geom_sf(data=filter(grid.long, predictor=="SST (growing season)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c("°C", option="C") + 
  ggtitle("SST (growing season)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
fetch.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="ln(fetch)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c("", option="B", end=0.9) + 
  ggtitle(expression(log[10](wave~fetch))) + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
fetchCat.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="Exposure category"), colour=NA, aes(fill=value)) +
  pub_theme + 
  scale_fill_viridis_c("", option="B", end=0.75, guide=guide_legend(), breaks=1:2, labels=c("low", "high")) + 
  ggtitle("Exposure category") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR0.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (surface)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (surface)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR2.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (2m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (2m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR5.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (5m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (5m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR10.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (10m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (10m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR15.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (15m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (15m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR20.p <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (20m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day)), limits=c(0, 29)) +
  ggtitle("PAR (20m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))

PAR0.p_ <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (surface)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  ggtitle("PAR (surface)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR2.p_ <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (2m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  ggtitle("PAR (2m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR5.p_ <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (5m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  ggtitle("PAR (5m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR10.p_ <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (10m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  ggtitle("PAR (10m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR15.p_ <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (15m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  ggtitle("PAR (15m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))
PAR20.p_ <- map_base.gg + 
  geom_sf(data=filter(grid.long, predictor=="PAR (20m)"), colour=NA, aes(fill=value)) +
  pub_theme + scale_fill_viridis_c(expression(mu*plain(mol)/plain(m)^2/plain(day))) +
  ggtitle("PAR (20m)") + theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"))

fig_cov_a <- ggarrange(sst.p, fetch.p, fetchCat.p, nrow=1)
fig_cov_b <- ggarrange(PAR0.p, PAR2.p, PAR5.p, PAR10.p, PAR15.p, PAR20.p, nrow=2, ncol=3,
                       common.legend=T, legend="bottom")
fig_cov_c <- ggarrange(PAR0.p_, PAR2.p_, PAR5.p_, PAR10.p_, PAR15.p_, PAR20.p_, nrow=2, ncol=3,
                       common.legend=T, legend="bottom")
ggsave(glue("{fig.dir}Fig_A6a.png"), fig_cov_a, width=7, height=4, dpi=300)
ggsave(glue("{fig.dir}Fig_A6b.png"), fig_cov_b, width=7.5, height=8, dpi=300)
ggsave(glue("{fig.dir}Fig_A6c.png"), fig_cov_c, width=7.5, height=8, dpi=300)



# * Fig A.7 ---------------------------------------------------------------

fetch_thresh <- 4.02 # based on Pedersen sites
grid.i %>%
  mutate(fetchCat=factor(fetchCat, labels=c("Low", "High"))) %>%
  ggplot(aes(fetch, fill=fetchCat)) + 
  geom_vline(xintercept=fetch_thresh, linetype=2, colour="grey30") +
  geom_histogram(bins=80) +
  scale_fill_viridis_d("Exposure\ncategory", option="B", end=0.75) +
  pub_theme +
  labs(x=expression(log[10](fetch)))
ggsave(glue("{fig.dir}Fig_A7.png"), height=4, width=6, dpi=300)






# Supp B ------------------------------------------------------------------


# * Fig B.1 ---------------------------------------------------------------

figb1a <- mass.sum %>%
  filter(month==7) %>% filter(!stochParams) %>% 
  ungroup %>% 
  mutate(depth=paste0(depth, "m"),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20), "m"))) %>%
  ggplot(aes(biomass_mn, biomass_sd/biomass_mn, colour=fetch)) + 
  scale_colour_viridis_c(option="B") + 
  geom_point(shape=1, alpha=0.5) + 
  facet_grid(depth~.) +
  pub_theme +
  labs(x=bquote(Biomass~mean~(kg/m^2)), y="Biomass interannual CV")
figb1b <- mass.sum %>%
  filter(month==7) %>% filter(!stochParams) %>% 
  ungroup %>% 
  mutate(depth=paste0(depth, "m"),
         depth=factor(depth, levels=paste0(c(2,5,10,15,20), "m"))) %>%
  ggplot(aes(biomass_mn, biomass_sd, colour=fetch)) + 
  scale_colour_viridis_c(option="B") + 
  geom_point(shape=1, alpha=0.5) + 
  facet_grid(depth~.) +
  pub_theme +
  labs(x=bquote(Biomass~mean~(kg/m^2)), y=bquote(Biomass~interannual~sd~(kg/m^2)))
ggsave(glue("{fig.dir}Fig_B1.png"), 
       ggarrange(figb1a, figb1b, common.legend=T, legend="right", labels="auto"),
       width=8, height=9, dpi=300)




# * Fig B.2 ---------------------------------------------------------------

fig_b2a <- mass.sum %>% filter(month==7) %>%
  filter(!stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_mn, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day))),
       y="Mean July canopy biomass (kg)") + 
  pub_theme +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
fig_b2b <- mass.sum %>% filter(month==7) %>%
  filter(!stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  labs(y="SD July canopy biomass (kg)") + 
  pub_theme +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())
fig_b2c <- mass.sum %>% filter(month==7) %>%
  filter(!stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd/biomass_mn, colour=SST_mn)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("SST (°C)", option="C") +
  pub_theme +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day))),
       y="CV(July canopy biomass)")

fig_b2d <- mass.sum %>% filter(month==7) %>%
  filter(!stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_mn, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  pub_theme +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank())
fig_b2e <- mass.sum %>% filter(month==7) %>%
  filter(!stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  pub_theme +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank())
fig_b2f <- mass.sum %>% filter(month==7) %>%
  filter(!stochParams) %>%
  ggplot(aes(PAR_atDepth, biomass_sd/biomass_mn, colour=fetch)) +
  geom_point(shape=1) + 
  scale_colour_viridis_c("log(fetch)") +
  pub_theme +
  labs(x=expression(plain(PAR)~(mu*plain(mol)/plain(m)^2/plain(day)))) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

fig.sst <- ggarrange(plotlist=list(fig_b2a, fig_b2b, fig_b2c), nrow=3, 
                     common.legend=T, legend="bottom", 
                     labels=c("a.", "b.", "c."), label.x=0.85, label.y=0.95)
fig.fetch <- ggarrange(plotlist=list(fig_b2d, fig_b2e, fig_b2f), nrow=3, 
                       common.legend=T, legend="bottom", 
                       labels=c("d.", "e.", "f."), label.x=0.85, label.y=0.95)
fig.sst_fetch <- ggarrange(fig.sst, fig.fetch, widths=c(1, 0.9))
ggsave(glue("{fig.dir}Fig_B2.png"), 
       fig.sst_fetch, width=6, height=9, dpi=300)



# * Fig B.3 ---------------------------------------------------------------
fig_b3 <- map(c("recruits", "subcanopy", "canopy"), 
              ~readRDS(glue("out/processed/b_lagEffects_{.x}.rds")) %>%
                filter(depth %in% c(2,5,10,15)) %>%
                mutate(depth=factor(paste0(depth, "m"), levels=paste0(c(2,5,10,15,20), "m")),
                       covar=factor(covar, levels=c("strm", "PAR", "SST"),
                                    labels=c("Storms", "PAR", "SST"))) %>%
                ggplot(aes(lag, b_md, ymin=b_q10, ymax=b_q90, colour=covar, fill=covar)) +
                geom_hline(yintercept=0, colour="grey80") +
                geom_ribbon(alpha=0.3, colour=NA) +
                geom_line() +
                scale_colour_manual("", values=c("grey40", "goldenrod", "red4")) +
                scale_fill_manual("", values=c("grey40", "goldenrod", "red4")) +
                facet_grid(depth~.) + 
                scale_x_continuous(breaks=c(5, 10, 15)) +
                pub_theme + theme(legend.position="bottom") +
                labs(x="Lag (years)", 
                     y=glue("Effect on {.x}\n(standardized slope: median + middle 80%)")))
ggsave(glue("{fig.dir}Fig_B3.png"),
       ggarrange(plotlist=fig_b3, ncol=3, common.legend=T, labels="auto", legend="bottom"),
       width=9, height=6, dpi=300)



# * Fig B.4 ---------------------------------------------------------------

fig_b4 <- vector("list")
for(i in c("recruits", "subcanopy", "canopy")) {
 b.sum.id_jul <- readRDS(glue("out/processed/b_lagEffects_id_{i}.rds"))
 b_rng <- range(b.sum.id_jul$b_md)
 
 b.amplitude <- b.sum.id_jul %>%
   filter(depth < 15) %>%
   group_by(covar, id, depth) %>%
   summarise(b_amplitude=max(b_mn)-min(b_mn)) %>%
   mutate(depth=factor(paste0(depth, "m"), 
                       levels=paste0(c(2,5,10,15,20), "m")),
          covar=factor(covar, levels=c("strm", "PAR", "SST"),
                       labels=c("Storms", "PAR", "SST"))) %>%
   full_join(grid.sf, ., by="id")
 
 fig_b4[[i]] <- map_base.gg +
   geom_sf(data=b.amplitude, aes(fill=b_amplitude), colour=NA) +
   scale_fill_viridis_c("Effect amplitude", option="inferno", limits=c(0, NA),
                        guide=guide_colorbar(title.position="top")) + 
   facet_grid(depth~covar) +
   pub_theme +
   theme(legend.position="bottom", 
         legend.key.height=unit(0.2, "cm"),
         legend.key.width=unit(1.2, "cm"),
         title=element_text(size=9))
}
ggsave(glue("{fig.dir}Fig_B4.png"),
       ggarrange(plotlist=fig_b4, ncol=3, common.legend=T, labels="auto", legend="bottom"),
       width=12, height=6, dpi=300)




# * Fig B.5 ---------------------------------------------------------------

ri.df %>%
  filter(var != "growStipe.canopy") %>%
  filter(month=="july") %>%
  right_join(grid.i %>% select(id, fetchCat), by="id") %>%
  mutate(rel.inf=rel.inf*100) %>%
  group_by(var, month, depth, response, fetchCat) %>%
  summarise(md=median(rel.inf),
            q10=quantile(rel.inf, probs=0.1),
            q25=quantile(rel.inf, probs=0.25),
            q75=quantile(rel.inf, probs=0.75),
            q90=quantile(rel.inf, probs=0.9)) %>%
  ungroup %>%
  mutate(response=factor(response, levels=c("biomass_mn", "biomass_sd"),
                         labels=c("Mean", "Interannual sd")),
         depth=as.numeric(str_remove(depth, "m")),
         param=sens.param$eq[match(var, sens.param$messy)],
         fetchCat=factor(fetchCat, labels=c("Low", "High"))) %>% 
  ggplot(aes(depth, md, colour=fetchCat, fill=fetchCat)) +
  geom_hline(yintercept=0, colour="grey90") +
  geom_ribbon(aes(ymin=q25, ymax=q75), colour=NA, alpha=0.4) + 
  geom_ribbon(aes(ymin=q10, ymax=q90), colour=NA, alpha=0.3) + 
  geom_point(position=position_dodge(width=0.5)) +
  geom_line() +
  facet_grid(response~param, 
             labeller=labeller(response=label_value, depth=label_value, param=label_parsed)) + 
  scale_colour_viridis_d("Exposure", option="B", end=0.75) +
  scale_fill_viridis_d("Exposure", option="B", end=0.75) +
  scale_x_continuous(breaks=c(0, 5, 10, 15), limits=c(0,16)) +
  scale_y_continuous(breaks=c(0, 25, 50)) + 
  pub_theme +
  theme(legend.position="bottom", legend.key.height=unit(0.2, "cm"),
        legend.title=element_text(size=11),
        axis.text.y=element_text(size=12), 
        strip.text=element_text(size=11), 
        plot.title=element_text(size=12)) +
  labs(y="Relative influence (%)", x="Depth (m)", 
       title="Canopy biomass sensitivity")
ggsave(glue("{fig.dir}Fig_B5.png"), width=11, height=4, dpi=300)

# * Fig B.6 ---------------------------------------------------------------

bind_rows(
  obs.ls$N_canopy.lm %>%
    group_by(id, PAR_atDepth) %>%
    summarise(N_mn=mean(log(N+1)),
              N_sd=sd(log(N+1))) %>%
    ungroup %>%
    mutate(type="Observed"),
  pop.sum_canopy %>% ungroup %>%
    filter(month=="July" & id %in% obs.ls$N_canopy.lm$id) %>%
    select(PAR_atDepth_mn, N_mn, logN_sd) %>%
    rename(PAR_atDepth=PAR_atDepth_mn,
           N_sd=logN_sd) %>%
    mutate(N_mn=log(N_mn+1)) %>%
    mutate(type="Modelled")
) %>%
  ggplot(aes(PAR_atDepth, N_mn, ymin=N_mn-N_sd, ymax=N_mn+N_sd, colour=type)) +
  geom_point(shape=1) + 
  geom_linerange(size=0.1) +
  stat_smooth(method="lm", se=F, linetype=2) +
  scale_colour_manual("", values=c("gray30", "red")) +
  labs(x=expression(PAR~(mu*plain(mol)/plain(m)^2/plain(day))),
       y=expression(Mean~ln~canopy~density~(N/m^2))) +
  guides(colour=guide_legend(override.aes=list(linetype=1))) +
  pub_theme +
  theme(legend.position=c(0.85, 0.15),
        legend.background=element_blank())
ggsave(glue("{fig.dir}Fig_B6.png"), width=4, height=4, dpi=300)






