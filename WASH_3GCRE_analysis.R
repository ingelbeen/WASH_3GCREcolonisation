##########################################################################
# WASH & 3GCRE colonisation                                              #
# Figures 1 to 5 and Table 1 — association WASH & 3GCR-E colonisation    #
##########################################################################
# install/load packages
pacman::p_load(readxl, dplyr, tidyr, stringr, tibble, purrr, ggplot2, patchwork, geepack, brms, posterior, ggnewscale, here)

# import data
setwd("C:/Users/binge/Documents/community-level AMR/scripts")
raw <- read_excel("ESBL_CRE_FecalColonization_HealthyIndividuals_GlobalEvidence_v12.xlsx", sheet = "Primary Studies") |> rename_with(str_squish)

# prep data
num <- function(x) suppressWarnings(as.numeric(x))
max_year <- function(s){ y <- as.integer(unlist(regmatches(s, gregexpr("(19|20)\\d{2}", s)))); if(length(y)) max(y) else NA_real_ }
cont_map <- c(`Eastern Africa`="Africa",`Western Africa`="Africa",`Middle Africa`="Africa",`Southern Africa`="Africa",`Northern Africa`="Africa",
              `Western Europe`="Europe",`Northern Europe`="Europe",`Southern Europe`="Europe",`Eastern Europe`="Europe",
              `Eastern Asia`="Asia",`South-eastern Asia`="Asia",`Southern Asia`="Asia",`Western Asia`="Asia",`Central Asia`="Asia",
              `South America`="Americas",`Central America`="Americas",`Northern America`="Americas",`Caribbean`="Americas",
              `Australia and New Zealand`="Oceania")

dat <- raw |>
  mutate(
    N      = num(`N Stool/Swab Collected`),
    npos   = num(`N Positive ESBL/3GCR`),
    # recover 3 multi-country arms (Parra et al.) with per-country prevalence but no count
    npos   = if_else(`Ref #` %in% c("57a","57c","57e"),
                     round(N * num(`Prevalence ESBL (%)`) / 100), npos),
    excluded = !is.na(`Excluded? (Reason)`),
    year   = map_dbl(`Year(s) Sample Collection`, max_year),
    year   = if_else(is.na(year), num(`Year Published`), year),
    year_c = year - 2015,
    amu    = num(`Antibiotic use (DDD/1000/day)`),
    w_sm=num(`Safely Managed Water (granular, stratum-matched) (%)`),
    w_bas=num(`At Least Basic Water (granular, stratum-matched) (%)`),
    w_lim=num(`Limited Water Service (>30 min) (%)`),
    w_surf=num(`Surface Water (No Service) (%)`), w_unimp=num(`Unimproved Water (%)`),
    w_basic=w_bas - w_sm,
    s_basal=num(`At Least Basic Sanitation (%) [ladder]`),
    s_lim=num(`Limited Sanitation (shared/unimproved) (%)`),
    s_unimp=num(`Unimproved Sanitation (%)`), s_ns=num(`No Sanitation Service (open defecation) (%)`),
    san_sm=num(`Safely managed sanitation (%)`),
    subreg = `UN M49 Subregion`,
    continent = unname(cont_map[subreg]),
    npos_cre=num(`N Positive CRE/CPE`), npos_ec=num(`N Positive ESBL E. coli`),
    w_ladder = (w_unimp + 2*w_lim + 3*w_basic + 4*w_sm)/(w_surf+w_unimp+w_lim+w_basic+w_sm),
    s_ladder = (s_unimp + 2*s_lim + 3*s_basal)/(s_ns+s_unimp+s_lim+s_basal))

# analysis set: drop flagged, invalid counts/year, and traveler cohorts
ana <- dat |>
  filter(!excluded, !is.na(N), !is.na(npos), npos<=N, !is.na(year),
         !str_detect(coalesce(Country,""), regex("travell|abroad", ignore_case=TRUE))) |>
  mutate(nneg = N - npos, prev = npos/N)
rare <- ana |> count(subreg) |> filter(n < 4) |> pull(subreg)
ana  <- ana |> mutate(subreg_g = if_else(is.na(subreg) | subreg %in% rare, "Other", subreg))
cat(sprintf("Analysis set: %d studies, %d participants, years %d-%d\n",
            nrow(ana), sum(ana$N), min(ana$year), max(ana$year)))

# harmonise colour palette between Figures
# continent hue families: light→mid = dots, dark = trend lines
fam <- list(
  Africa   = c(light="#FEC44F", mid="#EC7014", dark="#8C2D04"),
  Asia     = c(light="#FC9272", mid="#DE2D26", dark="#67000D"),
  Europe   = c(light="#A1D99B", mid="#41AB5D", dark="#00441B"),
  Americas = c(light="#9ECAE1", mid="#4292C6", dark="#08306B"),
  Oceania  = c(light="#BCBDDC", mid="#807DBA", dark="#3F007D"))
cont_cols    <- sapply(fam, `[[`, "dark")   # continent trend lines (Fig 3)
cont_cols_pt <- sapply(fam, `[[`, "mid")    # continent-coloured points (Fig 4C)
cont_lty     <- c(Africa="solid", Asia="dashed", Europe="dotdash", Americas="twodash", Oceania="longdash")

# subregion palette: light→mid shades within each continent family (matches Fig 3 dots)
subreg_levels <- sort(unique(na.omit(ana$subreg)))
subreg_cols <- unlist(lapply(names(fam), function(cc){
  srs <- sort(subreg_levels[unname(cont_map[subreg_levels]) == cc])
  if (!length(srs)) return(NULL)
  setNames(colorRampPalette(c(fam[[cc]]["light"], fam[[cc]]["mid"]))(length(srs)), srs)}))

# JMP ladder order (lowest → highest service)
water_order <- c("w_surf","w_unimp","w_lim","w_basic","w_sm")
san_order   <- c("s_ns","s_unimp","s_lim","s_basal")

# JMP ladder-rung colours (shared by Fig 2 stacked bars and Fig 4 violins)
water_rung_cols <- c(w_sm="#08589e", w_basic="#7bccc4", w_lim="#fee391", w_unimp="#fe9929", w_surf="#cc4c02")
san_rung_cols   <- c(s_basal="#7bccc4", s_lim="#fee391", s_unimp="#fe9929", s_ns="#bd0026")
water_rung_lab  <- c(w_sm="Safely managed", w_basic="Basic", w_lim="Limited", w_unimp="Unimproved", w_surf="Surface")
san_rung_lab    <- c(s_basal="At-least-basic", s_lim="Limited", s_unimp="Unimproved", s_ns="Open defecation")

# effect colours (Fig 5 forest): highlighted when 95% CrI excludes 1
effect_cols <- c(`TRUE`="#01665E", `FALSE`="grey65")

# beta-binomial GLMM helper (Bayesian MCMC via brms/Stan) + OR extractor
bb <- beta_binomial()
glmm <- function(rhs, data=ana, chains=4, iter=2000)
  brm(bf(as.formula(paste("npos | trials(N) ~", rhs))), data=data, family=bb,
      chains=chains, iter=iter, cores=chains, seed=1, refresh=0,
      control=list(adapt_delta=0.95))
OR_fix <- function(fit, p, sc=1) exp(sc*fixef(fit)[p, c("Estimate","Q2.5","Q97.5")])

#### FIGURE 2 — JMP water & sanitation ladder coverage, by UN subregion ####
keep <- ana |> count(subreg) |> filter(n >= 3) |> pull(subreg)
ord  <- ana |> filter(subreg %in% keep) |> group_by(subreg) |>
        summarise(m = median(prev)) |> arrange(m) |> pull(subreg)
ladder_cov <- function(cols)
  ana |> filter(subreg %in% keep) |> group_by(subreg) |>
    summarise(across(all_of(cols), ~mean(.x, na.rm=TRUE)), .groups="drop") |>
    rowwise() |> mutate(across(all_of(cols), ~100*.x/sum(c_across(all_of(cols))))) |> ungroup() |>
    pivot_longer(all_of(cols), names_to="rung", values_to="pct")
water_cov <- ladder_cov(c("w_sm","w_basic","w_lim","w_unimp","w_surf"))
san_cov   <- ladder_cov(c("s_basal","s_lim","s_unimp","s_ns"))
fig2A <- ggplot(water_cov, aes(pct, factor(subreg, ord), fill=factor(rung, water_order))) + geom_col() +
  scale_fill_manual(values=water_rung_cols, labels=water_rung_lab, breaks=water_order, name=NULL,
                    guide=guide_legend(reverse=TRUE)) +
  labs(title="A  Water ladder", x="Population (%)", y=NULL)
fig2B <- ggplot(san_cov, aes(pct, factor(subreg, ord), fill=factor(rung, san_order))) + geom_col() +
  scale_fill_manual(values=san_rung_cols, labels=san_rung_lab, breaks=san_order, name=NULL,
                    guide=guide_legend(reverse=TRUE)) +
  labs(title="B  Sanitation ladder", x="Population (%)", y=NULL)
figure2 <- fig2A + fig2B
ggsave("figure2.jpeg", plot = figure2, dpi = 300, height = 5, width = 10)

#### TABLE 1 — time-standardised (2015) and projected (2025) prevalence ####
ana_t <- ana |> filter(subreg %in% keep, !is.na(year_c), !is.na(continent)) |>
         mutate(subreg=factor(subreg), continent=factor(continent)) |> arrange(Country)
m2015 <- geeglm(cbind(npos,nneg) ~ subreg + year_c,             id=Country, data=ana_t, family=binomial, corstr="exchangeable")
mproj <- geeglm(cbind(npos,nneg) ~ subreg + year_c:continent,   id=Country, data=ana_t, family=binomial, corstr="exchangeable")
pci <- function(fit, newdata){                                   # delta-method CI on prevalence scale
  X   <- model.matrix(delete.response(terms(fit)), newdata, xlev=fit$xlevels)
  eta <- as.numeric(X %*% coef(fit)); se <- sqrt(diag(X %*% vcov(fit) %*% t(X)))
  round(100*plogis(eta + c(0,-1.96,1.96)*se), 1)
}
sub_cont <- ana_t |> distinct(subreg, continent) |> deframe()
table1 <- map_dfr(levels(ana_t$subreg), function(s){
  nd15 <- tibble(subreg=factor(s, levels(ana_t$subreg)), year_c=0)
  nd25 <- tibble(subreg=factor(s, levels(ana_t$subreg)),
                 continent=factor(sub_cont[[s]], levels(ana_t$continent)), year_c=10)
  a <- pci(m2015, nd15); p <- pci(mproj, nd25)
  tibble(subreg=s, n=sum(ana_t$subreg==s),
         p2015=a[1], lo2015=a[2], hi2015=a[3], p2025=p[1], lo2025=p[2], hi2025=p[3])
}) |> arrange(desc(p2015))
print(table1)

#### FIGURE 3 — 3GCR-E prevalence (A,B) and CRE (C) temporal trend (year) ####
# continent hue families: light→mid = subregion dots, dark = trend line 
fam <- list(
  Africa   = c(light="#FEC44F", mid="#EC7014", dark="#8C2D04"),  # oranges
  Asia     = c(light="#FC9272", mid="#DE2D26", dark="#67000D"),  # reds
  Europe   = c(light="#A1D99B", mid="#41AB5D", dark="#00441B"),  # greens
  Americas = c(light="#9ECAE1", mid="#4292C6", dark="#08306B"),  # blues
  Oceania  = c(light="#BCBDDC", mid="#807DBA", dark="#3F007D"))   # purples
cont_cols <- sapply(fam, `[[`, "dark")                            # trend-line colours
cont_lty  <- c(Africa="solid", Asia="dashed", Europe="dotdash",
               Americas="twodash", Oceania="longdash")

# subregion dot palette: distinguishable shades within each continent family
subreg_levels <- sort(unique(na.omit(ana$subreg)))
subreg_cols <- unlist(lapply(names(fam), function(cc){
  srs <- sort(subreg_levels[unname(cont_map[subreg_levels]) == cc])
  if (!length(srs)) return(NULL)
  setNames(colorRampPalette(c(fam[[cc]]["light"], fam[[cc]]["mid"]))(length(srs)), srs)
}))

# Panel A — 3GCR prevalence over time, each dot is a study in a colour by subregion with trendlines by conttinent
fig3A <- ggplot(filter(ana, !is.na(year)), aes(year, 100*prev)) +
  geom_point(aes(colour = subreg, size = sqrt(N)), alpha = 0.75) +
  scale_colour_manual(values = subreg_cols, limits = names(subreg_cols),
                      drop = FALSE, name = "Subregion", na.translate = FALSE) +
  scale_size_continuous(guide = "none") +
  new_scale_colour() +
  geom_line(data = cont_pred, aes(year, prev, colour = continent, linetype = continent),
            linewidth = 1, inherit.aes = FALSE) +
  scale_colour_manual(values = cont_cols, name = "Continent trend") +
  scale_linetype_manual(values = cont_lty, name = "Continent trend") +
  labs(title = "A  3GCR-E colonisation vs year",
       x = "Year of sample collection", y = "3GCR-E prevalence (%)")

# Panel B — boxplot, same subregion colours (no legend)
fig3B <- ggplot(filter(ana, subreg %in% keep),
                aes(100*prev, factor(subreg, ord), fill = subreg)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = subreg_cols) +
  labs(title = "B  Prevalence by UN M49 subregion", x = "3GCR-E prevalence (%)", y = NULL)

# Panel C — CR colonization by year, dots = subregion (same palette), lines = continent 
fig3C <- ggplot(filter(ana, !is.na(npos_cre), !is.na(year)), aes(year, 100*npos_cre/N)) +
  geom_point(aes(colour = subreg, size = sqrt(N)), alpha = 0.75) +
  scale_colour_manual(values = subreg_cols, limits = names(subreg_cols),
                      drop = FALSE, name = "Subregion", na.translate = FALSE) +
  scale_size_continuous(guide = "none") +
  new_scale_colour() +
  geom_smooth(aes(colour = continent, linetype = continent),
              method = "lm", se = FALSE, linewidth = 1) +
  scale_colour_manual(values = cont_cols, name = "Continent trend") +
  scale_linetype_manual(values = cont_lty, name = "Continent trend") +
  coord_cartesian(ylim = c(0, 15)) +
  labs(title = "C  Carbapenem-resistant Enterobacterales vs year",
       x = "Year of sample collection", y = "CRE prevalence (%)") + 
  theme(legend.position = "none")

figure3 <- fig3A + fig3B + fig3C + plot_layout(nrow = 1, guides = "collect")
figure3
ggsave("figure3.jpeg", plot = figure3, dpi = 300, height = 7, width = 16)

#### FIGURE 4 — distribution across WASH ladders and antibiotic use ####
r_ladders <- cor(ana$w_ladder, ana$s_ladder, use="complete.obs")
cat(sprintf("[Figure 4] ladder-score correlation r = %.2f\n", r_ladders))

# dominant rung per study (ties broken to first) for panels A/B
ana <- ana |> mutate(
  dom_water = if_else(rowSums(is.na(across(all_of(water_order))))==0,
                      water_order[max.col(across(all_of(water_order)), ties.method="first")], NA_character_),
  dom_san   = if_else(rowSums(is.na(across(all_of(san_order))))==0,
                      san_order[max.col(across(all_of(san_order)), ties.method="first")], NA_character_))

# A — colonisation by dominant water rung (violins, ladder order low→high)
fig4A <- ggplot(filter(ana, !is.na(dom_water)),
                aes(factor(dom_water, water_order), 100*prev, fill=dom_water)) +
  geom_violin(scale="width", alpha=0.75, colour=NA) +
  geom_jitter(width=0.12, size=0.7, alpha=0.5) +
  scale_fill_manual(values=water_rung_cols, guide="none") +
  scale_x_discrete(limits=water_order, labels=water_rung_lab) +
  labs(title="A  By dominant water rung", x=NULL, y="3GCR-E prevalence (%)")

# B — colonisation by dominant sanitation rung
fig4B <- ggplot(filter(ana, !is.na(dom_san)),
                aes(factor(dom_san, san_order), 100*prev, fill=dom_san)) +
  geom_violin(scale="width", alpha=0.75, colour=NA) +
  geom_jitter(width=0.12, size=0.7, alpha=0.5) +
  scale_fill_manual(values=san_rung_cols, guide="none") +
  scale_x_discrete(limits=san_order, labels=san_rung_lab) +
  labs(title="B  By dominant sanitation rung", x=NULL, y="3GCR-E prevalence (%)")

# C — colonisation vs antibiotic use, coloured by continent
fig4C <- ggplot(filter(ana, !is.na(amu)), aes(amu, 100*prev)) +
  geom_point(aes(colour=continent, size=sqrt(N)), alpha=0.75) +
  geom_smooth(method="lm", se=FALSE, colour="black", linetype=2, linewidth=0.8) +
  scale_colour_manual(values=cont_cols_pt, name="Continent") +
  scale_size_continuous(guide="none") +
  theme(legend.position = "bottom") +
  labs(title="C  Colonisation vs antibiotic use", x="Antibiotic use (DDD/1000/day)", y="3GCR-E prevalence (%)")

# D — water vs sanitation ladder jointly, point colour = prevalence, size = √N
fig4D <- ggplot(filter(ana, !is.na(w_ladder), !is.na(s_ladder)), aes(w_ladder, s_ladder)) +
  geom_point(aes(colour=100*prev, size=sqrt(N)), alpha=0.85) +
  scale_colour_viridis_c(option="magma", direction=-1, name="Prevalence (%)") +
  scale_size_continuous(guide="none") +
  labs(title="D  Water vs sanitation ladder jointly",
       x="Water-ladder score (0–4)", y="Sanitation-ladder score (0–3)") +
  theme(legend.position = "bottom") 

figure4 <- (fig4A + fig4B) / (fig4C + fig4D)
figure4
ggsave("figure4.jpeg", plot = figure4, dpi = 300, height = 10, width = 10)

#### FIGURE 5 — adjusted associations (beta-binomial GLMMs) ####
# fits 
fit_water        <- glmm("w_ladder + amu + year_c + (1|subreg_g) + (1|Country)")
fit_water_sanadj <- glmm("w_ladder + s_ladder + amu + year_c + (1|subreg_g) + (1|Country)")
fit_san          <- glmm("s_ladder + amu + year_c + (1|subreg_g) + (1|Country)")
fit_amu          <- glmm("amu + year_c + (1|subreg_g) + (1|Country)")
fit_san_sm       <- glmm("san_sm + amu + year_c + (1|subreg_g) + (1|Country)")

res <- list(
  water        = OR_fix(fit_water,        "w_ladder"),
  water_sanadj = OR_fix(fit_water_sanadj, "w_ladder"),
  sanitation   = OR_fix(fit_san,          "s_ladder"),
  amu          = OR_fix(fit_amu,          "amu", sc = 10),
  san_sm       = OR_fix(fit_san_sm,       "san_sm", sc = 10))

stopifnot(is.list(res)); str(res)   # confirm: $water Named num [1:3], etc.

# A — adjusted associations forest
fig5A <- tibble(
  term = c("Water ladder (adj. sanitation, AMU, year)","Water ladder (adj. AMU, year)",
           "Sanitation ladder (adj. AMU, year)","Antibiotic use (per 10 DDD)"),
  est = c(res$water_sanadj[["Estimate"]], res$water[["Estimate"]], res$sanitation[["Estimate"]], res$amu[["Estimate"]]),
  lo  = c(res$water_sanadj[["Q2.5"]],     res$water[["Q2.5"]],     res$sanitation[["Q2.5"]],     res$amu[["Q2.5"]]),
  hi  = c(res$water_sanadj[["Q97.5"]],    res$water[["Q97.5"]],    res$sanitation[["Q97.5"]],    res$amu[["Q97.5"]])) |>
  mutate(sig = lo>1 | hi<1, term = factor(term, levels=rev(term))) |>
  ggplot(aes(est, term, colour=sig)) +
  geom_vline(xintercept=1, linetype=2, colour="grey50") +
  geom_pointrange(aes(xmin=lo, xmax=hi), linewidth=0.9, fatten=3) +
  scale_x_log10() + scale_colour_manual(values=effect_cols, guide="none") +
  labs(title="A  Adjusted associations", x="OR for 3GCR-E colonisation (log scale)", y=NULL)

# if that gave an error, then
# rebuild the forest input directly from the fits (order-safe)
# orx <- function(fit, p, sc=1){
#   v <- exp(sc * fixef(fit)[p, c("Estimate","Q2.5","Q97.5")])
#   tibble(est = v[["Estimate"]], lo = v[["Q2.5"]], hi = v[["Q97.5"]])
# }
# 
# forestA <- bind_rows(
#   cbind(term = "Water ladder (adj. sanitation, AMU, year)", orx(fit_water_sanadj, "w_ladder")),
#   cbind(term = "Water ladder (adj. AMU, year)",             orx(fit_water,        "w_ladder")),
#   cbind(term = "Sanitation ladder (adj. AMU, year)",        orx(fit_san,          "s_ladder")),
#   cbind(term = "Antibiotic use (per 10 DDD)",               orx(fit_amu,          "amu", sc = 10))
# ) |>
#   mutate(sig = lo > 1 | hi < 1, term = factor(term, levels = rev(term)))
# 
# fig5A <- ggplot(forestA, aes(est, term, colour = sig)) +
#   geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
#   geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = 0.9, fatten = 3) +
#   scale_x_log10() + scale_colour_manual(values = effect_cols, guide = "none") +
#   labs(title = "A  Adjusted associations", x = "OR for 3GCR-E colonisation (log scale)", y = NULL)

# B — compositional substitution, OR per +10% up each adjacent rung (top rung as reference)
dw <- as_draws_df(glmm("w_surf + w_unimp + w_lim + w_basic + s_ladder + amu + year_c + (1|subreg_g) + (1|Country)"))
ds <- as_draws_df(glmm("s_ns + s_unimp + s_lim + amu + year_c + (1|subreg_g) + (1|Country)"))
mk <- function(label, x){ q <- quantile(x, c(.025,.975)); tibble(step=label, est=exp(mean(x)), lo=exp(q[[1]]), hi=exp(q[[2]])) }
# listed in JMP ladder order (lowest step first)
subst <- bind_rows(
  mk("Water: surface → unimproved",          10*(dw$b_w_unimp - dw$b_w_surf)),
  mk("Water: unimproved → limited",          10*(dw$b_w_lim   - dw$b_w_unimp)),
  mk("Water: limited → basic",               10*(dw$b_w_basic - dw$b_w_lim)),
  mk("Water: basic → safely managed",        10*(0            - dw$b_w_basic)),
  mk("Sanitation: open defec. → unimproved", 10*(ds$b_s_unimp - ds$b_s_ns)),
  mk("Sanitation: unimproved → limited",     10*(ds$b_s_lim   - ds$b_s_unimp)),
  mk("Sanitation: limited → at-least-basic", 10*(0            - ds$b_s_lim))) |>
  mutate(sig = lo>1 | hi<1,
         step = factor(step, levels=rev(step)))   # rev() → lowest step at top, top-of-ladder at bottom

fig5B <- ggplot(subst, aes(est, step, colour=sig)) +
  geom_vline(xintercept=1, linetype=2, colour="grey50") +
  geom_pointrange(aes(xmin=lo, xmax=hi), linewidth=0.9, fatten=3) +
  scale_x_log10() + scale_colour_manual(values=effect_cols, guide="none") +
  labs(title="B  Substitution: OR per +10% up one rung", x="OR (log scale)", y=NULL)

figure5 <- fig5A + fig5B
figure5
ggsave("figure5.jpeg", plot = figure5, dpi = 300, height = 3.5, width = 10)

#### SENSITIVITY ANALYSES ####
sens_fit <- function(data, var)
  OR_fix(glmm(paste(var, "+ amu + year_c + (1|subreg_g) + (1|Country)"), data=data, chains=2, iter=1600), var)
# (1) 3GCR-E. coli outcome
ec   <- ana |> filter(!is.na(npos_ec)) |> mutate(npos=npos_ec) |> filter(npos<=N)
# (2) LMIC only
HIC  <- c("Australia","Belgium","Cyprus","Denmark","France","Germany","Hungary","Japan","Netherlands",
          "Norway","Portugal","Qatar","Republic of Korea","Saudi Arabia","Singapore","Spain","Sweden",
          "Switzerland","United States of America","USA","French Guiana")
lmic <- ana |> filter(!Country %in% HIC)
# (3) leave-out flagged outliers (Datta #36, Moremi #113)
noout<- ana |> filter(!`Ref #` %in% c("36","113"))
sensitivity <- list(
  EC_water        = sens_fit(ec,   "w_ladder"), EC_sanitation      = sens_fit(ec,   "s_ladder"),
  LMIC_water      = sens_fit(lmic, "w_ladder"), LMIC_sanitation    = sens_fit(lmic, "s_ladder"),
  outlier_water   = sens_fit(noout,"w_ladder"), outlier_sanitation = sens_fit(noout,"s_ladder"))
print(lapply(sensitivity, round, 2))
