################################################################################
# Access to safe water or sanitation and gut colonisation with 3GC-resistant   #
# Enterobacterales (3GCR-E): global ecological analysis in healthy populations #
################################################################################

# Last edited by Brecht Ingelbeen, 13/10/2026
# script follows the order of the manuscript, but does not contain Figure 1 (created in ppt)

### 0. PACKAGES AND IMPORT DATA ####
# packages
pacman::p_load(readxl, dplyr, tidyr, stringr, purrr, tibble, ggplot2, patchwork,
               ggnewscale, scales, geepack, brms, posterior, readr)
# import data
raw <- read_excel("ESBL_CRE_FecalColonization_HealthyIndividuals_GlobalEvidence_v16.xlsx", sheet = "Primary Studies") |> rename_with(str_squish)
jmp_w <- read_excel("JMP_2025_WLD.xlsx", sheet = "wat")
jmp_s <- read_excel("JMP_2025_WLD.xlsx", sheet = "san") # JMP only needed for the WASH time-trend sentence (section 8), individual matching JMP data is already added in the study spreadsheet

# run settings for MCMC
QUICK    <- FALSE    # TRUE = 3 imputations, 1000 iterations: fast smoke test only
M_IMP    <- if (QUICK) 3 else 20     # number of imputed datasets (manuscript: 20)
CHAINS   <- 4                        # manuscript: four chains
ITER     <- if (QUICK) 1000 else 2000
SEED     <- 2026
RUN_LMIC <- FALSE    # section 15 (supplementary, not in the text)
set.seed(SEED)

# functions
fmt  <- function(x, d = 2) formatC(x, format = "f", digits = d)
fmt_or <- function(est, lo, hi, d = 2) sprintf("%s (%s-%s)", fmt(est, d), fmt(lo, d), fmt(hi, d))
num  <- function(x) suppressWarnings(as.numeric(x))
say  <- function(...) cat(sprintf(...), "\n")
hdr  <- function(x) cat("\n", strrep("=", 70), "\n", x, "\n", strrep("=", 70), "\n", sep = "")

### 1. STUDY IDENTIFICATION AND INCLUSION ####
# exclude studies from the main 3GCR analyses, after full text screening
# ref 42 (Yen et al), 53 (Arum et al) and 50 (Li et al) only report CRE; ref 175 (Hashemizadeh et al) reports on patients with digestive symptoms
raw <- raw |> filter(`Include in Analysis` == "Include")

say("Include-tagged rows in database: %d", nrow(raw))

# if data collection over multiple years, pick the latest one
max_year <- function(s) {
  y <- as.integer(unlist(regmatches(s, gregexpr("(19|20)\\d{2}", s))))
  if (length(y)) max(y) else NA_real_}

# assign continents to UN subregions (which have been added to the database)
cont_map <- c(`Eastern Africa`="Africa", `Western Africa`="Africa", `Middle Africa`="Africa",
              `Southern Africa`="Africa", `Northern Africa`="Africa",
              `Western Europe`="Europe", `Northern Europe`="Europe", `Southern Europe`="Europe",
              `Eastern Europe`="Europe",
              `Eastern Asia`="Asia", `South-eastern Asia`="Asia", `Southern Asia`="Asia",
              `Western Asia`="Asia", `Central Asia`="Asia",
              `South America`="Americas", `Central America`="Americas", `Northern America`="Americas",
              `Caribbean`="Americas", `Australia and New Zealand`="Oceania")

# clean data
dat <- raw |>  mutate(
    ref      = as.character(`Ref #`),
    study_id = str_extract(ref, "^\\d+"),          # sub-rows 57a, 57c ... belong to one study
    N        = num(`N Stool/Swab Collected`),
    npos     = num(`N Positive ESBL/3GCR`),    # checked, no studies with only percentages reported
    npos_ec  = num(`N Positive ESBL E. coli`),
    npos_cre = num(`N Positive CRE/CPE`),
    prev_pct_num = num(str_remove(`Prevalence ESBL (%)`, "%")),
    year     = map_dbl(`Year(s) Sample Collection`, max_year),
    year     = if_else(is.na(year), num(`Year Published`), year),
    year_c   = year - 2015,
    amu      = num(`Antibiotic use (DDD/1000/day)`),
    gdp      = num(`GDP per capita (current int$, matched year)`),
    # ---- WATER ladder (surface 0, unimproved 1, limited 2, basic 3, safely managed 4)
    w_sm_db  = num(`Safely Managed Water (granular, stratum-matched) (%)`),   # contains the imputed values
    w_sm_imp = !is.na(num(`Safely Managed Water - imputation SD (%)`)),       # TRUE = value was imputed
    w_sm_obs = if_else(w_sm_imp, NA_real_, w_sm_db),                          # observed JMP value only
    w_bas    = num(`At Least Basic Water (granular, stratum-matched) (%)`),   # cumulative (basic + SM)
    w_pip    = num(`Piped Water Supply (%)`),
    w_lim    = num(`Limited Water Service (>30 min) (%)`),
    w_unimp  = num(`Unimproved Water (%)`),
    w_surf   = num(`Surface Water (No Service) (%)`),
    # ---- SANITATION ladder (open defecation 0, unimproved 1, limited 2, basic 3, safely managed 4)
    s_sm     = num(`Safely managed sanitation (%)`),
    s_basic  = num(`Basic only sanitation (%)`),
    wash_ok  = !is.na(w_bas) & !is.na(s_basic),               # FALSE for Taiwan (Refs 135, 136): no JMP series
    s_lim    = num(`Limited Sanitation (shared/unimproved) (%)`),
    s_unimp  = num(`Unimproved Sanitation (%)`),
    s_ns     = num(`No Sanitation Service (open defecation) (%)`),
    subreg   = `UN M49 Subregion`,
    continent = unname(cont_map[subreg]),
    # normalise country names (strip "(region)" qualifiers; unify USA) so one country = one cluster
    Country  = str_squish(str_remove(Country, "\\s*\\(.*?\\)\\s*")),
    Country  = if_else(Country == "USA", "United States of America", Country),
    cid      = as.integer(factor(Country)),   # numeric country id for geeglm
    setting  = factor(case_when(str_to_lower(str_squish(`Urban/Rural Setting`)) == "urban" ~ "Urban",
                                str_to_lower(str_squish(`Urban/Rural Setting`)) == "rural" ~ "Rural",
                                TRUE ~ "Unspecified"),
                      levels = c("Urban", "Rural", "Unspecified")),
    .row = row_number())

# keep only observations with complete data for modelling: N, year present, counts <= N (outcome present for all-species OR E. coli)
dat_ok <- dat |>  filter(!is.na(N), !is.na(year),
         (!is.na(npos) & npos <= N) | (!is.na(npos_ec) & npos_ec <= N)) #none dropped
dat_w <- dat_ok |> filter(wash_ok)          # rows that can enter any WASH analysis / imputation
message("Rows without JMP WASH linkage (kept for prevalence/time trend only): ",
        paste(dat_ok$ref[!dat_ok$wash_ok], collapse = ", "))

# primary analysis set (all-species 3GCR-E)
ana <- dat_ok |> filter(!is.na(npos), npos <= N) |> mutate(nneg = N - npos, prev = npos / N)
# create a new var bringing together infrequent subregions so a single random intercept can be used, as a precaution for model stability
rare <- ana |> count(subreg) |> filter(n < 4) |> pull(subreg)
ana  <- ana |> mutate(subreg_g = if_else(is.na(subreg) | subreg %in% rare, "Other", subreg))

# E. coli-specific set (built independently: keeps studies with an E. coli count but no all-species count)
ana_ec <- dat_ok |> filter(!is.na(npos_ec), npos_ec <= N) |> mutate(nneg = N - npos_ec, prev = npos_ec / N)
rare_ec <- ana_ec |> count(subreg) |> filter(n < 4) |> pull(subreg) # similarly, combining rare subregions
ana_ec  <- ana_ec |> mutate(subreg_g = if_else(is.na(subreg) | subreg %in% rare_ec, "Other", subreg))

# counts for Fig 1 (study selection flow chart) "107 observations ... from 94 unique studies, covering 59804 participants in 56 countries"
say("[MS] observations = %d | unique studies = %d | participants = %s | countries = %d | subregions = %d",
    nrow(ana), n_distinct(ana$study_id), format(sum(ana$N), big.mark = ","),
    n_distinct(ana$Country), n_distinct(ana$subreg))

### 2. OUTCOME  (3GCR-E all species; 3GCR-E. coli for the sensitivity analysis) ####
# Primary outcome  : npos      / N   (all-species 3GCR-E)
# Sensitivity      : npos_ec   / N   (species-confirmed E. coli) -> section 12

### 3. WASH EXPOSURES - ladders and multiple imputation of safely-managed water ####
# --- 3a. Multiple imputation of safely-managed drinking water ------------------
# Missing because JMP has no water-quality data for some countries; treated as
# MAR given at-least-basic and piped coverage (both complete). Logit-scale linear
# regression with Smithson-Verkuilen squeeze; 20 posterior-predictive draws
# (Bayesian linear regression draw: sigma^2, then coefficients, then residual noise - the same
# algorithm as mice method "norm"). The outcome is NOT used.
n_i <- nrow(dat_w)
sv     <- function(p)  (p * (n_i - 1) + 0.5) / n_i      # squeeze (0,1) boundary values
sv_inv <- function(z)  (z * n_i - 0.5) / (n_i - 1)
imp_in <- tibble(l_sm  = qlogis(sv(dat_w$w_sm_obs / 100)),
                 l_bas = qlogis(sv(dat_w$w_bas   / 100)),
                 l_pip = qlogis(sv(dat_w$w_pip   / 100)))
# Rows used to fit the imputation model: observed SM water AND both predictors present.
# (Ref 140, Senegal rural, has observed SM water but no piped-water value -> not used for fitting.)
fit_rows  <- which(!is.na(imp_in$l_sm) & !is.na(imp_in$l_bas) & !is.na(imp_in$l_pip))
miss_rows <- which(is.na(imp_in$l_sm))                 # rows that need an imputed value
stopifnot(!anyNA(imp_in$l_bas[miss_rows]), !anyNA(imp_in$l_pip[miss_rows]))   # predictors needed for these only
imp_fit <- lm(l_sm ~ l_bas + l_pip, data = imp_in[fit_rows, ])
Xmis    <- cbind(1, imp_in$l_bas[miss_rows], imp_in$l_pip[miss_rows])
draw_l_sm <- function() {                              # one posterior-predictive draw of logit(SM water)
  s2 <- sum(resid(imp_fit)^2) / rchisq(1, df.residual(imp_fit))
  b  <- MASS::mvrnorm(1, coef(imp_fit), vcov(imp_fit) / sigma(imp_fit)^2 * s2)
  z  <- imp_in$l_sm
  z[miss_rows] <- as.numeric(Xmis %*% b) + rnorm(length(miss_rows), 0, sqrt(s2))
  z
}
l_sm_draws <- replicate(M_IMP, draw_l_sm(), simplify = FALSE)

derive_water <- function(d) {               # recompute exclusive basic rung + ladder score
  d |> mutate(w_basic  = pmax(w_bas - w_sm, 0),
              w_ladder = (w_unimp + 2*w_lim + 3*w_basic + 4*w_sm) /
                (w_surf + w_unimp + w_lim + w_basic + w_sm))
}
imps_ok <- map(seq_len(M_IMP), function(k) {
  z <- l_sm_draws[[k]]
  w <- pmin(100 * pmin(pmax(sv_inv(plogis(z)), 0), 1), dat_w$w_bas)   # SM cannot exceed at-least-basic
  dat_w |> mutate(w_sm = w) |> derive_water()
})
sub_imps <- function(d_set) map(imps_ok, ~ .x |> semi_join(d_set |> select(.row), by = ".row") |>
                                  left_join(d_set |> select(.row, subreg_g, nneg, prev), by = ".row"))
imps_ana <- sub_imps(ana)          # list of M_IMP datasets, primary set
imps_ec  <- sub_imps(ana_ec)       # list of M_IMP datasets, E. coli set
# (nneg/prev for ana_ec are recomputed from npos_ec inside glmm via the response argument)

# --- 3b. Sanitation ladder (direct JMP linkage; no imputation) ----------------
ana    <- ana    |> mutate(s_ladder = (s_unimp + 2*s_lim + 3*s_basic + 4*s_sm) / (s_ns + s_unimp + s_lim + s_basic + s_sm))
ana_ec <- ana_ec |> mutate(s_ladder = (s_unimp + 2*s_lim + 3*s_basic + 4*s_sm) / (s_ns + s_unimp + s_lim + s_basic + s_sm))
imps_ana <- map(imps_ana, ~ .x |> mutate(s_ladder = (s_unimp + 2*s_lim + 3*s_basic + 4*s_sm) / (s_ns + s_unimp + s_lim + s_basic + s_sm)))
imps_ec  <- map(imps_ec,  ~ .x |> mutate(s_ladder = (s_unimp + 2*s_lim + 3*s_basic + 4*s_sm) / (s_ns + s_unimp + s_lim + s_basic + s_sm)))

# WASH-analysis subsets: same rows and order as imps_ana / imps_ec (Taiwan, Refs 135-136, has no JMP series)
ana_w    <- ana    |> filter(wash_ok)
ana_ec_w <- ana_ec |> filter(wash_ok)

# validation: five mutually exclusive rungs sum to 100% (tolerance 1 percentage point for JMP rounding)
rung_check <- function(d, cols, label) {
  tot <- rowSums(d[, cols]); bad <- which(is.na(tot) | abs(tot - 100) >= 1)
  if (length(bad)) {
    say("WARNING: %s rungs do not sum to 100%% in %d row(s):", label, length(bad))
    print(d[bad, ] |> mutate(rung_sum = tot[bad]) |>
            select(ref, `First Author`, Country, setting, all_of(cols), rung_sum))
  } else say("%s rungs sum to 100%% in all rows", label)
  invisible(bad)
}
bad_w <- rung_check(imps_ana[[1]], c("w_surf","w_unimp","w_lim","w_basic","w_sm"), "Water")
bad_s <- rung_check(ana_w, c("s_ns","s_unimp","s_lim","s_basic","s_sm"), "Sanitation")

# For descriptives and plots only: values averaged over the imputed datasets
ana_mean <- ana_w |> select(-any_of(c("w_sm","w_basic","w_ladder"))) |>
  left_join(tibble(.row = imps_ana[[1]]$.row,
                   w_sm = rowMeans(sapply(imps_ana, `[[`, "w_sm"))), by = ".row") |>
  derive_water()

say("Missing water rungs (primary set, after imputation): %d", sum(is.na(ana_mean$w_ladder)))
say("Missing sanitation ladder (primary set): %d", sum(is.na(ana_w$s_ladder)))

### 4. ANTIBIOTIC USE AND DEMOGRAPHIC COVARIATES ####
# AMU (DDD/1000 inhabitants/day) was assembled in the database from Wang et al. S5, GLASS 2024,
# Browne et al., MAAP and Wirtz et al. (WLS interpolation for 16 countries); it is read as-is.
say("Missing AMU (primary set): %d", sum(is.na(ana$amu)))
# MS: "For 16 countries with no estimated antibiotic consumption values were imputed"
# (the database carries no flag column for imputed AMU, so this count is not recomputed here)

### 5. ANALYSIS PARAMETERS for GEE, beta-binomial GLMM, Rubin pooling ####
# ---- shared palette (Figures 2-6) -----------------------------------------
fam <- list(Africa   = c(light="#FEC44F", mid="#EC7014", dark="#8C2D04"),
            Asia     = c(light="#FC9272", mid="#DE2D26", dark="#67000D"),
            Europe   = c(light="#A1D99B", mid="#41AB5D", dark="#00441B"),
            Americas = c(light="#9ECAE1", mid="#4292C6", dark="#08306B"),
            Oceania  = c(light="#BCBDDC", mid="#807DBA", dark="#3F007D"))
cont_cols <- sapply(fam, `[[`, "dark")
cont_lty  <- c(Africa="solid", Asia="dashed", Europe="dotdash", Americas="twodash", Oceania="longdash")
subreg_levels <- sort(unique(na.omit(ana$subreg)))
subreg_cols <- unlist(lapply(names(fam), function(cc) {
  srs <- sort(subreg_levels[unname(cont_map[subreg_levels]) == cc])
  if (!length(srs)) return(NULL)
  setNames(colorRampPalette(c(fam[[cc]]["light"], fam[[cc]]["mid"]))(length(srs)), srs)
}))
water_order <- c("w_surf","w_unimp","w_lim","w_basic","w_sm")      # JMP ladder, bottom -> top
san_order   <- c("s_ns","s_unimp","s_lim","s_basic","s_sm")
water_rung_cols <- c(w_surf="#cc4c02", w_unimp="#fe9929", w_lim="#fee391", w_basic="#7bccc4", w_sm="#08589e")
san_rung_cols   <- c(s_ns="#bd0026",  s_unimp="#fe9929", s_lim="#fee391", s_basic="#7bccc4", s_sm="#08589e")
water_rung_lab  <- c(w_surf="Surface", w_unimp="Unimproved", w_lim="Limited", w_basic="Basic", w_sm="Safely managed")
san_rung_lab    <- c(s_ns="Open defecation", s_unimp="Unimproved", s_lim="Limited", s_basic="Basic", s_sm="Safely managed")
effect_cols <- c(`TRUE`="#01665E", `FALSE`="grey65")

# ---- GEE (population-averaged, clustered by country) ------------------------
GEE_WEIGHTED <- FALSE   # FALSE = every observation counts equally; TRUE = weighted by sample size
gee <- function(formula, d, weighted = GEE_WEIGHTED) {
  d <- d |> arrange(cid)
  if (!weighted) formula <- update(formula, prev ~ .)     # proportion response, unit weights
  suppressWarnings(                                       # "non-integer #successes" is harmless here
    geeglm(formula, id = cid, data = d, family = binomial(), corstr = "exchangeable"))
}
or_gee <- function(m, term) {
  b <- coef(m)[term]; se <- sqrt(diag(vcov(m)))[term]
  c(OR = unname(exp(b)), lo = unname(exp(b - 1.96*se)), hi = unname(exp(b + 1.96*se)))
}

# ---- beta-binomial GLMM (brms), fixed effects + random intercepts -----------
# Weakly informative N(0,2) prior on log-odds coefficients; brms defaults elsewhere.
# Bambi (used for the Python analysis) has different default priors.
glmm <- function(rhs, data, resp = "npos") {
  brm(bf(as.formula(paste0(resp, " | trials(N) ~ ", rhs))), data = data, family = beta_binomial(),
      prior = prior(normal(0, 2), class = b),
      chains = CHAINS, iter = ITER, cores = CHAINS, seed = SEED, refresh = 0, init = 0,
      control = list(adapt_delta = 0.95))
  }
# Fit once (compiles), then re-use the compiled model for the other imputed datasets.
# Returns the posterior draws of every fit so contrasts can be pooled with Rubin's rules.
fit_models <- function(rhs, data_list, resp = "npos") {
  if (!is.data.frame(data_list)) d1 <- data_list[[1]] else { d1 <- data_list; data_list <- list(d1) }
  f1 <- glmm(rhs, d1, resp)
  draws <- list(as.data.frame(as_draws_df(f1))); rh <- max(brms::rhat(f1), na.rm = TRUE)
  if (length(data_list) > 1) for (k in 2:length(data_list)) {
    fk <- update(f1, newdata = data_list[[k]], recompile = FALSE, refresh = 0, seed = SEED + k, init = 0)    
    draws[[k]] <- as.data.frame(as_draws_df(fk)); rh <- max(rh, brms::rhat(fk), na.rm = TRUE)
  }
  list(draws = draws, m = length(draws), rhat = rh, rhs = rhs)
}
# linear combination of coefficients, e.g. lc(c(b_w_basic = -1), sc = 10)
lc <- function(w, sc = 1) function(d) sc * as.numeric(as.matrix(d[, names(w), drop = FALSE]) %*% unname(w))
# Rubin's rules on the log-OR scale (reduces to the posterior summary when m = 1)
pool_ctr <- function(fit, fun, label) {
  q <- map_dbl(fit$draws, ~ mean(fun(.x))); u <- map_dbl(fit$draws, ~ var(fun(.x))); m <- fit$m
  if (m == 1) {
    x <- fun(fit$draws[[1]]); est <- mean(x); ci <- quantile(x, c(.025, .975)); fmi <- 0
  } else {
    B <- var(q); U <- mean(u); Tt <- U + (1 + 1/m) * B
    dfr <- (m - 1) * (1 + U / ((1 + 1/m) * B))^2
    est <- mean(q); ci <- est + c(-1, 1) * qt(.975, dfr) * sqrt(Tt); fmi <- (B + B/m) / Tt
  }
  tibble(term = label, OR = exp(est), lo = exp(ci[[1]]), hi = exp(ci[[2]]), fmi = fmi,
         text = fmt_or(exp(est), exp(ci[[1]]), exp(ci[[2]])))
}
# absolute prevalence difference for a given OR (draws stacked over imputations), baseline p0
abs_diff <- function(fit, fun, p0) {
  x <- unlist(map(fit$draws, fun)); dlt <- 100 * (plogis(qlogis(p0) + x) - p0)
  tibble(baseline_pct = 100 * p0, diff_pp = median(dlt),
         lo = quantile(dlt, .025), hi = quantile(dlt, .975))
}

BASE <- "amu + year_c + setting + (1|subreg_g) + (1|Country)"    # all primary models
# per-10-DDD scaling for AMU is applied when pooling (coefficient is per 1 DDD)

### 6. RESULTS - STUDY CHARACTERISTICS (Figure 2) ####
imp_avg <- function(f) mean(map_dbl(imps_ana, f))    # average of a statistic across imputed datasets

say("[MS] %d observations (%d countries), %s participants (median %.0f per observation, IQR %.0f-%.0f)",
    nrow(ana), n_distinct(ana$Country), format(sum(ana$N), big.mark = ","),
    median(ana$N), quantile(ana$N, .25), quantile(ana$N, .75))
say("[MS] collection years %d-%d (median %.0f)", min(ana$year), max(ana$year), median(ana$year))
sr <- ana |> count(subreg, sort = TRUE) |> mutate(pct = round(100 * n / nrow(ana), 1))
print(sr)   # MS: Eastern Africa 17 (15.9%), Eastern Asia 16 (15.0%), Western Europe & SE Asia 12 (11.2%) each
say("[MS] median 3GCR-E prevalence %.1f%% (IQR %.1f-%.1f); crude pooled %.1f%%",
    100*median(ana$prev), 100*quantile(ana$prev, .25), 100*quantile(ana$prev, .75), 100*sum(ana$npos)/sum(ana$N))

# subregion gradient (median prevalence by subregion)
print(ana |> group_by(subreg) |> summarise(n = n(), median_prev = round(100*median(prev), 1),
                                           pooled_prev = round(100*sum(npos)/sum(N), 1)) |> arrange(median_prev))

# pooled prevalence by period (abstract): TOTAL cases / TOTAL participants, sample-size weighted
period <- function(y) cut(y, c(-Inf, 2009, 2019, Inf), labels = c("<=2009", "2010-2019", "2020-2022"))
print(ana |> group_by(period = period(year)) |>
        summarise(n_obs = n(), cases = sum(npos), N = sum(N), pooled_pct = round(100*cases/N, 1),
                  mean_of_studies_pct = round(100*mean(prev), 1), median_pct = round(100*median(prev), 1)))

# ladder scores and rung coverage, averaged across imputations (water) / direct (sanitation)
wsc <- rowMeans(sapply(imps_ana, `[[`, "w_ladder"))      # water-ladder score averaged over imputations
say("[MS] median water-ladder score %.2f (IQR %.2f-%.2f); median sanitation-ladder score %.2f (IQR %.2f-%.2f)",
    median(wsc), quantile(wsc, .25), quantile(wsc, .75),
    median(ana_w$s_ladder, na.rm = TRUE), quantile(ana_w$s_ladder, .25, na.rm = TRUE), quantile(ana_w$s_ladder, .75, na.rm = TRUE))
wm <- sapply(water_order, function(r) imp_avg(~ mean(.x[[r]], na.rm = TRUE)))
sm <- colMeans(ana_w[, san_order], na.rm = TRUE)
say("[MS] water rungs (mean %% of population): %s", paste(sprintf("%s %.1f", water_rung_lab, wm[water_order]), collapse = "; "))
say("[MS] sanitation rungs: %s", paste(sprintf("%s %.1f", san_rung_lab, sm[san_order]), collapse = "; "))
say("[MS] water vs sanitation ladder-score correlation r = %.2f", cor(wsc, ana_w$s_ladder, use = "complete.obs"))

# ---- Figure 2 ---------------------------------------------------------------
keep <- ana |> count(subreg) |> filter(n >= 3) |> pull(subreg)
ord  <- ana |> filter(subreg %in% keep) |> group_by(subreg) |> summarise(m = median(prev)) |> arrange(m) |> pull(subreg)
ladder_cov <- function(d, cols)
  d |> filter(subreg %in% keep) |> group_by(subreg) |>
  summarise(across(all_of(cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  rowwise() |> mutate(across(all_of(cols), ~ 100 * .x / sum(c_across(all_of(cols))))) |> ungroup() |>
  pivot_longer(all_of(cols), names_to = "rung", values_to = "pct") |> mutate(rung = factor(rung, levels = cols))
fig2A <- ggplot(ladder_cov(ana_mean, water_order), aes(pct, factor(subreg, ord), fill = rung)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = water_rung_cols, labels = water_rung_lab, breaks = water_order, name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  labs(title = "A  Water ladder distribution", x = "Population (%)", y = NULL)
fig2B <- ggplot(ladder_cov(ana_mean, san_order), aes(pct, factor(subreg, ord), fill = rung)) +
  geom_col(position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = san_rung_cols, labels = san_rung_lab, breaks = san_order, name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  labs(title = "B  Sanitation ladder distribution", x = "Population (%)", y = NULL)
figure2 <- fig2A + fig2B
figure2
ggsave("figure2.jpg", plot = figure2, width = 12, height = 4, dpi = 300)

### 7. RESULTS - MULTIPLE IMPUTATION OF SAFELY-MANAGED WATER COVERAGE ####
n_miss <- sum(ana_w$w_sm_imp)
say("[MS] safely-managed water missing for %d of %d observations (%.0f%%)", n_miss, nrow(ana_w), 100*n_miss/nrow(ana_w))

obs_idx <- fit_rows
fit_all <- imp_fit
loo <- map_dfr(obs_idx, function(i) {
  f <- lm(l_sm ~ l_bas + l_pip, data = imp_in[setdiff(obs_idx, i), ])
  tibble(i = i, obs = 100*sv_inv(plogis(imp_in$l_sm[i])),
         pred = 100*sv_inv(plogis(predict(f, imp_in[i, ]))))
})
say("[MS] imputation model: R2 = %.2f (logit scale); LOOCV r = %.2f, MAE = %.1f percentage points (n = %d observed)",
    summary(fit_all)$r.squared, cor(loo$obs, loo$pred), mean(abs(loo$obs - loo$pred)), length(obs_idx))
say("     R2 on the percentage scale (obs vs LOOCV prediction) = %.2f", cor(loo$obs, loo$pred)^2)

imp_rows <- dat_w |> mutate(w_sm_mean = rowMeans(sapply(imps_ok, `[[`, "w_sm"))) |> filter(w_sm_imp)
say("[MS] imputed values range %.1f%% (%s) to %.1f%% (%s)",
    min(imp_rows$w_sm_mean), imp_rows$Country[which.min(imp_rows$w_sm_mean)],
    max(imp_rows$w_sm_mean), imp_rows$Country[which.max(imp_rows$w_sm_mean)])

### 8. RESULTS - TIME TREND (Table 1, Figure 3) ####
mt <- gee(cbind(npos, nneg) ~ year_c, ana)
mq <- gee(cbind(npos, nneg) ~ year_c + I(year_c^2), ana)
o  <- or_gee(mt, "year_c")
say("[MS] overall OR per year %s; quadratic term p = %.2f",
    fmt_or(o["OR"], o["lo"], o["hi"]), summary(mq)$coefficients["I(year_c^2)", "Pr(>|W|)"])
say("     = %.1f%% per year on the odds scale",
    100*(o["OR"] - 1))
for (cc in c("Africa", "Asia", "Europe")) {
  oc <- or_gee(gee(cbind(npos, nneg) ~ year_c, filter(ana, continent == cc)), "year_c")
  say("     %-7s OR per year %s", cc, fmt_or(oc["OR"], oc["lo"], oc["hi"]))
}

mq <- gee(cbind(npos, nneg) ~ year_c + I(year_c^2), ana)
print(summary(mq)$coefficients)                       # sign, estimate, SE, p
print(c(n_obs = nrow(ana), n_countries = n_distinct(ana$cid)))
o <- or_gee(gee(cbind(npos, nneg) ~ year_c, ana), "year_c"); print(o)   # should be ~1.14 (1.07-1.21)
b <- coef(mq); say("turning point at %.0f", 2015 - b["year_c"] / (2 * b["I(year_c^2)"]))
gr <- tibble(year_c = -13:7); print(gr |> mutate(year = year_c + 2015,
                                                 pct = 100 * plogis(b[1] + b[2]*year_c + b[3]*year_c^2)))

# ---- Table 1: standardised to 2015 (interpolation) and projected to 2025 -------
ana_t <- ana |> filter(subreg %in% keep, !is.na(continent)) |>
  mutate(subreg = factor(subreg), continent = factor(continent)) |> arrange(cid)
m2015 <- gee(cbind(npos, nneg) ~ subreg + year_c, ana_t)
mproj <- gee(cbind(npos, nneg) ~ subreg + year_c:continent, ana_t)
pci <- function(fit, newdata) {
  X <- model.matrix(delete.response(terms(fit)), newdata,
                    xlev = list(subreg = levels(ana_t$subreg), continent = levels(ana_t$continent)))
  ok <- !is.na(coef(fit)); b <- coef(fit)[ok]; X <- X[, ok, drop = FALSE]
  V <- vcov(fit); if (nrow(V) == length(ok)) V <- V[ok, ok, drop = FALSE]
  eta <- as.numeric(X %*% b); se <- sqrt(diag(X %*% V %*% t(X)))
  round(100 * plogis(eta + c(0, -1.96, 1.96) * se), 1)
}
sub_cont <- ana_t |> distinct(subreg, continent); sub_cont <- setNames(as.character(sub_cont$continent), as.character(sub_cont$subreg))
table1 <- map_dfr(levels(ana_t$subreg), function(s) {
  a <- pci(m2015, tibble(subreg = factor(s, levels(ana_t$subreg)), year_c = 0))
  p <- pci(mproj, tibble(subreg = factor(s, levels(ana_t$subreg)),
                         continent = factor(sub_cont[[s]], levels(ana_t$continent)), year_c = 10))
  tibble(subregion = s, n = sum(ana_t$subreg == s),
         p2015 = a[1], lo2015 = a[2], hi2015 = a[3], p2025 = p[1], lo2025 = p[2], hi2025 = p[3])
}) |> arrange(desc(p2015))
print(table1)   # MS: 2015 from 4.3% (3.1-5.8) Southern Europe to 57.9% (38.0-75.5) South-eastern Asia
# 2025 values assume the log-linear trend simply continues (trend-based extrapolation, not a forecast)
write_csv(table1, "table1_reanalysis.csv")

# ---- WASH coverage over time (JMP, matched to the study country/stratum) ------
jmp_nm <- function(x, nm) {
  fix <- c(Vietnam="Viet Nam", Laos="Lao People's Democratic Republic", Iran="Iran (Islamic Republic of)",
           Tanzania="United Republic of Tanzania", `The Gambia`="Gambia", Venezuela="Venezuela (Bolivarian Republic of)")
  if (x %in% nm) x else if (x %in% names(fix)) fix[[x]] else { h <- nm[startsWith(nm, x)]; if (length(h)) h[1] else NA_character_ }
}
jmp_get <- function(tab, nm, y, col) { v <- tab[[col]][tab$name == nm & tab$year == y]; if (length(v)) v[1] else NA_real_ }
wash_trend <- function(tab, stem, yr) {
  # unweighted mean over the UNIQUE study countries (one value per country), national total stratum
  nm_all <- na.omit(unique(tab$name))
  nms <- unique(na.omit(map_chr(ana$Country, jmp_nm, nm = nm_all)))
  mean(map_dbl(nms, ~ jmp_get(tab, .x, yr, paste0(stem, "_t"))), na.rm = TRUE)
}
say("[MS] mean WASH coverage at the study country/stratum, 2000 -> 2022 (JMP; Taiwan has no JMP series):")
say("     safely managed water   %.1f -> %.1f", wash_trend(jmp_w, "wat_sm", 2000),    wash_trend(jmp_w, "wat_sm", 2022))
say("     at-least-basic water   %.1f -> %.1f", wash_trend(jmp_w, "wat_basal", 2000), wash_trend(jmp_w, "wat_basal", 2022))
say("     safely managed sanit.  %.1f -> %.1f", wash_trend(jmp_s, "san_sm", 2000),    wash_trend(jmp_s, "san_sm", 2022))
say("     at-least-basic sanit.  %.1f -> %.1f", wash_trend(jmp_s, "san_basal", 2000), wash_trend(jmp_s, "san_basal", 2022))

# ---- Figure 3 -----------------------------------------------------------------
int_years <- scale_x_continuous(breaks = breaks_width(2), labels = label_number(accuracy = 1, big.mark = ""))
cont_pred <- map_dfr(names(fam), function(cc) {
  d <- ana |> filter(continent == cc)
  if (n_distinct(d$Country) < 2) return(NULL)
  m <- tryCatch(gee(cbind(npos, nneg) ~ year_c, d), error = function(e) NULL)
  if (is.null(m) || anyNA(coef(m))) return(NULL)            # Americas / Oceania: too few observations
  gr <- tibble(year_c = seq(min(d$year_c), max(d$year_c), length.out = 50))
  gr |> mutate(continent = cc, year = year_c + 2015, prev = 100 * plogis(predict(m, newdata = gr)))
})
fig3A <- ggplot(ana, aes(year, 100*prev)) +
  geom_point(aes(colour = subreg, size = sqrt(N)), alpha = 0.75) +
  scale_colour_manual(values = subreg_cols, limits = names(subreg_cols), drop = FALSE, name = "Subregion", na.translate = FALSE) +
  scale_size_continuous(guide = "none") + new_scale_colour() +
  geom_line(data = cont_pred, aes(year, prev, colour = continent, linetype = continent), linewidth = 1, inherit.aes = FALSE) +
  scale_colour_manual(values = cont_cols, name = "Continent trend") +
  scale_linetype_manual(values = cont_lty, name = "Continent trend") + int_years +
  labs(title = "A  3GCR-E prevalence by year", x = "Year of sample collection", y = "3GCR-E prevalence (%)")
fig3B <- ggplot(filter(ana, subreg %in% keep), aes(100*prev, factor(subreg, ord), fill = subreg)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) + scale_fill_manual(values = subreg_cols) +
  labs(title = "B  3GCR-E prevalence by UN M49 subregion", x = "3GCR-E prevalence (%)", y = NULL)
fig3C <- ggplot(filter(ana, !is.na(npos_cre)), aes(year, 100*npos_cre/N)) +
  geom_point(aes(colour = subreg, size = sqrt(N)), alpha = 0.75) +
  scale_colour_manual(values = subreg_cols, limits = names(subreg_cols), drop = FALSE, name = "Subregion", na.translate = FALSE) +
  scale_size_continuous(guide = "none") + new_scale_colour() +
  geom_smooth(aes(colour = continent, linetype = continent), method = "lm", se = FALSE, linewidth = 1) +
  scale_colour_manual(values = cont_cols, name = "Continent trend") +
  scale_linetype_manual(values = cont_lty, name = "Continent trend") + int_years +
  coord_cartesian(ylim = c(0, 15)) +
  labs(title = "C  Carbapenem-resistant Enterobacterales prevalence", x = "Year of sample collection", y = "CRE prevalence (%)")
figure3 <- fig3A + fig3B + (fig3C + theme(legend.position = "none")) +  plot_layout(nrow = 1, guides = "collect")
figure3
ggsave("figure3.jpg", plot = figure3, dpi = 300, width = 16, height = 7)

### 9. RESULTS - ASSOCIATION: COLONISATION, WASH AND ANTIBIOTIC USE (Bayesian beta-binomial GLMMs & Fig 4) ####
# Water models are fitted in every imputed dataset and pooled with Rubin's rules;
# sanitation has no imputation (single fit). Adjusted for year, setting, AMU.
fit_w     <- fit_models(paste("w_ladder +", BASE), imps_ana)                  # water ladder score
fit_s     <- fit_models(paste("s_ladder +", BASE), ana_w)                       # sanitation ladder score
fit_mutual<- fit_models(paste("w_ladder + s_ladder +", BASE), imps_ana)       # collinearity check

res_w <- pool_ctr(fit_w, lc(c(b_w_ladder = 1)), "Water ladder, per rung")
res_s <- pool_ctr(fit_s, lc(c(b_s_ladder = 1)), "Sanitation ladder, per rung")
res_amu <- pool_ctr(fit_w, lc(c(b_amu = 1), sc = 10), "AMU, per +10 DDD/1000/day")
say("[MS] one step up the water ladder        OR %s", res_w$text)
say("[MS] one step up the sanitation ladder   OR %s", res_s$text)
say("[MS] antibiotic use per 10 DDD           OR %s", res_amu$text)
say("[MS] fraction of missing information (water ladder) = %.1f%%", 100*res_w$fmi)
say("     max R-hat: water %.3f | sanitation %.3f", fit_w$rhat, fit_s$rhat)
say("[MS] water vs sanitation ladder r = %.2f (collinearity)", cor(wsc, ana_w$s_ladder, use = "complete.obs"))
print(bind_rows(pool_ctr(fit_mutual, lc(c(b_w_ladder = 1)), "Water ladder, mutually adjusted for sanitation"),
                pool_ctr(fit_mutual, lc(c(b_s_ladder = 1)), "Sanitation ladder, mutually adjusted for water")) |>
        select(term, text))

# complete-case: only observations with directly observed safely-managed water
ana_cc <- ana_w |> filter(!w_sm_imp) |>
  mutate(w_sm = w_sm_obs, w_basic = pmax(w_bas - w_sm, 0)) |> derive_water()
fit_cc <- fit_models(paste("w_ladder +", BASE), ana_cc)
say("[MS] complete case (n = %d): water ladder OR %s",
    nrow(ana_cc), pool_ctr(fit_cc, lc(c(b_w_ladder = 1)), "cc")$text)

# ---- Figure 4 ------------------------------------------------------------------
ana_mean <- ana_mean |> mutate(s_ladder = ana_w$s_ladder[match(.row, ana_w$.row)])
ladder_panel <- function(score, xlab, title) {
  d  <- ana_mean |> filter(!is.na(.data[[score]]))
  m  <- gee(as.formula(paste("cbind(npos, nneg) ~", score)), d)
  gr <- tibble(x = seq(min(d[[score]]), max(d[[score]]), length.out = 100))
  X  <- cbind(1, gr$x); eta <- as.numeric(X %*% coef(m)); se <- sqrt(diag(X %*% vcov(m) %*% t(X)))
  gr <- gr |> mutate(fit = 100*plogis(eta), lo = 100*plogis(eta - 1.96*se), hi = 100*plogis(eta + 1.96*se))
  ggplot(d, aes(.data[[score]], 100*prev)) +
    geom_ribbon(data = gr, aes(x, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey30", alpha = 0.12) +
    geom_point(aes(colour = subreg, size = sqrt(N)), alpha = 0.75) +
    geom_line(data = gr, aes(x, fit), inherit.aes = FALSE, colour = "grey15", linewidth = 1) +
    geom_vline(xintercept = 1:4, colour = "grey80", linetype = "dotted", linewidth = 0.3) +
    scale_colour_manual(values = subreg_cols, limits = names(subreg_cols), name = "UN M49 subregion") +
    scale_size_continuous(guide = "none") + xlim(0, 4) +
    labs(title = title, x = xlab, y = "3GCR-E prevalence (%)")
}
fig4A <- ladder_panel("w_ladder", "Water-ladder score (0 = surface water ... 4 = safely managed)", "A  Colonisation vs water-ladder score")
fig4B <- ladder_panel("s_ladder", "Sanitation-ladder score (0 = open defecation ... 4 = safely managed)", "B  Colonisation vs sanitation-ladder score")
fig4C <- ggplot(filter(ana_mean, !is.na(amu)), aes(amu, 100*prev)) +
  geom_point(aes(colour = subreg, size = sqrt(N)), alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, colour = "black", linetype = 2, linewidth = 0.8) +
  scale_colour_manual(values = subreg_cols, limits = names(subreg_cols), name = "UN M49 subregion") +
  scale_size_continuous(guide = "none") +
  labs(title = "C  Colonisation vs antibiotic use", x = "Antibiotic use (DDD/1000/day)", y = "3GCR-E prevalence (%)")
fig4D <- ggplot(filter(ana_mean, !is.na(w_ladder), !is.na(s_ladder)), aes(w_ladder, s_ladder)) +
  geom_point(aes(colour = 100*prev, size = sqrt(N)), alpha = 0.85) +
  scale_colour_viridis_c(option = "magma", direction = -1, name = "Prevalence (%)") +
  scale_size_continuous(guide = "none") + xlim(0, 4) + ylim(0, 4) +
  labs(title = "D  Water vs sanitation ladder jointly", x = "Water-ladder score (0-4)", y = "Sanitation-ladder score (0-4)")
figure4 <- (fig4A + fig4B) / (fig4C + fig4D) + plot_layout(guides = "collect") & theme(legend.position = "right")
figure4
ggsave("figure4.jpg", plot = figure4, width = 12, height = 9, dpi = 300)

### 10. RESULTS - GAINS FROM SPECIFIC STEPS UP THE LADDER (Compositional (rung-level) models & Figure 5) ####
# Rung percentages enter the model with the TOP rung (safely managed) as reference, so the
# effect of moving 10% of the population from rung k to k+1 is 10 * (b[k+1] - b[k]), with b[top] = 0.
fit_wc <- fit_models(paste("w_surf + w_unimp + w_lim + w_basic +", BASE), imps_ana)
fit_sc <- fit_models(paste("s_ns + s_unimp + s_lim + s_basic +", BASE), ana_w)
fit_wsm <- fit_models(paste("w_sm +", BASE), imps_ana)      # continuous safely-managed water, per +10 pp
fit_ssm <- fit_models(paste("s_sm +", BASE), ana_w)          # continuous safely-managed sanitation, per +10 pp

w_steps <- list(
  "Surface → unimproved"       = lc(c(b_w_unimp = 1, b_w_surf = -1), 10),
  "Unimproved → limited"       = lc(c(b_w_lim = 1, b_w_unimp = -1), 10),
  "Limited → basic"            = lc(c(b_w_basic = 1, b_w_lim = -1), 10),
  "Basic → safely managed"     = lc(c(b_w_basic = -1), 10))
s_steps <- list(
  "Open defecation → unimproved" = lc(c(b_s_unimp = 1, b_s_ns = -1), 10),
  "Unimproved → limited"         = lc(c(b_s_lim = 1, b_s_unimp = -1), 10),
  "Limited → basic"              = lc(c(b_s_basic = 1, b_s_lim = -1), 10),
  "Basic → safely managed"       = lc(c(b_s_basic = -1), 10))
res_wsteps <- imap_dfr(w_steps, ~ pool_ctr(fit_wc, .x, .y))
res_ssteps <- imap_dfr(s_steps, ~ pool_ctr(fit_sc, .x, .y))
res_wsm <- pool_ctr(fit_wsm, lc(c(b_w_sm = 1), 10), "Safely-managed water, per +10 pp")
res_ssm <- pool_ctr(fit_ssm, lc(c(b_s_sm = 1), 10), "Safely-managed sanitation, per +10 pp")
cat("\nWATER transitions (OR for moving 10% of the population up one rung)\n");      print(res_wsteps |> select(term, text))
cat("\nSANITATION transitions\n");                                                    print(res_ssteps |> select(term, text))

# as a check, any rung (not just basic) to safely managed (not withheld in the text to avoid confusion)
say("[MS] safely-managed sanitation per +10pp OR %s", res_ssm$text)
say("[MS] safely-managed water per +10pp      OR %s", res_wsm$text)

# absolute prevalence differences for the basic -> safely managed step at baseline 20% and 40%
cat("\nAbsolute prevalence change (percentage points) for 10% moving basic -> safely managed:\n")
print(bind_rows(
  map_dfr(c(.2, .4), ~ abs_diff(fit_wc, w_steps[[4]], .x) |> mutate(ladder = "water")),
  map_dfr(c(.2, .4), ~ abs_diff(fit_sc, s_steps[[4]], .x) |> mutate(ladder = "sanitation"))))

# ---- Figure 5 (top = highest rung; bottom row = continuous safely-managed coverage) ----
fmt_ci <- function(or, lo, hi)
  gsub(".", "\u00B7", sprintf("%.2f (%.2f\u2013%.2f)", or, lo, hi), fixed = TRUE)

mk_fig5 <- function(overall, steps, title) {
  lev <- c(steps$term, "Whole ladder, per step")                           # bottom -> top
  d <- bind_rows(overall |> mutate(term = "Whole ladder, per step"), steps) |>
    mutate(sig   = lo > 1 | hi < 1,
           kind  = if_else(term %in% steps$term, "Transition (10% of population)", "Whole ladder"),
           term  = factor(term, levels = lev),
           label = fmt_ci(OR, lo, hi))
  ggplot(d, aes(OR, term, colour = sig, shape = kind)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_hline(yintercept = length(lev) - 0.5, colour = "grey85") +        # separates whole ladder from transitions
    geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = 0.9, fatten = 3) +
    geom_text(aes(label = label), nudge_y = -0.32, size = 2.7,             # small label, centred under the estimate
              colour = "grey25", show.legend = FALSE) +
    scale_x_log10() +
    scale_y_discrete(expand = expansion(add = c(0.8, 0.5))) +               # room for the label under the bottom row
    scale_colour_manual(values = effect_cols, guide = "none") +
    scale_shape_manual(values = c(`Whole ladder` = 15, `Transition (10% of population)` = 16), name = NULL) +
    labs(title = title, x = "Adjusted OR for 3GCR-E colonisation (log scale)", y = NULL)
}
fig5A <- mk_fig5(res_w, res_wsteps, "A  Water ladder")
fig5B <- mk_fig5(res_s, res_ssteps, "B  Sanitation ladder")
figure5 <- fig5A + fig5B + plot_layout(guides = "collect") & theme(legend.position = "none")
figure5
ggsave("figure5.jpg", plot = figure5, width = 12, height = 3, dpi = 300)

### 11. RESULTS - RURAL VS URBAN SETTING (Figure 6) ####

ur <- ana_w |> filter(setting != "Unspecified")
say("[MS] %d observations specified urban (n=%d) or rural (n=%d); %d unspecified (national estimates)",
    nrow(ur), sum(ur$setting == "Urban"), sum(ur$setting == "Rural"), sum(ana_w$setting == "Unspecified"))
pooled_by <- ana_w |> group_by(setting) |> summarise(pooled_pct = 100*sum(npos)/sum(N))
say("[MS] crude pooled prevalence rural %.1f%% vs urban %.1f%%",
    pooled_by$pooled_pct[pooled_by$setting == "Rural"], pooled_by$pooled_pct[pooled_by$setting == "Urban"])
wsm_mean <- rowMeans(sapply(imps_ana, `[[`, "w_sm"))
say("[MS] median safely-managed water: rural %.1f%% vs urban %.1f%%",
    median(wsm_mean[ana_w$setting == "Rural"]), median(wsm_mean[ana_w$setting == "Urban"]))
imps_ur <- map(imps_ana, ~ .x |> filter(setting != "Unspecified") |> mutate(setting = droplevels(setting)))
fit_set  <- fit_models("setting + year_c + (1|subreg_g)",        imps_ur)
fit_setw <- fit_models("setting + w_sm + year_c + (1|subreg_g)", imps_ur)
say("[MS] rural vs urban (year-adjusted)            OR %s", pool_ctr(fit_set,  lc(c(b_settingRural = 1)), "r")$text)
say("[MS] rural vs urban (+ safely-managed water)   OR %s", pool_ctr(fit_setw, lc(c(b_settingRural = 1)), "r")$text)

fig6A <- ggplot(filter(ana_mean, setting != "Unspecified"), aes(setting, w_sm, fill = setting)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) + geom_jitter(width = .12, alpha = .5, size = 1) +
  scale_fill_manual(values = c(Urban = "#9ECAE1", Rural = "#FDAE6B")) +
  labs(title = "A  Safely-managed water", x = NULL, y = "Population with safely-managed water (%)")
fig6B <- ggplot(filter(ana_mean, setting != "Unspecified"), aes(setting, 100*prev, fill = setting)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) + geom_jitter(width = .12, alpha = .5, size = 1) +
  scale_fill_manual(values = c(Urban = "#9ECAE1", Rural = "#FDAE6B")) +
  labs(title = "B  3GCR-E prevalence", x = NULL, y = "3GCR-E prevalence (%)")
figure6 <- fig6A + fig6B
figure6
ggsave("figure6.jpg", plot = figure6, width = 8, height = 3, dpi = 300)

### 12. RESULTS - SENSITIVITY ANALYSIS: 3GCR-E. COLI ONLY ####
n_ec_in_ana <- sum(!is.na(ana$npos_ec))
ec_in_ana <- ana |> filter(!is.na(npos_ec))
say("[MS] %d of %d observations reported species-confirmed 3GCR-E. coli counts", n_ec_in_ana, nrow(ana))
say("[MS] pooled 3GCR-E. coli %.1f%% (median %.1f%%) vs all-species %.1f%% (median %.1f%%)",
    100*sum(ec_in_ana$npos_ec)/sum(ec_in_ana$N), 100*median(ec_in_ana$npos_ec/ec_in_ana$N),
    100*sum(ana$npos)/sum(ana$N), 100*median(ana$prev))
fit_ecw <- fit_models(paste("w_ladder +", BASE), imps_ec, resp = "npos_ec")
fit_ecs <- fit_models(paste("s_ladder +", BASE), ana_ec_w,  resp = "npos_ec")
say("[MS] E. coli, one step up water ladder       OR %s (manuscript 0.80, 0.47-1.41)", pool_ctr(fit_ecw, lc(c(b_w_ladder = 1)), "w")$text)
say("[MS] E. coli, one step up sanitation ladder  OR %s (manuscript 0.82, 0.56-1.21)", pool_ctr(fit_ecs, lc(c(b_s_ladder = 1)), "s")$text)

### 13. COLLINEARITY OF JMP LADDER SCORES WITH GDP PER CAPITA ####
gd <- ana_mean |> filter(!is.na(gdp))
say("GDP available for %d of %d observations (World Bank; IMF for Taiwan)", nrow(gd), nrow(ana_w))
say("Pearson r (log GDP vs water-ladder score)       = %.2f", cor(log(gd$gdp), gd$w_ladder, use = "complete.obs"))
say("Pearson r (log GDP vs sanitation-ladder score)  = %.2f", cor(log(gd$gdp), gd$s_ladder, use = "complete.obs"))
say("Pearson r (water vs sanitation ladder score)    = %.2f", cor(gd$w_ladder, gd$s_ladder, use = "complete.obs"))
# GDP is NOT added to the primary models (the analysis is kept to WASH vs colonisation).

### 15. OPTIONAL: LMIC-ONLY SENSITIVITY (supplementary; not in the current text) ####
if (RUN_LMIC) {
  HIC <- c("Australia","Belgium","Cyprus","Denmark","France","Germany","Hungary","Japan","Netherlands",
           "Norway","Portugal","Qatar","Republic of Korea","Saudi Arabia","Singapore","Spain","Sweden",
           "Switzerland","United States of America")
  imps_lmic <- map(imps_ana, ~ filter(.x, !(Country %in% HIC)))
  say("LMIC-only set: %d observations, %d countries", nrow(imps_lmic[[1]]), n_distinct(imps_lmic[[1]]$Country))
  f_lw <- fit_models(paste("w_ladder +", BASE), imps_lmic)
  f_ls <- fit_models(paste("s_ladder +", BASE), filter(ana_w, !(Country %in% HIC)))
  print(bind_rows(pool_ctr(f_lw, lc(c(b_w_ladder = 1)), "LMIC water ladder"),
                  pool_ctr(f_ls, lc(c(b_s_ladder = 1)), "LMIC sanitation ladder")) |> select(term, text))
}
