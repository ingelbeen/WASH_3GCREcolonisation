## figure6_final.R
## ---------------------------------------------------------------------------
## Figure 6 — final analyses only.
##   Panel A: study-level 3GCR-E colonisation (each study's own observed/crude
##            prevalence, from its own positives/N) vs that study's own
##            country's nearest-year estimate of E. coli 3GC-R bloodstream (or
##            clinical, where noted) resistance from WHO GLASS, supplemented
##            for GLASS non-reporters by EARS-Net (Spain), CDC (USA), CARSS
##            (China). One point per included study; matched year is nearest
##            to that study's own sample-collection year (cols 74-77 of the
##            v15 workbook).
##   Panel B: GBD-region 3GCR-E colonisation (pooled directly from included
##            study counts, standardised to 2019) vs the GBD 2019 E. coli /
##            3GC attributable death rate for the same GBD region.
##
## Superseded material removed from this version: the earlier GBD SUPER-region
## (n=7) trend (odds-change-per-decade) and size (2021 levels) comparisons are
## no longer part of the analysis and are not reproduced here.
## ---------------------------------------------------------------------------

library(readxl); library(dplyr); library(tidyr); library(stringr)
library(purrr);  library(ggplot2); library(patchwork)

V15   <- "ESBL_CRE_FecalColonization_HealthyIndividuals_GlobalEvidence_v15_imputed_missing_safely_managed.xlsx"
GLASS <- "GLASS_percentages_3GCRE.xlsx"
GBD_RATES_2019 <- "figure6_gbd_v15_data.csv"    # region, prev2019, rate, ... (GBD 2019 attributable death rate by GBD region; see prior extraction from IHME_AMR_BURDEN_2019_RATE)

m49_colors <- c(
  'Eastern Africa'='#FEC44F','Middle Africa'='#F9AF40','Northern Africa'='#F59A31',
  'Southern Africa'='#F08523','Western Africa'='#EC7014',
  'Eastern Asia'='#FC9272','South-eastern Asia'='#F27059','Southern Asia'='#E84F3F','Western Asia'='#DE2D26',
  'Eastern Europe'='#A1D99B','Northern Europe'='#81CA86','Southern Europe'='#61BA72','Western Europe'='#41AB5D',
  'Central America'='#9ECAE1','Northern America'='#70AED3','South America'='#4292C6',
  'Australia and New Zealand'='#BCBDDC')

mid_year <- function(s) {
  yrs <- regmatches(as.character(s), gregexpr("(19|20)[0-9]{2}", as.character(s)))
  vapply(yrs, function(y) if (length(y)) mean(as.numeric(y)) else NA_real_, numeric(1))
}

## ---- shared: clean included studies, reconstruct pos/N, resolve a year ----
raw <- read_excel(V15, sheet = "Primary Studies") |> rename_with(str_squish)
pick <- function(p) names(raw)[startsWith(names(raw), p)][1]

inc <- raw |>
  filter(.data[[pick("Include in Analysis")]] == "Include") |>
  transmute(
    country_raw = .data[[pick("Country")]],
    N   = suppressWarnings(as.numeric(.data[[pick("N Stool")]])),
    pos = suppressWarnings(as.numeric(.data[[pick("N Positive ESBL/3GCR")]])),
    prev_pub = suppressWarnings(as.numeric(.data[[pick("Prevalence ESBL")]])),
    yr_sample = .data[[pick("Year(s) Sample")]],
    yr_pub = .data[[pick("Year Published")]],
    m49 = .data[[pick("UN M49")]]) |>
  mutate(
    pos = if_else(is.na(pos) & !is.na(prev_pub) & !is.na(N), round(prev_pub/100*N), pos),
    year = mid_year(yr_sample),
    year = if_else(is.na(year) & !is.na(yr_pub), yr_pub - 1, year),
    country = recode(country_raw,
      "China (Central South)"="China", "Ethiopia (Somali region)"="Ethiopia",
      "Tanzania (Mwanza)"="Tanzania", "United States of America"="USA")) |>
  filter(!is.na(N), !is.na(pos), !is.na(year), N > 0, pos >= 0, pos <= N)

## GBD region for panel B (Caribbean has no matching 2019 death rate and is
## dropped downstream; kept here for completeness of the crosswalk)
gbd_region_of <- function(country) {
  REG <- list(
    "Australasia"=c("Australia"), "Central Europe"=c("Hungary"),
    "Central Latin America"=c("Mexico","Guatemala","Venezuela"), "Caribbean"=c("Guyana"),
    "Central Sub-Saharan Africa"=c("Angola","Central African Republic","Cameroon","Chad"),
    "East Asia"=c("China","China (Central South)"),
    "Eastern Sub-Saharan Africa"=c("Ethiopia","Ethiopia (Somali region)","Kenya","Tanzania",
      "Tanzania (Mwanza)","Uganda","Rwanda","Malawi","Mozambique","Madagascar"),
    "High-income Asia Pacific"=c("Japan","Republic of Korea","Singapore","Taiwan"),
    "High-income North America"=c("USA","United States of America"),
    "North Africa and Middle East"=c("Iran","Saudi Arabia","Tunisia","Egypt","Lebanon","Morocco"),
    "South Asia"=c("India","Bangladesh","Nepal"),
    "Southeast Asia"=c("Thailand","Vietnam","Cambodia","Laos","Indonesia"),
    "Southern Latin America"=c("Chile"), "Southern Sub-Saharan Africa"=c("Botswana"),
    "Western Europe"=c("Belgium","Cyprus","Denmark","France","Germany","Netherlands",
      "Norway","Portugal","Spain","Sweden","Switzerland"),
    "Western Sub-Saharan Africa"=c("Ghana","Senegal","Burkina Faso","The Gambia","Nigeria"))
  hit <- names(REG)[map_lgl(REG, ~ country %in% .x)]
  if (length(hit)) hit[1] else NA_character_
}
inc <- inc |> mutate(gbd_region = map_chr(country_raw, gbd_region_of))

## ---- shared pooling machinery: WLS logit meta-regression, DL heterogeneity ----
wls <- function(X, y, w) { cov <- solve(t(X) %*% (X*w)); list(b = cov %*% (t(X)%*%(y*w)), cov = cov) }
pool_logit <- function(d, group_var, ref_year) {
  d <- d |> mutate(lp = qlogis((pos+.5)/(N+1)), v = 1/(pos+.5) + 1/(N-pos+.5), yc = year - ref_year)
  grp <- sort(unique(d[[group_var]]))
  X <- cbind(sapply(grp, \(g) as.numeric(d[[group_var]] == g)), yc = d$yc)
  y <- d$lp; v <- d$v
  b0 <- wls(X, y, 1/v)$b; res <- y - X %*% b0; w <- 1/v
  tau2 <- max(0, (sum(res^2/v) - (length(y)-ncol(X))) / (sum(w)-sum(w^2)/sum(w)))
  fit <- wls(X, y, 1/(v+tau2)); beta <- as.numeric(fit$b); se <- sqrt(diag(fit$cov))
  k <- length(grp)
  tibble(group = grp, prev = plogis(beta[1:k]),
         lo = plogis(beta[1:k]-1.96*se[1:k]), hi = plogis(beta[1:k]+1.96*se[1:k]),
         n_studies = as.integer(table(factor(d[[group_var]], grp))),
         sumN = map_dbl(grp, \(g) sum(d$N[d[[group_var]]==g])))
}

## =============================================================================
## PANEL A — study-level colonisation (each study's own crude prevalence) vs
##           that study's own nearest-year GLASS / EARS-Net (Spain) / CDC
##           (USA) / CARSS (China) estimate for its country
## =============================================================================
## The v15 workbook already carries this per-study match (cols 74-77, written
## in an earlier pass): "AMR match (%)", "AMR match year", "AMR match source".
## Panel A simply reads those columns back alongside each study's own crude
## prevalence — no re-pooling or re-matching is needed here.
amr_col <- names(raw)[str_starts(names(raw), "AMR match \\(%\\)")]

colA <- raw |>
  filter(.data[[pick("Include in Analysis")]] == "Include") |>
  transmute(
    ref = .data[[pick("Ref #")]], country = .data[[pick("Country")]],
    m49 = .data[[pick("UN M49")]],
    N   = suppressWarnings(as.numeric(.data[[pick("N Stool")]])),
    pos = suppressWarnings(as.numeric(.data[[pick("N Positive ESBL/3GCR")]])),
    prev_pub = suppressWarnings(as.numeric(.data[[pick("Prevalence ESBL")]])),
    match_pct = suppressWarnings(as.numeric(.data[[amr_col]])),
    match_year = .data[["AMR match year"]], match_source = .data[["AMR match source"]]) |>
  mutate(pos = if_else(is.na(pos) & !is.na(prev_pub) & !is.na(N), round(prev_pub/100*N), pos),
         crude_prev = 100 * pos / N) |>
  filter(!is.na(crude_prev), !is.na(N), !is.na(match_pct))

cat(sprintf("Panel A: n = %d studies with both a crude prevalence and a matched estimate\n", nrow(colA)))
cat(sprintf("  Spearman rho = %.2f (p = %.2g) | N-weighted r = %.2f\n",
            cor(colA$crude_prev, colA$match_pct, method="spearman"),
            cor.test(colA$crude_prev, colA$match_pct, method="spearman")$p.value,
            { mx<-weighted.mean(colA$crude_prev,colA$N); my<-weighted.mean(colA$match_pct,colA$N)
              sum(colA$N*(colA$crude_prev-mx)*(colA$match_pct-my)) /
                sqrt(sum(colA$N*(colA$crude_prev-mx)^2)*sum(colA$N*(colA$match_pct-my)^2)) }))

## =============================================================================
## PANEL B — GBD-region colonisation (std. 2019) vs GBD 2019 attributable
##           death rate for E. coli / 3GC (drug-susceptible-infection
##           counterfactual), by GBD region
## =============================================================================
incB <- inc |> filter(!is.na(gbd_region))
colB <- pool_logit(incB, "gbd_region", ref_year = 2019) |> rename(region = group, prev2019 = prev)
colB$se_p <- (colB$hi - colB$lo) / (2*1.96); colB$w_prec <- 1/colB$se_p^2

gbd_rates <- read_csv(GBD_RATES_2019, show_col_types = FALSE) |> select(region, rate)
colB <- colB |> inner_join(gbd_rates, by = "region")   # Caribbean dropped: no 2019 GBD rate available

wcorr <- function(x, y, w) {
  mx <- weighted.mean(x, w); my <- weighted.mean(y, w)
  sum(w*(x-mx)*(y-my)) / sqrt(sum(w*(x-mx)^2) * sum(w*(y-my)^2))
}
cat(sprintf("\nPanel B: n = %d GBD regions\n", nrow(colB)))
cat(sprintf("  Spearman rho = %.2f | precision-weighted r = %.2f\n",
            cor(colB$prev2019, colB$rate, method="spearman"),
            wcorr(colB$prev2019, colB$rate, colB$w_prec)))

## ---- figures, database export, etc. would follow (plotting code omitted
##      here; see figure6_combined.png for the rendered two-panel figure) ----
