##########################################################################
# WASH & 3GCRE colonisation                                              #
# Figure 6 — 3GCR-E. coli community colonisation vs GBD 2021 E. coli/3GC #
##########################################################################

# install/load packages
pacman::p_load(readxl, dplyr, tidyr, stringr, purrr, ggplot2, patchwork)

#import data
raw <- read_excel("ESBL_CRE_FecalColonization_HealthyIndividuals_GlobalEvidence_v12.xlsx", sheet = "Primary Studies") |> rename_with(str_squish)

# add GBD super-regions
SR <- c("Central Europe, Eastern Europe, and Central Asia","High-income",
        "Latin America and Caribbean","North Africa and Middle East","South Asia",
        "Southeast Asia, East Asia, and Oceania","Sub-Saharan Africa")
gbd_sr <- function(c) {
  sets <- list(
    "High-income"=c("Australia","Belgium","Cyprus","Denmark","Finland","France","Germany",
      "Netherlands","Norway","Portugal","Spain","Sweden","Switzerland","Japan",
      "Republic of Korea","Singapore","Taiwan","USA","Chile","Multiple (Europe/USA)"),
    "Central Europe, Eastern Europe, and Central Asia"=c("Hungary","Turkey"),
    "North Africa and Middle East"=c("Iran","Qatar","Saudi Arabia","Tunisia"),
    "South Asia"=c("Bangladesh","India","Nepal"),
    "Southeast Asia, East Asia, and Oceania"=c("China","China (Central South)","Cambodia",
      "Indonesia","Laos","Thailand","Vietnam"),
    "Latin America and Caribbean"=c("Brazil","Colombia","Ecuador","French Guiana","Guatemala",
      "Mexico","Venezuela"))
  hit <- names(sets)[map_lgl(sets, ~ c %in% .x)]
  if (length(hit)) hit[1] else "Sub-Saharan Africa"        # default
}
mid_year <- function(s) {                                 # mean of any 4-digit years in the cell
  yrs <- regmatches(as.character(s), gregexpr("(19|20)[0-9]{2}", as.character(s)))
  vapply(yrs, function(y) if (length(y)) mean(as.numeric(y)) else NA_real_, numeric(1))
}


## ---- load, clean, exclude to the 147-study core ---------------------------
pick <- function(p) names(raw)[startsWith(names(raw), p)][1]

d <- tibble(
    country = raw[[pick("Country")]],
    year    = mid_year(raw[[pick("Year(s) Sample")]]),
    N       = suppressWarnings(as.numeric(raw[[pick("N Stool")]])),
    pos     = suppressWarnings(as.numeric(raw[[pick("N Positive ESBL/3GCR")]])),
    prev    = suppressWarnings(as.numeric(raw[[pick("Prevalence ESBL")]])),
    excl    = raw[[pick("Excluded")]]) |>
  mutate(pos = if_else(is.na(pos) & !is.na(prev), round(prev/100*N), pos)) |>
  filter(is.na(excl),
         !str_detect(country, regex("travel|abroad|healthcare students", TRUE)),
         !is.na(year), !is.na(N), !is.na(pos), N > 0, pos >= 0, pos <= N) |>
  mutate(super_region = map_chr(country, gbd_sr),
         p_hat = (pos + .5)/(N + 1),
         lp = qlogis(p_hat),                                # empirical logit
         v  = 1/(pos + .5) + 1/(N - pos + .5))              # logit variance

## ---- colonisation pooled to super-region, standardised to 2021 ------------
## WLS logit meta-regression: 7 super-region intercepts + shared year slope,
## DerSimonian-Laird between-study heterogeneity.
d$yc <- d$year - 2021
X <- cbind(sapply(SR, \(r) as.numeric(d$super_region == r)), yc = d$yc)
y <- d$lp; v <- d$v
wls <- function(X, y, w) { cov <- solve(t(X) %*% (X * w)); list(b = cov %*% (t(X) %*% (y * w)), cov = cov) }
b0  <- wls(X, y, 1/v)$b
res <- y - X %*% b0; w <- 1/v
tau2 <- max(0, (sum(res^2/v) - (length(y) - ncol(X))) / (sum(w) - sum(w^2)/sum(w)))
fit <- wls(X, y, 1/(v + tau2)); beta <- as.numeric(fit$b); se <- sqrt(diag(fit$cov))

col <- tibble(super_region = SR,
              prev2021 = plogis(beta[1:7]),
              p_lo = plogis(beta[1:7] - 1.96*se[1:7]),
              p_hi = plogis(beta[1:7] + 1.96*se[1:7]),
              n = as.integer(table(factor(d$super_region, SR))))

## per-super-region decadal odds trend: DL-weighted logit slope
## (NA if <6 studies or <8-year span)
slope_sr <- function(r) {
  g <- filter(d, super_region == r)
  if (nrow(g) < 6 || diff(range(g$year)) < 8) return(NA_real_)
  Xg <- cbind(1, g$year - 2021); yg <- g$lp; vg <- g$v
  b  <- solve(t(Xg) %*% (Xg/vg)) %*% (t(Xg) %*% (yg/vg))
  rr <- yg - Xg %*% b; wg <- 1/vg
  t2 <- max(0, (sum(rr^2/vg) - (nrow(g) - 2)) / (sum(wg) - sum(wg^2)/sum(wg)))
  b2 <- solve(t(Xg) %*% (Xg/(vg + t2))) %*% (t(Xg) %*% (yg/(vg + t2)))
  exp(b2[2] * 10)
}
col$odds_dec <- map_dbl(SR, slope_sr)

## ---- GBD 2021 AMR: E. coli/3GC attributable death rate --------------------
# extracted from Supplement to: GBD 2021 Antimicrobial Resistance Collaborators. Global burden of
# bacterial antimicrobial resistance 1990–2021: a systematic analysis with forecasts
# to 2050. Lancet 2024; published online Sept 16. https://doi.org/10.1016/S0140-6736(24)01867-1
## Tables S13 (1990) / S14 (2021) attributable death counts & populations back-derived
## from Table 3 associated deaths (count / rate * 1e5).
gbd <- tibble(super_region = SR,
  d1990 = c(1470,1500,1010,1030,10200,3130,3410),
  d2021 = c(2170,3480,2050,1710, 9650,7980,6050),
  a90c = c(285,477,247,264,1400,1110,990), a90r = c(67.8,52.5,63.3,78.0,128,65.8,202),
  a21c = c(265,553,322,226,1260,1150,923),  a21r = c(63.4,50.7,54.2,36.3,68.5,52.7,81.5)) |>
  mutate(rate1990 = d1990/(a90c*1e3/a90r*1e5)*1e5,
         rate2021 = d2021/(a21c*1e3/a21r*1e5)*1e5,
         rate_chg = (rate2021 - rate1990)/rate1990)

f <- col |> left_join(gbd, "super_region") |>
  mutate(w_prec = (2*1.96/(p_hi - p_lo))^2)              # prevalence-scale precision weight

## ---- correlations: precision-weighted (primary) + Spearman (sensitivity) --
wcorr <- function(x, y, w) { mx <- weighted.mean(x, w); my <- weighted.mean(y, w)
  sum(w*(x-mx)*(y-my)) / sqrt(sum(w*(x-mx)^2) * sum(w*(y-my)^2)) }
size_w  <- with(f, wcorr(prev2021, rate2021, w_prec))
size_s  <- with(f, cor(prev2021, rate2021, method = "spearman"))
tt <- filter(f, !is.na(odds_dec))
trend_w <- with(tt, wcorr(odds_dec, rate_chg, w_prec))
trend_s <- with(tt, cor(odds_dec, rate_chg, method = "spearman"))
cat(sprintf("SIZE  weighted r=%.2f (Spearman rho=%.2f)\n", size_w, size_s))
cat(sprintf("TREND weighted r=%.2f (Spearman rho=%.2f)\n", trend_w, trend_s))
write.csv(f, "figure6_data.csv", row.names = FALSE)

## ---- Figure 6 -------------------------------------------------------------
ab <- setNames(c("C/E Eur & C Asia","High-income","Lat Am & Carib","N Afr & M East",
                 "South Asia","SE/E Asia & Ocean","Sub-Saharan Afr"), SR)
f$lab <- paste0(ab[f$super_region], " (n=", f$n, ")"); tt$lab <- paste0(ab[tt$super_region], " (n=", tt$n, ")")

f6a <- ggplot(f, aes(prev2021*100, rate2021)) +
  geom_errorbarh(aes(xmin = p_lo*100, xmax = p_hi*100), colour = "grey80", height = 0) +
  geom_point(aes(colour = super_region, size = n)) + geom_text(aes(label = lab), size = 2.7, hjust = -.1, vjust = -.3) +
  scale_size(range = c(2, 11)) + guides(colour = "none", size = "none") +
  labs(x = "3GCR-E colonisation prevalence, 2021 (%)",
       y = "E. coli/3GC attributable death rate, 2021 (/100k)",
       subtitle = sprintf("A  SIZE   weighted r=%.2f (Spearman rho=%.2f)", size_w, size_s)) +
  theme_minimal()
f6b <- ggplot(tt, aes(odds_dec, rate_chg*100)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(aes(colour = super_region, size = n)) + geom_text(aes(label = lab), size = 2.7, hjust = -.1, vjust = -.3) +
  scale_size(range = c(2, 11)) + guides(colour = "none", size = "none") +
  labs(x = "Colonisation trend (odds x per decade; >1 = rising)", y = "Death-rate change 1990-2021 (%)",
       subtitle = sprintf("B  TREND   weighted r=%.2f (Spearman rho=%.2f)", trend_w, trend_s)) +
  theme_minimal()
fig6 <- f6a | f6b
ggsave("figure6.png", fig6, width = 13.5, height = 5.8, dpi = 200)
