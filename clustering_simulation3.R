#=================================================================
#  dental_simulation_v2.R
#
#  R 4.5.2 – Übersetzung von dental_simulation_v2.do (Stata ≥ 16)
#
#  REVIDIERTE VERSION
#
#  Änderungen:
#  1. Systematische Variation: Praxen (10–100),
#     Patienten/Praxis (5–1000), Zeitpunkte (2–10)
#  2. Grafische + tabellarische Darstellung der Szenarien
#  3. Exposition: ausschließlich car_new
#  4. 50 Bootstrap-Replikationen (Testlauf)
#  5. Alle Modelle: Gamma-GLM, log-link
#     → exp(β) = relative Änderung im Outcome (PPD)
#
#  Drei-Ebenen-Clusterung:
#    Ebene 3: Zahnarztpraxis
#    Ebene 2: Patient
#    Ebene 1: Zahn × Zeit
#=================================================================

# ── Arbeitsverzeichnis ──
setwd("c:/arbeit/fachgesellschaften/DGoEV/2026/Impuls_ClusterkorrelierterDaten/Daten/R/")

# ── Logging (≈ Stata log using) ──
log_con <- file("log2.txt", open = "wt")
sink(log_con, type = "output", split = TRUE)

# ── Pakete laden ──
required_pkgs <- c("dplyr", "tidyr", "ggplot2", "lme4",
                   "geepack", "sandwich", "lmtest", "patchwork")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("R Version:", R.version.string, "\n")
cat("Datum:", format(Sys.time()), "\n\n")


# ═══════════════════════════════════════════════════════════════
#  TEIL 1 – SIMULATIONSFUNKTION
# ═══════════════════════════════════════════════════════════════

sim_dental <- function(practices, patients, timepoints, seed, teeth = 28L) {
  
  set.seed(seed)
  
  # --- Simulationsparameter ---
  b_car_ppd <- 0.60;  b_car_cal <- 0.80
  b_dur_ppd <- 0.20;  b_dur_cal <- 0.25
  b_car_tl  <- 0.70
  sd_prac   <- 0.40;  sd_pat  <- 0.70;  sd_tooth <- 0.30
  N <- practices * patients
  
  # ── A. PRAXIS-EBENE ──
  prac_df <- data.frame(
    practice_id = 1L:practices,
    re_prac_c   = rnorm(practices, 0, sd_prac),
    re_prac_p   = rnorm(practices, 0, sd_prac * 0.85),
    re_prac_t   = rnorm(practices, 0, sd_prac * 0.75),
    prac_urban  = rbinom(practices, 1, 0.60),
    prac_prev   = rbinom(practices, 1, 0.50)
  )
  
  # ── B. PATIENTEN-EBENE ──
  z1 <- rnorm(N); z2 <- rnorm(N); z3 <- rnorm(N)
  u_smk <- runif(N)
  
  pat_df <- data.frame(
    patient_id  = 1L:N,
    practice_id = as.integer(ceiling(seq_len(N) / patients)),
    re_pat_c    = sd_pat * z1,
    re_pat_p    = sd_pat * 0.90 * (0.500 * z1 + 0.866 * z2),
    re_pat_t    = sd_pat * 0.70 * (0.300 * z1 + 0.520 * z2 + 0.800 * z3),
    age_bl      = pmax(18, pmin(75, rnorm(N, 40, 12))),
    female      = rbinom(N, 1, 0.52),
    smoking     = ifelse(u_smk < 0.50, 0L, ifelse(u_smk < 0.75, 1L, 2L)),
    diabetes    = rbinom(N, 1, 0.10),
    ses         = sample(1L:5L, N, replace = TRUE),
    ohyg        = sample(0L:3L, N, replace = TRUE),
    bmi         = pmax(18, pmin(45, rnorm(N, 26, 5)))
  )
  rm(z1, z2, z3, u_smk)
  
  pat_df <- merge(pat_df, prac_df, by = "practice_id", all.x = TRUE)
  
  # ── C. PATIENT × ZAHN ──
  df <- pat_df[rep(seq_len(nrow(pat_df)), each = teeth), ]
  rownames(df) <- NULL
  df$tooth_num <- rep(1L:teeth, times = N)
  
  df <- df %>%
    mutate(
      ttype    = case_when(tooth_num <= 8L  ~ 1L,
                           tooth_num <= 12L ~ 2L,
                           tooth_num <= 20L ~ 3L,
                           TRUE             ~ 4L),
      jaw      = ifelse(tooth_num <= 14L, 1L, 2L),
      molar    = as.integer(ttype == 4L),
      post     = as.integer(ttype >= 3L),
      re_tooth = rnorm(n(), 0, sd_tooth),
      tooth_id = (patient_id - 1L) * teeth + tooth_num
    )
  rm(pat_df)
  
  # ── D. PATIENT × ZAHN × ZEIT ──
  n_tooth_rows <- nrow(df)
  df <- df[rep(seq_len(n_tooth_rows), each = timepoints), ]
  rownames(df) <- NULL
  df$time <- rep(0L:(timepoints - 1L), times = n_tooth_rows)
  df$age  <- df$age_bl + df$time
  df <- df[order(df$patient_id, df$tooth_num, df$time), ]
  
  # ── E. KARIES (Discrete-Time Survival) ──
  df <- df %>%
    mutate(
      xb_c = ifelse(time == 0L, -20,
                    -3.00 + 0.010 * (age_bl - 40) + 0.080 * time
                    - 0.100 * female
                    + 0.350 * (smoking == 2L) + 0.100 * (smoking == 1L)
                    + 0.200 * diabetes - 0.100 * ses + 0.200 * ohyg
                    + 0.020 * (bmi - 26)
                    + 0.350 * molar + 0.150 * (ttype == 3L)
                    + 0.100 * (jaw == 2L)
                    - 0.150 * prac_prev + 0.100 * prac_urban
                    + re_prac_c + re_pat_c + re_tooth),
      h_c = plogis(xb_c)
    ) %>%
    group_by(patient_id, tooth_num) %>%
    arrange(patient_id, tooth_num, time) %>%
    mutate(
      logS_c = cumsum(log(1 - h_c)),
      S_c    = exp(logS_c),
      U_c    = runif(1)
    ) %>%
    ungroup() %>%
    mutate(has_car = as.integer(S_c <= U_c)) %>%
    group_by(patient_id, tooth_num) %>%
    arrange(patient_id, tooth_num, time) %>%
    mutate(
      car_new = as.integer(has_car == 1L &
                             (lag(has_car, default = 0L) == 0L)),
      car_t0  = ifelse(car_new == 1L, time, NA_integer_)
    ) %>%
    tidyr::fill(car_t0, .direction = "down") %>%
    mutate(
      car_t0  = ifelse(has_car == 1L, car_t0, NA_integer_),
      dur_car = ifelse(has_car == 1L, time - car_t0, 0L)
    ) %>%
    ungroup() %>%
    select(-xb_c, -h_c, -logS_c, -S_c, -U_c)
  
  # ── F. PPD / CAL ──
  df <- df %>%
    mutate(
      mu_ppd = pmax(0.5,
                    2.50 + b_car_ppd * has_car + b_dur_ppd * dur_car
                    + 0.008 * (age - 40)
                    + 0.300 * (smoking == 2L) + 0.120 * (smoking == 1L)
                    + 0.200 * diabetes - 0.040 * ses + 0.120 * ohyg
                    + 0.200 * molar + 0.030 * time
                    + 0.30 * re_prac_p + 0.40 * re_pat_p + 0.20 * re_tooth),
      ppd = pmax(1, pmin(12,
                         rgamma(n(), shape = 15, scale = mu_ppd / 15))),
      mu_cal = pmax(0.1,
                    1.20 + b_car_cal * has_car + b_dur_cal * dur_car
                    + 0.012 * (age - 40)
                    + 0.400 * (smoking == 2L) + 0.180 * (smoking == 1L)
                    + 0.250 * diabetes - 0.050 * ses + 0.150 * ohyg
                    + 0.250 * molar + 0.050 * time
                    + 0.30 * re_prac_p + 0.50 * re_pat_p + 0.25 * re_tooth),
      cal = pmax(0, pmin(15,
                         rgamma(n(), shape = 10, scale = mu_cal / 10)))
    ) %>%
    select(-mu_ppd, -mu_cal)
  
  # ── G. PARODONTITIS ──
  df <- df %>%
    mutate(
      perio     = as.integer(cal >= 3 & ppd >= 4),
      perio_sev = case_when(
        cal >= 6 & ppd >= 5 ~ 3L,
        cal >= 4 & ppd >= 5 ~ 2L,
        cal >= 3 & ppd >= 4 ~ 1L,
        TRUE                ~ 0L)
    ) %>%
    group_by(patient_id, time) %>%
    mutate(n_perio   = sum(perio, na.rm = TRUE),
           perio_pat = as.integer(n_perio >= 2L)) %>%
    ungroup() %>%
    arrange(patient_id, tooth_num, time) %>%
    group_by(patient_id, tooth_num) %>%
    mutate(perio_inc = as.integer(
      perio == 1L & (lag(perio, default = 0L) == 0L))) %>%
    ungroup()
  
  # ── H. ZAHNVERLUST (Discrete-Time Survival) ──
  df <- df %>%
    mutate(
      xb_tl = -5.50
      + b_car_tl * has_car + 0.250 * dur_car
      + 1.000 * perio + 0.500 * (perio_sev >= 2L)
      + 0.015 * (age - 40)
      + 0.250 * (smoking == 2L) + 0.100 * (smoking == 1L)
      + 0.150 * diabetes - 0.040 * ses + 0.100 * ohyg
      + 0.150 * molar
      + re_prac_t + re_pat_t + 0.30 * re_tooth,
      h_tl = plogis(xb_tl)
    ) %>%
    group_by(patient_id, tooth_num) %>%
    arrange(patient_id, tooth_num, time) %>%
    mutate(
      logS_tl = cumsum(log(1 - h_tl)),
      S_tl    = exp(logS_tl),
      U_tl    = runif(1)
    ) %>%
    ungroup() %>%
    mutate(lost = as.integer(S_tl <= U_tl)) %>%
    group_by(patient_id, tooth_num) %>%
    arrange(patient_id, tooth_num, time) %>%
    mutate(
      tl_new = as.integer(lost == 1L &
                            (lag(lost, default = 0L) == 0L)),
      lost_t = ifelse(tl_new == 1L, time, NA_integer_)
    ) %>%
    tidyr::fill(lost_t, .direction = "down") %>%
    mutate(lost_t = ifelse(lost == 1L, lost_t, NA_integer_)) %>%
    ungroup() %>%
    mutate(
      ppd       = ifelse(lost == 1L, NA_real_, ppd),
      cal       = ifelse(lost == 1L, NA_real_, cal),
      perio     = ifelse(lost == 1L, NA_integer_, perio),
      perio_sev = ifelse(lost == 1L, NA_integer_, perio_sev)
    ) %>%
    select(-xb_tl, -h_tl, -logS_tl, -S_tl, -U_tl, -tl_new)
  
  # ── I. ABGELEITETE VARIABLEN ──
  df <- df %>%
    mutate(
      atrisk = as.integer(lost == 0L),
      smk_f  = as.integer(smoking == 1L),
      smk_c  = as.integer(smoking == 2L)
    )
  
  return(df)
}


# ── Hilfsfunktion: Cluster-Bootstrap ──
cluster_bootstrap_glm <- function(data, formula, family,
                                  cluster_var, R = 50, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  orig_mod  <- glm(formula, family = family, data = data)
  orig_coef <- coef(orig_mod)
  
  clusters <- unique(data[[cluster_var]])
  n_clust  <- length(clusters)
  
  boot_coefs <- matrix(NA_real_, nrow = R, ncol = length(orig_coef))
  colnames(boot_coefs) <- names(orig_coef)
  
  for (r in seq_len(R)) {
    sampled <- sample(clusters, n_clust, replace = TRUE)
    boot_data <- do.call(rbind, lapply(seq_along(sampled), function(j) {
      d <- data[data[[cluster_var]] == sampled[j], , drop = FALSE]
      d[[cluster_var]] <- j
      d
    }))
    mod <- tryCatch(
      suppressWarnings(glm(formula, family = family, data = boot_data)),
      error = function(e) NULL)
    if (!is.null(mod) && mod$converged)
      boot_coefs[r, ] <- coef(mod)
  }
  
  list(coefficients = orig_coef,
       se           = apply(boot_coefs, 2, sd, na.rm = TRUE),
       boot_coefs   = boot_coefs)
}


# ═══════════════════════════════════════════════════════════════
#  TEIL 2 – SZENARIOSPEZIFIKATION & ÜBERSICHT
# ═══════════════════════════════════════════════════════════════

# ── 2A: Szenario-Tabelle ──
scenarios <- tibble(
  scen_id  = 1:16,
  n_prac   = c(10L, 20L, 30L, 50L, 100L,
               rep(30L, 7),
               rep(30L, 4)),
  n_pat    = c(rep(50L, 5),
               5L, 10L, 20L, 100L, 200L, 500L, 1000L,
               rep(50L, 4)),
  n_time   = c(rep(5L, 12),
               2L, 3L, 8L, 10L),
  vary_dim = c(rep("Praxen", 5),
               rep("Patienten", 7),
               rep("Zeitpunkte", 4))
) %>%
  mutate(
    n_total   = n_prac * n_pat * 28L * n_time,
    n_analyse = n_prac * n_pat * 28L * (n_time - 1L),
    do_mixed  = as.integer(n_analyse <= 500000),
    do_gee    = as.integer(n_analyse <= 800000),
    do_boot   = as.integer(n_analyse <= 300000)
  )

# ── 2B: Tabellarische Darstellung ──
cat("\n", strrep("=", 85), "\n")
cat("  SIMULATIONSSZENARIEN - UEBERSICHT\n")
cat(strrep("=", 85), "\n")
cat("  Szenario 3 (30 Praxen, 50 Pat./Praxis, 5 ZP) = Basisszenario\n")
cat("  do_mixed/gee/boot = 0: Methode bei diesem N uebersprungen\n")
cat(strrep("=", 85), "\n\n")
print(as.data.frame(scenarios), row.names = FALSE)

saveRDS(scenarios, "scenarios.rds")

# ── 2C: Grafische Darstellung ──

# Graph 1: N Beobachtungen pro Variationsdimension
gs1 <- scenarios %>% filter(vary_dim == "Praxen") %>%
  ggplot(aes(x = factor(n_prac), y = n_total)) +
  geom_col(fill = alpha("navy", 0.7)) +
  coord_flip() +
  labs(x = "Praxen", y = "N Beobachtungen",
       title = "Praxen variiert",
       subtitle = "(50 Pat./Praxis, 5 ZP)") +
  theme_minimal()

gs2 <- scenarios %>% filter(vary_dim == "Patienten") %>%
  ggplot(aes(x = factor(n_pat), y = n_total)) +
  geom_col(fill = alpha("firebrick", 0.7)) +
  coord_flip() +
  labs(x = "Pat./Praxis", y = "N Beobachtungen",
       title = "Patienten variiert",
       subtitle = "(30 Praxen, 5 ZP)") +
  theme_minimal()

gs3 <- scenarios %>% filter(vary_dim == "Zeitpunkte") %>%
  ggplot(aes(x = factor(n_time), y = n_total)) +
  geom_col(fill = alpha("forestgreen", 0.7)) +
  coord_flip() +
  labs(x = "Zeitpunkte", y = "N Beobachtungen",
       title = "Zeitpunkte variiert",
       subtitle = "(30 Praxen, 50 Pat./Praxis)") +
  theme_minimal()

p_nobs <- (gs1 | gs2 | gs3) +
  plot_annotation(
    title   = "Datensatzgroesse nach Szenario",
    caption = "28 Zaehne pro Patient")
ggsave("fig_szenarien_nobs.png", p_nobs, width = 14, height = 5, dpi = 150)

# Graph 2: Szenarioraum (Scatter)
p_raum <- scenarios %>%
  ggplot(aes(x = n_prac, y = n_pat,
             colour = vary_dim, shape = vary_dim)) +
  geom_point(size = 5, alpha = 0.6) +
  scale_y_log10(breaks = c(5, 10, 20, 50, 100, 200, 500, 1000)) +
  scale_x_continuous(breaks = c(10, 20, 30, 50, 100)) +
  scale_colour_manual(values = c(Praxen     = "navy",
                                 Patienten  = "firebrick",
                                 Zeitpunkte = "forestgreen")) +
  scale_shape_manual(values = c(Praxen = 16, Patienten = 15,
                                Zeitpunkte = 17)) +
  labs(x = "Anzahl Praxen", y = "Patienten / Praxis (log)",
       title   = "Szenarioraum",
       colour  = "", shape = "",
       caption = "Zeitpunkte: 2-10 (im Symbol nicht dargestellt)") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("fig_szenarien_raum.png", p_raum, width = 10, height = 6, dpi = 150)


# ═══════════════════════════════════════════════════════════════
#  TEIL 3 – HAUPTSCHLEIFE: SIMULATION & ANALYSE
# ═══════════════════════════════════════════════════════════════

nscen <- nrow(scenarios)
results_list <- vector("list", nscen)

fml <- ppd ~ car_new + age + smk_f + smk_c +
  diabetes + ses + ohyg + molar

fml_mix <- ppd ~ car_new + age + smk_f + smk_c +
  diabetes + ses + ohyg + molar +
  (1 | practice_id) + (1 | patient_id)

for (i in seq_len(nscen)) {
  
  sc <- scenarios[i, ]
  np <- sc$n_prac;  nn <- sc$n_pat;  nt <- sc$n_time
  
  cat(sprintf(
    "\n%s\n  SZENARIO %d/%d: %d Praxen x %d Pat. x %d ZP\n%s\n",
    strrep("=", 50), i, nscen, np, nn, nt,
    strrep("=", 50)))
  
  t0 <- Sys.time()
  
  # ── 1. Simulieren ──
  df_sim <- sim_dental(practices = np, patients = nn,
                       timepoints = nt, seed = 20240100L + i)
  
  # ── 2. Analysedatensatz ──
  df_a <- df_sim %>% filter(atrisk == 1L, time > 0L)
  rm(df_sim); gc(verbose = FALSE)
  nobs <- nrow(df_a)
  nexp <- sum(df_a$car_new == 1L, na.rm = TRUE)
  cat(sprintf("  N Analysebeobachtungen = %s\n",
              formatC(nobs, format = "d", big.mark = ",")))
  cat(sprintf("  N car_new = 1           = %s\n",
              formatC(nexp, format = "d", big.mark = ",")))
  
  # Initialisierung
  res <- tibble(
    scen_id  = i, n_prac = np, n_pat = nn, n_time = nt,
    n_obs    = nobs, n_exp = nexp,
    icc_prac = NA_real_, icc_pat  = NA_real_,
    b_naive  = NA_real_, se_naive  = NA_real_,
    b_clpat  = NA_real_, se_clpat  = NA_real_,
    b_clprac = NA_real_, se_clprac = NA_real_,
    b_mixed  = NA_real_, se_mixed  = NA_real_,
    b_gee    = NA_real_, se_gee    = NA_real_,
    b_boot   = NA_real_, se_boot   = NA_real_)
  
  if (nexp < 20) {
    cat("  !! Zu wenige Exponierte - Szenario uebersprungen\n")
    results_list[[i]] <- res
    next
  }
  
  # ── ICC (Nullmodell, linear) ──
  cat("  -> ICC ...\n")
  tryCatch({
    mod_icc <- lmer(ppd ~ 1 + (1 | practice_id) + (1 | patient_id),
                    data = df_a, REML = FALSE,
                    control = lmerControl(
                      optCtrl = list(maxfun = 2e4)))
    vc <- as.data.frame(VarCorr(mod_icc))
    v1 <- vc$vcov[vc$grp == "practice_id"]
    v2 <- vc$vcov[vc$grp == "patient_id"]
    ve <- vc$vcov[vc$grp == "Residual"]
    vt <- v1 + v2 + ve
    res$icc_prac <- v1 / vt
    res$icc_pat  <- v2 / vt
    cat(sprintf("    ICC(Praxis)  = %.4f\n", res$icc_prac))
    cat(sprintf("    ICC(Patient) = %.4f\n", res$icc_pat))
  }, error = function(e)
    cat("    ICC-Berechnung fehlgeschlagen:", e$message, "\n"))
  
  # ── M1: Naiv (Gamma log-link, keine Clusterung) ──
  cat("  -> M1 Naiv ...\n")
  tryCatch({
    mod_glm <- glm(fml, family = Gamma(link = "log"), data = df_a)
    b_car <- coef(mod_glm)["car_new"]
    
    res$b_naive  <- b_car
    res$se_naive <- sqrt(diag(vcov(mod_glm)))["car_new"]
  }, error = function(e)
    cat("    M1 fehlgeschlagen:", e$message, "\n"))
  
  # ── M2: Cluster-robust SE auf Patientenebene ──
  cat("  -> M2 Cluster-Patient ...\n")
  tryCatch({
    if (exists("mod_glm")) {
      res$b_clpat  <- coef(mod_glm)["car_new"]
      res$se_clpat <- sqrt(diag(
        sandwich::vcovCL(mod_glm,
                         cluster = df_a$patient_id)))["car_new"]
    }
  }, error = function(e)
    cat("    M2 fehlgeschlagen:", e$message, "\n"))
  
  # ── M3: Cluster-robust SE auf Praxenebene ──
  cat("  -> M3 Cluster-Praxis ...\n")
  tryCatch({
    if (exists("mod_glm")) {
      res$b_clprac  <- coef(mod_glm)["car_new"]
      res$se_clprac <- sqrt(diag(
        sandwich::vcovCL(mod_glm,
                         cluster = df_a$practice_id)))["car_new"]
    }
  }, error = function(e)
    cat("    M3 fehlgeschlagen:", e$message, "\n"))
  
  # ── M4: Mixed (Gamma log-link, Fallback: lmer log(PPD)) ──
  if (sc$do_mixed == 1L) {
    cat("  -> M4 Mixed ...\n")
    mixed_ok <- FALSE
    tryCatch({
      mod_mx <- glmer(fml_mix, family = Gamma(link = "log"),
                      data = df_a, nAGQ = 1,
                      control = glmerControl(
                        optCtrl = list(maxfun = 1e5)))
      res$b_mixed  <- fixef(mod_mx)["car_new"]
      res$se_mixed <- sqrt(diag(vcov(mod_mx)))["car_new"]
      mixed_ok <- TRUE
      cat("    (glmer konvergiert)\n")
    }, error = function(e)
      cat("    glmer gescheitert - Fallback: lmer ln(PPD)\n"))
    if (!mixed_ok) {
      tryCatch({
        df_a$log_ppd <- log(df_a$ppd)
        fml_log <- log_ppd ~ car_new + age + smk_f + smk_c +
          diabetes + ses + ohyg + molar +
          (1 | practice_id) + (1 | patient_id)
        mod_mx <- lmer(fml_log, data = df_a, REML = FALSE)
        res$b_mixed  <- fixef(mod_mx)["car_new"]
        res$se_mixed <- sqrt(diag(vcov(mod_mx)))["car_new"]
        df_a$log_ppd <- NULL
      }, error = function(e)
        cat("    Fallback auch gescheitert:", e$message, "\n"))
    }
  } else {
    cat("  -> M4 Mixed uebersprungen (N > 500.000)\n")
  }
  
  # ── M5: GEE (Gamma log-link, exchangeable, robust SE) ──
  if (sc$do_gee == 1L) {
    cat("  -> M5 GEE ...\n")
    tryCatch({
      df_a <- df_a[order(df_a$tooth_id, df_a$time), ]
      mod_gee <- geeglm(fml, id = tooth_id, data = df_a,
                        family = Gamma(link = "log"),
                        corstr = "exchangeable")
      res$b_gee  <- coef(mod_gee)["car_new"]
      res$se_gee <- summary(mod_gee)$coefficients["car_new",
                                                  "Std.err"]
    }, error = function(e)
      cat("    GEE fehlgeschlagen:", e$message, "\n"))
  } else {
    cat("  -> M5 GEE uebersprungen (N > 800.000)\n")
  }
  
  # ── M6: Cluster-Bootstrap auf Praxenebene (50 Reps) ──
  if (sc$do_boot == 1L) {
    cat("  -> M6 Bootstrap (50 Reps) ...\n")
    tryCatch({
      bres <- cluster_bootstrap_glm(
        data = df_a, formula = fml,
        family = Gamma(link = "log"),
        cluster_var = "practice_id",
        R = 50, seed = 99900L + i)
      res$b_boot  <- bres$coefficients["car_new"]
      res$se_boot <- bres$se["car_new"]
    }, error = function(e)
      cat("    Bootstrap fehlgeschlagen:", e$message, "\n"))
  } else {
    cat("  -> M6 Bootstrap uebersprungen (N > 300.000)\n")
  }
  
  results_list[[i]] <- res
  
  t_el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  Fertig  (Dauer: %.1f Sekunden)\n", t_el))
  
  # Aufraeumen
  if (exists("mod_glm")) rm(mod_glm)
  rm(df_a); gc(verbose = FALSE)
}

results_raw <- bind_rows(results_list)
saveRDS(results_raw, "results_raw.rds")

cat("\n=== Alle Szenarien abgeschlossen ===\n")


# ═══════════════════════════════════════════════════════════════
#  TEIL 4 – ERGEBNISSE: TABELLEN & GRAFIKEN
# ═══════════════════════════════════════════════════════════════

results <- results_raw %>%
  left_join(scenarios %>% select(scen_id, vary_dim),
            by = "scen_id")

# ── Abgeleitete Groessen ──
for (m in c("naive", "clpat", "clprac", "mixed", "gee", "boot")) {
  results[[paste0("exp_b_", m)]] <- exp(results[[paste0("b_", m)]])
  results[[paste0("ci_lo_", m)]] <-
    exp(results[[paste0("b_", m)]] - 1.96 * results[[paste0("se_", m)]])
  results[[paste0("ci_hi_", m)]] <-
    exp(results[[paste0("b_", m)]] + 1.96 * results[[paste0("se_", m)]])
}
for (m in c("clpat", "clprac", "mixed", "gee", "boot")) {
  results[[paste0("ratio_", m)]] <-
    results[[paste0("se_", m)]] / results$se_naive
}
results$deff <- (results$se_clprac / results$se_naive)^2

saveRDS(results, "results_final.rds")

# ─────────────────────────────────────────────────────
#  4A. TABELLARISCHE DARSTELLUNG
# ─────────────────────────────────────────────────────

cat("\n", strrep("=", 85), "\n")
cat("  ERGEBNISSE: Gamma-GLM (log-link)\n")
cat("  Exposition:  car_new  (neue Karies, zahnspezifisch)\n")
cat("  Outcome:     PPD  (Probing Pocket Depth)\n")
cat("  exp(beta):   relative Aenderung in PPD\n")
cat(strrep("=", 85), "\n")

cat("\nTABELLE 1: exp(beta) - Punktschaetzer nach Methode\n")
cat("  exp(beta) > 1: PPD bei neuer Karies hoeher\n")
print(results %>%
        select(scen_id, n_prac, n_pat, n_time, n_exp,
               exp_b_naive, exp_b_clpat, exp_b_mixed, exp_b_gee) %>%
        as.data.frame(),
      digits = 3, row.names = FALSE)

cat("\nTABELLE 2: Standardfehler von beta(car_new)\n")
print(results %>%
        select(scen_id, n_prac, n_pat, n_time,
               se_naive, se_clpat, se_clprac,
               se_mixed, se_gee, se_boot) %>%
        as.data.frame(),
      digits = 4, row.names = FALSE)

cat("\nTABELLE 3: SE-Verhaeltnis (Methode / Naiv) und DEFF\n")
cat("  Ratio > 1: Naiver SE ist zu klein\n")
cat("  DEFF = (SE_cluster-Praxis / SE_naiv)^2\n")
print(results %>%
        select(scen_id, n_prac, n_pat, n_time,
               ratio_clpat, ratio_clprac, ratio_mixed,
               ratio_gee, ratio_boot, deff) %>%
        as.data.frame(),
      digits = 2, row.names = FALSE)

cat("\nTABELLE 4: Intraklassenkorrelationen (PPD)\n")
print(results %>%
        select(scen_id, n_prac, n_pat, n_time,
               icc_prac, icc_pat) %>%
        as.data.frame(),
      digits = 4, row.names = FALSE)

cat("\nTABELLE 5: exp(beta) mit 95%-KI - Basisszenario (Szenario 3)\n")
r3 <- results %>% filter(scen_id == 3)
for (m in c("naive", "clpat", "clprac", "mixed", "gee", "boot")) {
  cat(sprintf("  %-12s  exp(beta) = %.3f  [%.3f, %.3f]  SE = %.4f\n",
              m,
              r3[[paste0("exp_b_", m)]],
              r3[[paste0("ci_lo_", m)]],
              r3[[paste0("ci_hi_", m)]],
              r3[[paste0("se_", m)]]))
}


# ─────────────────────────────────────────────────────
#  4B. GRAFISCHE DARSTELLUNG
# ─────────────────────────────────────────────────────

# Hilfsfunktion: Ratio-Daten in Langformat
make_ratio_data <- function(df) {
  df %>%
    select(scen_id, n_prac, n_pat, n_time, vary_dim,
           ratio_clpat, ratio_clprac, ratio_mixed,
           ratio_gee, ratio_boot) %>%
    pivot_longer(cols = starts_with("ratio_"),
                 names_to = "method", values_to = "ratio",
                 names_prefix = "ratio_") %>%
    mutate(method = factor(method,
                           levels = c("clpat", "clprac", "mixed", "gee", "boot"),
                           labels = c("Cluster-Pat.", "Cluster-Praxis",
                                      "Mixed", "GEE", "Bootstrap")))
}

method_colours <- c("Cluster-Pat."   = "navy",
                    "Cluster-Praxis" = "firebrick",
                    "Mixed"          = "forestgreen",
                    "GEE"            = "orange",
                    "Bootstrap"      = "purple")
method_shapes  <- c("Cluster-Pat." = 16, "Cluster-Praxis" = 15,
                    "Mixed" = 17, "GEE" = 18, "Bootstrap" = 4)
method_lines   <- c("Cluster-Pat."   = "solid",
                    "Cluster-Praxis" = "dashed",
                    "Mixed"          = "twodash",
                    "GEE"            = "longdash",
                    "Bootstrap"      = "dotdash")

# === Graph 1: SE-Ratio vs Anzahl Praxen ===
g_r1 <- results %>%
  filter(vary_dim == "Praxen") %>%
  make_ratio_data() %>%
  ggplot(aes(x = n_prac, y = ratio,
             colour = method, shape = method, linetype = method)) +
  geom_hline(yintercept = 1, colour = "grey70", linetype = "dotted") +
  geom_line() + geom_point(size = 3) +
  scale_x_continuous(breaks = c(10, 20, 30, 50, 100)) +
  scale_colour_manual(values = method_colours) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_lines) +
  labs(x = "Anzahl Zahnarztpraxen",
       y = "SE-Verhaeltnis (Methode / Naiv)",
       title    = "A: Variation der Praxenanzahl",
       subtitle = "(50 Pat./Praxis, 5 Zeitpunkte)",
       colour = "", shape = "", linetype = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7))

# === Graph 2: SE-Ratio vs Patienten/Praxis ===
g_r2 <- results %>%
  filter(vary_dim == "Patienten" | scen_id == 3) %>%
  make_ratio_data() %>%
  ggplot(aes(x = n_pat, y = ratio,
             colour = method, shape = method, linetype = method)) +
  geom_hline(yintercept = 1, colour = "grey70", linetype = "dotted") +
  geom_line() + geom_point(size = 3) +
  scale_x_log10(breaks = c(5, 10, 20, 50, 100, 200, 500, 1000)) +
  scale_colour_manual(values = method_colours) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_lines) +
  labs(x = "Patienten pro Praxis (log-Skala)",
       y = "SE-Verhaeltnis (Methode / Naiv)",
       title    = "B: Variation der Clustergroesse",
       subtitle = "(30 Praxen, 5 Zeitpunkte)",
       colour = "", shape = "", linetype = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7))

# === Graph 3: SE-Ratio vs Zeitpunkte ===
g_r3 <- results %>%
  filter(vary_dim == "Zeitpunkte" | scen_id == 3) %>%
  make_ratio_data() %>%
  ggplot(aes(x = n_time, y = ratio,
             colour = method, shape = method, linetype = method)) +
  geom_hline(yintercept = 1, colour = "grey70", linetype = "dotted") +
  geom_line() + geom_point(size = 3) +
  scale_x_continuous(breaks = c(2, 3, 5, 8, 10)) +
  scale_colour_manual(values = method_colours) +
  scale_shape_manual(values = method_shapes) +
  scale_linetype_manual(values = method_lines) +
  labs(x = "Anzahl Messwiederholungen",
       y = "SE-Verhaeltnis (Methode / Naiv)",
       title    = "C: Variation der Zeitpunkte",
       subtitle = "(30 Praxen, 50 Pat./Praxis)",
       colour = "", shape = "", linetype = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 7))

# === Kombinierter Graph ===
p_se <- g_r1 / g_r2 / g_r3 +
  plot_annotation(
    title    = "Design-Effekte auf Standardfehler",
    subtitle = "Gamma-GLM (log-link); Exposition = car_new; Outcome = PPD",
    caption  = paste(
      "SE-Verhaeltnis > 1: Naives Modell unterschaetzt den SE |",
      "50 Bootstrap-Replikationen (Testlauf)"))
ggsave("fig_se_ratio_combined.png", p_se,
       width = 8, height = 14, dpi = 150)

# === Graph 4: Design-Effekt vs Clustergroesse ===
p_deff <- results %>%
  filter(vary_dim == "Patienten" | scen_id == 3) %>%
  ggplot(aes(x = n_pat, y = deff)) +
  geom_hline(yintercept = 1, colour = "grey70", linetype = "dotted") +
  geom_line(colour = "firebrick", linewidth = 0.8) +
  geom_point(colour = "firebrick", size = 3) +
  scale_x_log10(breaks = c(5, 10, 20, 50, 100, 200, 500, 1000)) +
  labs(x = "Patienten pro Praxis (log-Skala)",
       y = "Design-Effekt (DEFF)",
       title    = "Design-Effekt vs. Clustergroesse",
       subtitle = expression(DEFF == (SE[cluster] / SE[naiv])^2),
       caption  = "DEFF > 1: Effektive Stichprobe kleiner als nominal") +
  theme_minimal()
ggsave("fig_deff_clustersize.png", p_deff,
       width = 10, height = 6, dpi = 150)

# === Graph 5: ICC nach Szenario ===
p_icc <- results %>%
  filter(vary_dim == "Patienten" | scen_id == 3) %>%
  select(n_pat, icc_prac, icc_pat) %>%
  pivot_longer(-n_pat, names_to = "level", values_to = "icc",
               names_prefix = "icc_") %>%
  mutate(level = factor(level,
                        levels = c("prac", "pat"),
                        labels = c("ICC(Praxis)", "ICC(Patient)"))) %>%
  ggplot(aes(x = n_pat, y = icc, colour = level, shape = level)) +
  geom_line() + geom_point(size = 3) +
  scale_x_log10(breaks = c(5, 10, 20, 50, 100, 200, 500, 1000)) +
  scale_colour_manual(values = c("ICC(Praxis)"  = "navy",
                                 "ICC(Patient)" = "firebrick")) +
  labs(x = "Patienten pro Praxis (log-Skala)", y = "ICC",
       title  = "Intraklassenkorrelation vs. Clustergroesse",
       colour = "", shape = "") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("fig_icc_clustersize.png", p_icc,
       width = 10, height = 6, dpi = 150)

# === Graph 6: exp(beta) Forest Plot - Basisszenario ===
forest_df <- r3 %>%
  select(starts_with("exp_b_"), starts_with("ci_lo_"),
         starts_with("ci_hi_"), starts_with("se_")) %>%
  pivot_longer(
    everything(),
    names_to  = c(".value", "method"),
    names_pattern = "^(exp_b|ci_lo|ci_hi|se)_(.+)$"
  ) %>%
  rename(exp_b = exp_b, ci_lo = ci_lo, ci_hi = ci_hi, se = se) %>%
  filter(!is.na(exp_b)) %>%
  mutate(method_label = factor(method,
                               levels = c("boot", "gee", "mixed", "clprac", "clpat", "naive"),
                               labels = c("Cluster-Bootstrap (Praxis)",
                                          "GEE (exchangeable)",
                                          "Mixed Model (glmer)",
                                          "Cluster-robust (Praxis)",
                                          "Cluster-robust (Patient)",
                                          "Naiv (keine Clusterung)")))

p_forest <- ggplot(forest_df,
                   aes(x = exp_b, y = method_label,
                       xmin = ci_lo, xmax = ci_hi)) +
  geom_vline(xintercept = 1, colour = "grey70", linetype = "dotted") +
  geom_errorbarh(height = 0.2, colour = "navy") +
  geom_point(colour = "firebrick", size = 4) +
  labs(x = "exp(beta): Relative Aenderung in PPD", y = "",
       title    = "Basisszenario: exp(beta) mit 95%-KI nach Methode",
       subtitle = "30 Praxen, 50 Pat./Praxis, 5 Zeitpunkte",
       caption  = "Gestrichelt: exp(beta) = 1 (kein Effekt)") +
  theme_minimal()
ggsave("fig_forest_base.png", p_forest,
       width = 10, height = 5, dpi = 150)


# ═══════════════════════════════════════════════════════════════
#  TEIL 5 – ICC & DESIGN-EFFEKTE: BASISSZENARIO VERTIEFT
# ═══════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("  VERTIEFUNG: ICC & DESIGN-EFFEKTE (Basisszenario)\n")
cat(strrep("=", 70), "\n")

# Basisszenario neu simulieren
df_base <- sim_dental(practices = 30, patients = 50,
                      timepoints = 5, seed = 20240103)
saveRDS(df_base, "sim_base.rds")

df_b <- df_base %>% filter(atrisk == 1L, time > 0L)
rm(df_base); gc(verbose = FALSE)

# ── 5A: ICC aus Nullmodell ──
cat("\n--- ICC fuer PPD (Nullmodell) ---\n")
mod_null <- lmer(ppd ~ 1 + (1 | practice_id) + (1 | patient_id),
                 data = df_b, REML = FALSE)
print(summary(mod_null))

vc     <- as.data.frame(VarCorr(mod_null))
v_prac <- vc$vcov[vc$grp == "practice_id"]
v_pat  <- vc$vcov[vc$grp == "patient_id"]
v_res  <- vc$vcov[vc$grp == "Residual"]
v_tot  <- v_prac + v_pat + v_res

icc_prac <- v_prac / v_tot
icc_pat  <- v_pat  / v_tot
icc_both <- (v_prac + v_pat) / v_tot

# ── 5B: Varianzkomponenten & Design-Effekt ──

# Mittlere Clustergroessen
m_prac <- df_b %>% group_by(practice_id) %>%
  summarise(n = n(), .groups = "drop") %>% pull(n) %>% mean()
m_pat  <- df_b %>% group_by(patient_id) %>%
  summarise(n = n(), .groups = "drop") %>% pull(n) %>% mean()

# Design-Effekte: DEFF = 1 + (m - 1) * ICC
deff_prac <- 1 + (m_prac - 1) * icc_prac
deff_pat  <- 1 + (m_pat  - 1) * icc_pat
deff_both <- 1 + (m_prac - 1) * icc_both

cat("\n--- Varianzkomponenten ---\n")
cat(sprintf("  Var(Praxis)   = %.4f\n", v_prac))
cat(sprintf("  Var(Patient)  = %.4f\n", v_pat))
cat(sprintf("  Var(Residual) = %.4f\n", v_res))

cat("\n--- ICC ---\n")
cat(sprintf("  ICC(Praxis)          = %.4f\n", icc_prac))
cat(sprintf("  ICC(Patient)         = %.4f\n", icc_pat))
cat(sprintf("  ICC(Praxis+Patient)  = %.4f\n", icc_both))

cat("\n--- Mittlere Clustergroessen ---\n")
cat(sprintf("  m(Praxis)  = %.0f\n", m_prac))
cat(sprintf("  m(Patient) = %.0f\n", m_pat))

cat("\n--- Design-Effekte ---\n")
cat(sprintf("  DEFF(Praxis)          = %.2f\n", deff_prac))
cat(sprintf("  DEFF(Patient)         = %.2f\n", deff_pat))
cat(sprintf("  DEFF(Praxis+Patient)  = %.2f\n", deff_both))

cat("\n--- SE-Inflationsfaktoren ---\n")
cat(sprintf("  sqrt(DEFF(Praxis))         = %.3f\n", sqrt(deff_prac)))
cat(sprintf("  sqrt(DEFF(Praxis+Patient)) = %.3f\n", sqrt(deff_both)))
cat(sprintf("\n  -> Naive SE um Faktor %.2f zu klein,\n", sqrt(deff_both)))
cat("    wenn Clusterung auf beiden Ebenen ignoriert wird.\n")

# ── 5C: Vollstaendige Modellergebnisse (Basisszenario) ──
cat("\n--- Detaillierte Modellvergleiche ---\n")

# M1: Naiv
mod_naive <- glm(fml, family = Gamma(link = "log"), data = df_b)
b_all  <- coef(mod_naive)["car_new"]
se_nv  <- sqrt(diag(vcov(mod_naive)))["car_new"]

# M2: Cluster-Patient SE
se_cp  <- sqrt(diag(
  sandwich::vcovCL(mod_naive, cluster = df_b$patient_id)))["car_new"]

# M3: Cluster-Praxis SE
se_cpr <- sqrt(diag(
  sandwich::vcovCL(mod_naive, cluster = df_b$practice_id)))["car_new"]

# M4: Mixed
mod_mixed <- tryCatch(
  glmer(fml_mix, family = Gamma(link = "log"),
        data = df_b, nAGQ = 1,
        control = glmerControl(optCtrl = list(maxfun = 1e5))),
  error = function(e) {
    cat("  glmer gescheitert - Fallback: lmer ln(PPD)\n")
    df_b$log_ppd <<- log(df_b$ppd)
    fml_log2 <- log_ppd ~ car_new + age + smk_f + smk_c +
      diabetes + ses + ohyg + molar +
      (1 | practice_id) + (1 | patient_id)
    lmer(fml_log2, data = df_b, REML = FALSE)
  })
b_mx  <- fixef(mod_mixed)["car_new"]
se_mx <- sqrt(diag(vcov(mod_mixed)))["car_new"]

# M5: GEE
df_b <- df_b[order(df_b$tooth_id, df_b$time), ]
mod_gee <- geeglm(fml, id = tooth_id, data = df_b,
                  family = Gamma(link = "log"),
                  corstr = "exchangeable")
b_ge  <- coef(mod_gee)["car_new"]
se_ge <- summary(mod_gee)$coefficients["car_new", "Std.err"]

# M6: Bootstrap
bres <- cluster_bootstrap_glm(
  df_b, fml, Gamma(link = "log"), "practice_id",
  R = 50, seed = 99903)
b_bt  <- bres$coefficients["car_new"]
se_bt <- bres$se["car_new"]

# Vergleichstabelle
comp <- tibble(
  Methode = c("Naiv", "Cluster-Pat.", "Cluster-Praxis",
              "Mixed", "GEE", "Bootstrap"),
  b       = c(b_all,  b_all,  b_all,  b_mx,  b_ge,  b_bt),
  SE      = c(se_nv,  se_cp,  se_cpr, se_mx, se_ge, se_bt),
  exp_b   = exp(c(b_all, b_all, b_all, b_mx, b_ge, b_bt)),
  ci_lo   = exp(c(b_all, b_all, b_all, b_mx, b_ge, b_bt) -
                  1.96 * c(se_nv, se_cp, se_cpr, se_mx, se_ge, se_bt)),
  ci_hi   = exp(c(b_all, b_all, b_all, b_mx, b_ge, b_bt) +
                  1.96 * c(se_nv, se_cp, se_cpr, se_mx, se_ge, se_bt))
)

cat("\nbeta-Koeffizienten (log-Skala) und exp(beta):\n")
print(as.data.frame(comp), digits = 4, row.names = FALSE)

cat("\nexp(beta) - Relative Aenderung in PPD bei neuer Karies:\n")
for (j in seq_len(nrow(comp))) {
  cat(sprintf("  %-16s  exp(beta) = %.3f  [%.3f, %.3f]  SE = %.4f\n",
              comp$Methode[j], comp$exp_b[j],
              comp$ci_lo[j], comp$ci_hi[j], comp$SE[j]))
}


# ═══════════════════════════════════════════════════════════════
#  TEIL 6 – ZUSAMMENFASSUNG & INTERPRETATION
# ═══════════════════════════════════════════════════════════════

cat("\n", strrep("=", 70), "\n")
cat("  ZUSAMMENFASSUNG & INTERPRETATION\n")
cat(strrep("=", 70), "\n")

cat("
Modell:
  Gamma-GLM mit Log-Link
  log(E[PPD | X]) = b0 + b1*car_new + b2*age + ...
  -> exp(b1) = relative Aenderung in E[PPD]

Exposition:
  car_new = 1 genau zum Zeitpunkt der neuen Karies
  (danach wieder 0 -> erfasst akuten Onset-Effekt)

Interpretation von exp(beta):
  exp(beta) = 1.05:  PPD ist 5% hoeher bei Zaehnen
                     mit neuer Karies (zum Onset-Zeitpunkt)
  exp(beta) = 1.00:  kein Effekt

Erwartete Beobachtungen:

  1. SE-INFLATION BEI IGNORIEREN DER CLUSTERUNG
     Ratio SE(Cluster) / SE(Naiv) > 1 in allen Szenarien.
     Je groesser die Cluster (mehr Patienten/Praxis),
     desto staerker die SE-Inflation.
     -> Naive p-Werte und KIs sind anti-konservativ!

  2. MIXED VS. GEE (Log-Link)
     Log-Link ist kollabierbar:
     -> beta(Mixed) ca. beta(GEE) fuer car_new
     -> Hauptunterschied in den SE, nicht in exp(beta)
     (Bei Logit-Link waere: OR_konditional > OR_marginal)

  3. ANZAHL CLUSTER (Praxen)
     Wenige Praxen -> Sandwich-SE nach unten verzerrt
     -> Cluster-Bootstrap liefert konservativere SE
     Faustregel: g >= 30-50 Cluster fuer zuverlaessige
     Sandwich-Schaetzer

  4. DESIGN-EFFEKT
     DEFF = 1 + (m - 1) * ICC
     Grosse Cluster + hoher ICC -> grosser DEFF
     -> Effektive Stichprobe = N / DEFF

Wahre Simulationsparameter (DGP):
  Karies -> PPD: +0.60 mm  (identity-Skala)
  Kariesdauer -> PPD: +0.20 mm/Jahr
  SD: Praxis=0.40, Patient=0.70, Zahn=0.30

  Hinweis: car_new erfasst nur den akuten Onset.
  Im DGP wirkt has_car (kumulativ), daher ist
  exp(beta_car_new) nicht direkt mit beta_DGP=0.60
  vergleichbar.
")

cat("\n", strrep("=", 70), "\n")
cat("  DATEIEN\n")
cat(strrep("=", 70), "\n")
cat("
  scenarios.rds         - Szenariospezifikation
  results_raw.rds       - Rohergebnisse (Koeffizienten, SE)
  results_final.rds     - Ergebnisse + exp(beta), Ratios, DEFF
  sim_base.rds          - Basisszenario-Daten
  fig_szenarien_*.png   - Szenario-Uebersichtsgrafiken
  fig_se_ratio_*.png    - SE-Verhaeltnis-Grafiken
  fig_deff_*.png        - Design-Effekt-Grafik
  fig_icc_*.png         - ICC-Grafik
  fig_forest_*.png      - Forest Plot (Basisszenario)
")
cat(strrep("=", 70), "\n")

# ── Log schliessen ──
sink(type = "output")
close(log_con)

cat("Fertig. Log gespeichert in log2.txt\n")
```