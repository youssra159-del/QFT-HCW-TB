###############################################################################
#  Latent tuberculosis infection (LTBI) among healthcare workers
#  Complete statistical analysis of QuantiFERON (QFT) data – from A to Z
#
#  Research questions
#   Q1. What is the prevalence of LTBI (QFT positive) in this population?
#   Q2. Which factors are associated with QFT positivity (bivariate analysis)?
#   Q3. Which factors remain associated after adjustment (multivariable
#       logistic regression), and how reliable is the model?
#
#  Methods: a priori variable
#  selection, events per parameter, multiple imputation, Firth regression,
#  bootstrap internal validation.

###############################################################################


## 0. SETUP ====================================================================
pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "forcats", "ggplot2", "patchwork",
          "broom", "car", "pROC", "rms", "brglm2", "mice", "openxlsx",
          "flextable", "officer")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

data_file <- "Data_QFT_Youssra.xlsx"   # adapt the path if needed
out_dir   <- "results"
fig_dir   <- file.path(out_dir, "figures")
supp_dir  <- file.path(fig_dir, "supplementary")
tab_dir   <- file.path(out_dir, "tables")
for (p in c(out_dir, fig_dir, supp_dir, tab_dir)) dir.create(p, showWarnings = FALSE, recursive = TRUE)

SEED <- 2026        # fixed seed: Monte Carlo p-values, bootstrap and imputation
options(digits = 4)

# Formatting helpers (2 decimals for OR, 3 for p-values)
f1  <- function(x) formatC(x, format = "f", digits = 1)
f2  <- function(x) formatC(x, format = "f", digits = 2)
fp  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3)))
fOR <- function(or, lo, hi) sprintf("%s (%s–%s)", f2(or), f2(lo), f2(hi))
miss_txt <- function(m) if (m > 0) sprintf(" (%d missing)", m) else ""


## 1. IMPORT ===================================================================
raw <- read_excel(data_file, na = c("", "NA"))     # "NA" typed as text = missing
names(raw) <- c("id", "age", "sex", "occupation", "department", "tb_exposure",
                "experience", "mask", "biosafety_training", "bcg", "bcg_scar",
                "smoking", "qft_raw", "history", "practices")
cat("Records imported:", nrow(raw), "\n")


## 2. DATA CLEANING AND RECODING ===============================================
# Every recoding decision is explicit here, so that it can be checked.
sq <- function(x) str_squish(as.character(x))
yes_no <- function(x) {
  x <- str_to_lower(sq(x))
  case_when(x %in% c("yes", "oui", "rarement") ~ "Yes",   # one answer "rarement" = exposed
            x %in% c("no", "n", "non")         ~ "No",
            TRUE ~ NA_character_)                          # e.g. "RAS" = not interpretable
}

d <- raw %>%
  mutate(
    across(everything(), ~ if (is.character(.x)) sq(.x) else .x),
    age        = as.numeric(age),
    experience = suppressWarnings(as.numeric(experience)),

    # Outcome: QFT positive (1) / negative (0); indeterminate = missing
    qft_pos = case_when(qft_raw == "QFT+" ~ 1L, qft_raw == "QFT-" ~ 0L, TRUE ~ NA_integer_),
    sex = factor(sex, levels = c("Female", "Male")),

    # Occupation: 5 clinically homogeneous groups (reference = Nurse)
    occupation_grp = case_when(
      str_detect(occupation, regex("doctor|médcin|docteur|intern", TRUE))               ~ "Physician / intern",
      str_detect(occupation, regex("^nurse", TRUE))                                     ~ "Nurse",
      str_detect(occupation, regex("nursing assistant|ambulancier", TRUE))              ~ "Nursing assistant / EMT",
      str_detect(occupation, regex("laboratory|scientifique|technical|pharmacien", TRUE)) ~ "Laboratory / pharmacy",
      !is.na(occupation) ~ "Non-clinical staff",
      TRUE ~ NA_character_),
    occupation_grp = factor(occupation_grp, levels = c("Nurse", "Physician / intern",
      "Nursing assistant / EMT", "Laboratory / pharmacy", "Non-clinical staff")),

    # Department: 6 groups (reference = Intensive care)
    department_grp = case_when(
      str_detect(department, regex("intensive care", TRUE))                   ~ "Intensive care",
      str_detect(department, regex("pulmonology", TRUE))                      ~ "Pulmonology",
      str_detect(department, regex("emergency|resuscitation", TRUE))          ~ "Emergency / EMS",
      str_detect(department, regex("laboratory|genetics", TRUE))              ~ "Laboratory",
      str_detect(department, regex("consultation|internal|ophthalmo", TRUE))  ~ "Medicine / outpatient",
      !is.na(department) ~ "Administration / security",
      TRUE ~ NA_character_),
    department_grp = factor(department_grp, levels = c("Intensive care", "Pulmonology",
      "Emergency / EMS", "Laboratory", "Medicine / outpatient", "Administration / security")),

    tb_exposure        = factor(yes_no(tb_exposure),        levels = c("No", "Yes")),
    biosafety_training = factor(yes_no(biosafety_training), levels = c("Yes", "No")),
    bcg                = factor(yes_no(bcg),                levels = c("Yes", "No")),
    bcg_scar           = factor(yes_no(bcg_scar),           levels = c("Yes", "No"),
                                labels = c("Present", "Absent")),

    # Mask use: ordinal; "NeveRarely" = rarely; "P" = missing; never + rarely merged
    mask4 = case_when(mask %in% c("Never", "Rarely", "NeveRarely") ~ "Never / rarely",
                      mask == "Sometimes" ~ "Sometimes",
                      mask == "Often"     ~ "Often",
                      mask == "Always"    ~ "Always",
                      TRUE ~ NA_character_),
    mask4 = factor(mask4, levels = c("Always", "Often", "Sometimes", "Never / rarely")),
    mask_regular = factor(ifelse(is.na(mask4), NA,
                          ifelse(mask4 %in% c("Always", "Often"), "Regular", "Irregular")),
                          levels = c("Regular", "Irregular")),

    # Smoking: current or former smoker = "Yes"
    smoking = case_when(str_to_lower(smoking) == "no" ~ "No",
                        is.na(smoking) ~ NA_character_,
                        TRUE ~ "Yes"),
    smoking = factor(smoking, levels = c("No", "Yes")),

    # Medical history (free text; blank = none reported)
    history = ifelse(is.na(history) | str_to_upper(history) == "NA", "", history),
    tb_history    = factor(ifelse(str_detect(history, regex("tubercul", TRUE)), "Yes", "No"), c("No", "Yes")),
    covid_history = factor(ifelse(str_detect(history, regex("covid", TRUE)), "Yes", "No"), c("No", "Yes")),

    # High-risk practices (free text; blank = none reported)
    practices = ifelse(is.na(practices), "", practices),
    agp = factor(ifelse(str_detect(practices,
            regex("intubation|suction|aspiration|bronchoscop|sputum|expector|exprector|r[ée]animation", TRUE)),
            "Yes", "No"), c("No", "Yes")),                     # aerosol-generating procedures
    specimen = factor(ifelse(str_detect(practices, regex("sample|sampling|prél", TRUE)), "Yes", "No"), c("No", "Yes")),

    # Age and experience groups (description only)
    age_group = cut(age, c(-Inf, 24, 29, Inf), labels = c("< 25 years", "25–29 years", "≥ 30 years")),
    exp_group = cut(experience, c(-Inf, 1, 4, Inf), labels = c("≤ 1 year", "> 1 to 4 years", "> 4 years"))
  )


## 3. DATA QUALITY CHECKS ======================================================
cat("\nDuplicate identifiers:", sum(duplicated(d$id)), "\n")
cat("Indeterminate QFT results (excluded):", sum(is.na(d$qft_pos)), "\n")
d$prefix <- str_extract(d$id, "(?<=-)[A-Z]+")          # ITL = LTBI case, TS = healthy control
print(table(Prefix = d$prefix, QFT = d$qft_raw, useNA = "ifany"))
print(d %>% filter((prefix %in% c("ITL", "ILT") & qft_raw != "QFT+") |
                   (prefix == "TS" & qft_raw != "QFT-")) %>% select(id, qft_raw))
print(summary(d[, c("age", "experience")]))
cat("Experience > age - 16 years:", sum(d$experience > d$age - 16, na.rm = TRUE), "\n")

# Analysis population: interpretable QFT result
a <- d %>% filter(!is.na(qft_pos)) %>% mutate(QFT = factor(qft_pos, 0:1, c("QFT-", "QFT+")))
n_tot <- nrow(a); n_pos <- sum(a$qft_pos); n_neg <- n_tot - n_pos
cat("Analysis population: n =", n_tot, "(QFT+ =", n_pos, ")\n")

labels_en <- c(age = "Age (years)", experience = "Professional experience (years)",
  sex = "Sex", age_group = "Age group", occupation_grp = "Occupation", department_grp = "Department",
  tb_exposure = "Self-reported TB exposure", exp_group = "Professional experience",
  mask4 = "Mask use", mask_regular = "Regular mask use", biosafety_training = "Biosafety training",
  bcg = "BCG vaccination", bcg_scar = "BCG scar", smoking = "Current or former smoker",
  tb_history = "History of tuberculosis", covid_history = "History of COVID-19",
  agp = "Aerosol-generating procedures", specimen = "Specimen handling")

# Supplementary Table S1 – missing data by variable
vars_missing <- c("bcg_scar", "biosafety_training", "mask4", "experience", "bcg", "tb_exposure",
                  "smoking", "occupation_grp", "department_grp", "age", "sex")
table_S1 <- data.frame(Variable = unname(labels_en[vars_missing]),
                       `Missing, n` = sapply(a[vars_missing], function(x) sum(is.na(x))),
                       check.names = FALSE) %>%
  mutate(`Missing, %` = f1(100 * `Missing, n` / n_tot))
rownames(table_S1) <- NULL
print(table_S1)


## 4. DESCRIPTIVE ANALYSIS – TABLE 2 ==========================================
# Normality of continuous variables
for (v in c("age", "experience")) {
  x <- a[[v]][!is.na(a[[v]])]
  cat(sprintf("%s: mean %.1f, median %.1f, skewness %.2f, Shapiro-Wilk p = %.4f\n", v, mean(x),
              median(x), mean((x - mean(x))^3) / sd(x)^3, shapiro.test(x)$p.value))
}

desc_cat <- function(var, label, levels = NULL, single = NULL) {
  x <- a[[var]]; tab <- table(x); tot <- sum(tab); m <- sum(is.na(x))
  if (!is.null(single)) {   # binary variable: one line for the category of interest
    return(data.frame(Characteristic = paste0(label, miss_txt(m)),
                      Value = sprintf("%d (%s)", tab[single], f1(100 * tab[single] / tot)), type = "row"))
  }
  if (is.null(levels)) levels <- names(tab)
  rbind(data.frame(Characteristic = paste0(label, miss_txt(m)), Value = "", type = "variable"),
        data.frame(Characteristic = paste0("    ", levels),
                   Value = sprintf("%d (%s)", as.integer(tab[levels]), f1(100 * tab[levels] / tot)),
                   type = "row"))
}
desc_num <- function(var, label, k = 0) {
  x <- a[[var]]; m <- sum(is.na(x)); x <- x[!is.na(x)]
  fk <- function(z) formatC(z, format = "f", digits = k)
  data.frame(Characteristic = c(paste0(label, ", mean ± SD", miss_txt(m)),
                                paste0(label, ", median [IQR] (range)")),
             Value = c(sprintf("%s ± %s", f1(mean(x)), f1(sd(x))),
                       sprintf("%s [%s–%s] (%s–%s)", fk(median(x)), fk(quantile(x, .25)),
                               fk(quantile(x, .75)), fk(min(x)), fk(max(x)))),
             type = "row")
}
section <- function(t) data.frame(Characteristic = t, Value = "", type = "section")

table2 <- rbind(
  section("Sociodemographic characteristics"),
  desc_num("age", "Age (years)"),
  desc_cat("age_group", "Age group"),
  desc_cat("sex", "Sex"),
  desc_cat("smoking", "Current or former smoker", single = "Yes"),
  section("Occupational characteristics"),
  desc_cat("occupation_grp", "Occupation"),
  desc_cat("department_grp", "Department", levels = c("Intensive care", "Laboratory",
           "Medicine / outpatient", "Emergency / EMS", "Pulmonology", "Administration / security")),
  desc_num("experience", "Professional experience (years)", k = 1),
  desc_cat("exp_group", "Professional experience, groups"),
  section("Exposure and high-risk practices"),
  desc_cat("tb_exposure", "Self-reported exposure to tuberculosis", single = "Yes"),
  desc_cat("agp", "Aerosol-generating procedures", single = "Yes"),
  desc_cat("specimen", "Specimen handling", single = "Yes"),
  section("Preventive measures"),
  desc_cat("mask4", "Mask use"),
  desc_cat("biosafety_training", "Biosafety training", single = "Yes"),
  desc_cat("bcg", "BCG vaccinated", single = "Yes"),
  desc_cat("bcg_scar", "BCG scar present", single = "Present"),
  section("Medical history"),
  desc_cat("tb_history", "History of tuberculosis", single = "Yes"),
  desc_cat("covid_history", "History of COVID-19", single = "Yes"),
  section("QuantiFERON result"),
  desc_cat("QFT", "QuantiFERON result", levels = c("QFT+", "QFT-"))
)
rownames(table2) <- NULL


## 5. Q1 – PREVALENCE OF LTBI ==================================================
ci_wilson <- prop.test(n_pos, n_tot, correct = FALSE)$conf.int * 100
ci_exact  <- binom.test(n_pos, n_tot)$conf.int * 100
cat(sprintf("\nLTBI prevalence: %d/%d = %.1f%% (Wilson 95%% CI %.1f–%.1f; exact %.1f–%.1f)\n",
            n_pos, n_tot, 100 * n_pos / n_tot, ci_wilson[1], ci_wilson[2], ci_exact[1], ci_exact[2]))
cat(sprintf("Half-width of the Wilson CI: %.1f percentage points\n", diff(ci_wilson) / 2))

wilson <- function(x, n) suppressWarnings(prop.test(x, n, correct = FALSE)$conf.int)
prev_by <- function(var) a %>% filter(!is.na(.data[[var]])) %>% group_by(group = .data[[var]]) %>%
  summarise(n = n(), positive = sum(qft_pos), .groups = "drop") %>% rowwise() %>%
  mutate(prevalence = 100 * positive / n, ci_low = 100 * wilson(positive, n)[1],
         ci_high = 100 * wilson(positive, n)[2]) %>% ungroup()
prev_department <- prev_by("department_grp"); prev_occupation <- prev_by("occupation_grp")
print(prev_department); print(prev_occupation)


## 6. Q2 – BIVARIATE ANALYSIS – TABLE 3 =======================================
# Tests: chi-square if all expected counts >= 5, otherwise Fisher's exact test
# (Monte Carlo p-value for tables with > 2 rows); Mann-Whitney for continuous
# variables. Crude OR: univariable logistic regression, profile-likelihood CI.
crude_or <- function(var) {
  m  <- glm(reformulate(var, "qft_pos"), data = a, family = binomial)
  ci <- suppressMessages(exp(confint(m)))
  list(or = exp(coef(m))[-1], lo = ci[-1, 1], hi = ci[-1, 2])
}
biv_num <- function(var, label) {
  x <- a[[var]]; g <- a$QFT; m <- sum(is.na(x))
  med <- function(z) sprintf("%s [%s–%s]", f1(median(z, na.rm = TRUE)),
                             f1(quantile(z, .25, na.rm = TRUE)), f1(quantile(z, .75, na.rm = TRUE)))
  o <- crude_or(var)
  data.frame(Variable = paste0(label, miss_txt(m)), Category = "median [IQR]; OR per year",
             neg = med(x[g == "QFT-"]), pos = med(x[g == "QFT+"]), OR = fOR(o$or, o$lo, o$hi),
             p = fp(wilcox.test(x ~ g, exact = FALSE)$p.value), Test = "Mann–Whitney")
}
biv_cat <- function(var, label) {
  x <- a[[var]]; m <- sum(is.na(x)); tab <- table(x, a$QFT)
  expected <- suppressWarnings(chisq.test(tab)$expected)
  fisher <- any(expected < 5)
  p <- if (fisher) fisher.test(tab, simulate.p.value = nrow(tab) > 2, B = 1e5)$p.value
       else chisq.test(tab, correct = FALSE)$p.value
  o <- crude_or(var); k <- nrow(tab)
  data.frame(Variable = c(paste0(label, miss_txt(m)), rep("", k - 1)), Category = rownames(tab),
             neg = sprintf("%d (%s)", tab[, "QFT-"], f1(100 * tab[, "QFT-"] / sum(tab[, "QFT-"]))),
             pos = sprintf("%d (%s)", tab[, "QFT+"], f1(100 * tab[, "QFT+"] / sum(tab[, "QFT+"]))),
             OR = c("1 (ref.)", fOR(o$or, o$lo, o$hi)), p = c(fp(p), rep("", k - 1)),
             Test = c(ifelse(fisher, "Fisher", "Chi²"), rep("", k - 1)))
}

set.seed(SEED)   # Monte Carlo p-values (occupation, department)
table3 <- rbind(
  biv_num("age", "Age (years)"),
  biv_cat("age_group", "Age group"),
  biv_cat("sex", "Sex"),
  biv_cat("occupation_grp", "Occupation"),
  biv_cat("department_grp", "Department"),
  biv_num("experience", "Professional experience (years)"),
  biv_cat("tb_exposure", "Self-reported TB exposure"),
  biv_cat("agp", "Aerosol-generating procedures"),
  biv_cat("specimen", "Specimen handling"),
  biv_cat("mask4", "Mask use"),
  biv_cat("biosafety_training", "Biosafety training"),
  biv_cat("bcg", "BCG vaccination"),
  biv_cat("bcg_scar", "BCG scar"),
  biv_cat("smoking", "Current or former smoker"),
  biv_cat("tb_history", "History of tuberculosis"),
  biv_cat("covid_history", "History of COVID-19")
)
rownames(table3) <- NULL
names(table3) <- c("Variable", "Category", sprintf("QFT− (n = %d)", n_neg),
                   sprintf("QFT+ (n = %d)", n_pos), "Crude OR (95% CI)", "p", "Test")

# Cochran-Armitage test for trend across the ordered mask-use categories
tm <- table(a$mask4, a$QFT)
p_trend_mask <- prop.trend.test(tm[, "QFT+"], rowSums(tm))$p.value
cat("Cochran–Armitage test for trend (mask use): p =", round(p_trend_mask, 3), "\n")


## 7. Q3 – MULTIVARIABLE LOGISTIC REGRESSION ===================================
# 7.1 Variables chosen A PRIORI on clinical and epidemiological grounds
#     (not on p-values): age (cumulative exposure), sex, self-reported TB
#     exposure, aerosol-generating procedures, regular mask use, BCG vaccination.
#     QFT is not affected by BCG (ESAT-6/CFP-10 antigens are absent from BCG).
#     Experience: strongly related to age and 16.5% missing -> not included.
#     Occupation/department: 9 extra parameters -> not enough events.
model_vars <- c("age", "sex", "tb_exposure", "agp", "mask_regular", "bcg")
f_model <- reformulate(model_vars, "qft_pos")
n_param <- 6
cat("\nEvents per parameter (all participants):", round(n_pos / n_param, 1), "\n")

# 7.2 Complete-case model
cc <- a %>% select(qft_pos, all_of(model_vars)) %>% na.omit()
cat("Complete cases: n =", nrow(cc), "| events =", sum(cc$qft_pos),
    "| EPV =", round(sum(cc$qft_pos) / n_param, 1), "\n")
m_cc <- glm(f_model, data = cc, family = binomial)
res_cc <- broom::tidy(m_cc, conf.int = TRUE, exponentiate = TRUE, conf.method = "profile") %>%
  filter(term != "(Intercept)") %>%
  transmute(method = "Complete cases", term, OR = estimate, lo = conf.low, hi = conf.high, p = p.value)

# 7.3 Diagnostics
vif_values <- car::vif(m_cc); print(vif_values)
dd <- datadist(cc); options(datadist = "dd")
m_lin <- lrm(f_model, data = cc, x = TRUE, y = TRUE)
m_rcs <- lrm(qft_pos ~ rcs(age, 3) + sex + tb_exposure + agp + mask_regular + bcg, data = cc)
p_nonlin_age <- 1 - pchisq(m_rcs$stats["Model L.R."] - m_lin$stats["Model L.R."], 1)
n_influential <- sum(cooks.distance(m_cc) > 4 / nrow(cc))
cat("Non-linearity of age: p =", round(p_nonlin_age, 3),
    "| observations with Cook's distance > 4/n:", n_influential, "\n")

# 7.4 Performance: discrimination, calibration, internal validation
roc_cc <- roc(cc$qft_pos, fitted(m_cc), quiet = TRUE, ci = TRUE)       # DeLong CI
hosmer_lemeshow <- function(y, p, g = 10) {
  grp <- cut(p, unique(quantile(p, seq(0, 1, length.out = g + 1))), include.lowest = TRUE)
  obs <- tapply(y, grp, sum); expd <- tapply(p, grp, sum); n <- tapply(y, grp, length)
  stat <- sum((obs - expd)^2 / (expd * (1 - expd / n)))
  c(chi2 = stat, df = nlevels(grp) - 2, p = 1 - pchisq(stat, nlevels(grp) - 2))
}
hl  <- hosmer_lemeshow(cc$qft_pos, fitted(m_cc))
gof <- residuals(m_lin, type = "gof")                  # le Cessie–van Houwelingen
brier      <- mean((fitted(m_cc) - cc$qft_pos)^2)
brier_null <- mean(cc$qft_pos) * (1 - mean(cc$qft_pos))
val <- validate(m_lin, method = "boot", B = 500)       # bootstrap, 500 resamples
auc_corrected   <- 0.5 + val["Dxy", "index.corrected"] / 2
slope_corrected <- val["Slope", "index.corrected"]
cal <- calibrate(m_lin, method = "boot", B = 500)

# 7.5 Firth bias-reduced logistic regression (sensitivity analysis)
m_firth <- glm(f_model, data = cc, family = binomial, method = "brglmFit", type = "AS_mean")
res_firth <- suppressWarnings(broom::tidy(m_firth)) %>% filter(term != "(Intercept)") %>%
  transmute(method = "Firth penalised", term, OR = exp(estimate),
            lo = exp(estimate - 1.96 * std.error), hi = exp(estimate + 1.96 * std.error), p = p.value)

# 7.6 Multiple imputation by chained equations (PRIMARY analysis; MAR assumption)
imp_vars <- c("qft_pos", model_vars, "experience", "occupation_grp", "department_grp",
              "biosafety_training", "smoking")
imp <- mice(a %>% select(all_of(imp_vars)), m = 30, maxit = 20, printFlag = FALSE, seed = SEED)
fit_mi <- with(imp, glm(qft_pos ~ age + sex + tb_exposure + agp + mask_regular + bcg, family = binomial))
res_mi <- as.data.frame(summary(pool(fit_mi), conf.int = TRUE)) %>% filter(term != "(Intercept)") %>%
  transmute(method = "Multiple imputation", term = as.character(term), OR = exp(estimate),
            lo = exp(`2.5 %`), hi = exp(`97.5 %`), p = p.value)
res_multi <- bind_rows(res_mi, res_cc, res_firth)

# 7.7 Sensitivity analysis: adding department (likelihood ratio test, same participants)
d_dep <- a %>% select(qft_pos, all_of(model_vars), department_grp) %>% na.omit()
m_without <- glm(f_model, data = d_dep, family = binomial)
m_with    <- update(m_without, . ~ . + department_grp)
p_department <- anova(m_without, m_with, test = "LRT")[2, "Pr(>Chi)"]
cat("Adding department to the model (LRT, 5 df): p =", round(p_department, 3), "\n")

# Table 4 – adjusted odds ratios by estimation method
term_labels <- c(age = "Age (per year)", sexMale = "Male vs female",
  tb_exposureYes = "Self-reported TB exposure (yes vs no)",
  agpYes = "Aerosol-generating procedures (yes vs no)",
  mask_regularIrregular = "Irregular vs regular mask use", bcgNo = "Not BCG vaccinated vs vaccinated")
cell <- function(df, t) with(df[df$term == t, ], sprintf("%s; %s", fOR(OR, lo, hi), fp(p)))
table4 <- data.frame(
  Variable = unname(term_labels),
  `Multiple imputation (primary), n = 164` = sapply(names(term_labels), cell, df = res_mi),
  `Complete cases, n = 123` = sapply(names(term_labels), cell, df = res_cc),
  `Firth penalised, n = 123` = sapply(names(term_labels), cell, df = res_firth),
  check.names = FALSE)
names(table4)[2:4] <- c(sprintf("Multiple imputation (primary), n = %d", n_tot),
                        sprintf("Complete cases, n = %d", nrow(cc)),
                        sprintf("Firth penalised, n = %d", nrow(cc)))
rownames(table4) <- NULL

# Table 5 – model diagnostics and performance
table5 <- data.frame(
  Indicator = c("Events per parameter", "Maximum variance inflation factor",
                "Non-linearity of age (restricted cubic spline, 3 knots)",
                "Influential observations (Cook's distance > 4/n)", "Apparent AUC (95% CI, DeLong)",
                "Optimism-corrected AUC (500 bootstrap resamples)", "Optimism-corrected calibration slope",
                "Hosmer–Lemeshow test (10 groups)", "le Cessie–van Houwelingen test", "Brier score"),
  Value = c(sprintf("%s (%s with imputation)", f1(sum(cc$qft_pos) / n_param), f1(n_pos / n_param)),
            f2(max(vif_values)), sprintf("p = %s", f2(p_nonlin_age)), as.character(n_influential),
            sprintf("%s (%s–%s)", f2(roc_cc$auc), f2(roc_cc$ci[1]), f2(roc_cc$ci[3])),
            f2(auc_corrected), f2(slope_corrected), sprintf("p = %s", f2(hl["p"])),
            sprintf("p = %s", f2(gof["P"])),
            sprintf("%s (null model %s)", formatC(brier, format = "f", digits = 3),
                    formatC(brier_null, format = "f", digits = 3))),
  Interpretation = c("Acceptable for six parameters", "No multicollinearity", "Linear effect of age adequate",
                     "Expected number; results stable across methods", "Poor discrimination",
                     "Poor discrimination", "Moderate overfitting", "No significant lack of fit",
                     "No significant lack of fit", "Small improvement over no model"))
print(table4); print(table5)


## 8. FIGURES ==================================================================
col_qft <- c("QFT-" = "#4C72B0", "QFT+" = "#C44E52")
th <- theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), panel.grid.minor = element_blank())
save_fig <- function(g, file, w, h, dir = fig_dir)
  ggsave(file.path(dir, file), g, width = w, height = h, dpi = 300, bg = "white")

# Figure 1 – flow diagram (STROBE item 13c)
n_raw <- nrow(d); n_indet <- sum(is.na(d$qft_pos))
boxes <- data.frame(
  x = c(4, 8.3, 4, 2.2, 5.8, 4), y = c(9.1, 7.4, 5.85, 3.55, 3.55, 1.35),
  w = c(5.6, 3.0, 5.6, 3.2, 3.2, 6.8), h = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.3),
  fill = c("white", "grey95", "white", "#F8DCDC", "#DCE6F2", "white"),
  label = c(sprintf("Healthcare workers enrolled and tested\nwith QuantiFERON (n = %d)", n_raw),
            sprintf("Excluded:\nindeterminate QFT result\n(n = %d)", n_indet),
            sprintf("Participants analysed (n = %d)\nbivariate analysis and multiple imputation", n_tot),
            sprintf("QFT positive\n(LTBI) n = %d", n_pos), sprintf("QFT negative\nn = %d", n_neg),
            sprintf(paste0("Complete-case sensitivity analyses (n = %d; %d QFT+)\n",
                           "%d participants with ≥ 1 missing value among model variables\n",
                           "(mask use %d, BCG vaccination %d, TB exposure %d; overlapping)"),
                    nrow(cc), sum(cc$qft_pos), n_tot - nrow(cc), sum(is.na(a$mask_regular)),
                    sum(is.na(a$bcg)), sum(is.na(a$tb_exposure)))),
  bold = c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))
arrows_df <- data.frame(x = c(4, 6.0, 4, 4), y = c(8.6, 7.4, 5.35, 5.35),
                        xend = c(4, 6.8, 2.2, 5.8), yend = c(6.35, 7.4, 4.05, 4.05))
fig1 <- ggplot() +
  geom_rect(data = boxes, aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
                              fill = fill), colour = "grey20") +
  scale_fill_identity() +
  geom_text(data = boxes, aes(x, y, label = label, fontface = ifelse(bold, "bold", "plain")), size = 3.4) +
  geom_segment(aes(x = 4, y = 7.4, xend = 6.0, yend = 7.4), colour = "grey20") +
  geom_segment(data = arrows_df, aes(x, y, xend = xend, yend = yend), colour = "grey20",
               arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
  geom_segment(aes(x = 4, y = 3.05, xend = 4, yend = 2.0), colour = "grey20", linetype = 2) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0.5, 9.8)) + theme_void()
save_fig(fig1, "Figure1_flow_diagram.png", 7.5, 6.2)

# Figure 2 – descriptive overview by QFT result (panels A–F)
pA <- a %>% count(res = QFT) %>% mutate(pct = 100 * n / sum(n)) %>%
  ggplot(aes(res, n, fill = res)) + geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", n, pct)), vjust = -0.2, size = 3.4) +
  scale_fill_manual(values = col_qft, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = NULL, y = "Participants", title = "A. QuantiFERON result",
       subtitle = sprintf("LTBI prevalence: %.1f%% (95%% CI %.1f–%.1f)",
                          100 * n_pos / n_tot, ci_wilson[1], ci_wilson[2])) + th
pB <- ggplot(a, aes(age, fill = QFT)) +
  geom_histogram(binwidth = 2, boundary = 18, colour = "white") +
  geom_vline(xintercept = median(a$age), linetype = 2) + scale_fill_manual(values = col_qft) +
  labs(x = "Age (years)", y = "Participants", fill = NULL, title = "B. Age distribution",
       subtitle = sprintf("Median %.0f years [IQR %.0f–%.0f]; %.0f%% women", median(a$age),
                          quantile(a$age, .25), quantile(a$age, .75), 100 * mean(a$sex == "Female"))) +
  th + theme(legend.position = "top")
bars_by <- function(var, title) {
  df <- a %>% filter(!is.na(.data[[var]])) %>% count(grp = .data[[var]], QFT)
  tot <- df %>% group_by(grp) %>% summarise(pos = sum(n[QFT == "QFT+"]), n = sum(n), .groups = "drop") %>%
    mutate(lab = sprintf("n=%d · %.0f%% QFT+", n, 100 * pos / n))
  ggplot(df, aes(y = reorder(grp, n, sum), x = n, fill = QFT)) + geom_col(width = 0.7) +
    geom_text(data = tot, inherit.aes = FALSE, aes(y = grp, x = n, label = lab), hjust = -0.05, size = 2.8) +
    scale_fill_manual(values = col_qft, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.55))) +
    labs(x = "Participants", y = NULL, title = title) + th
}
pC <- bars_by("occupation_grp", "C. Occupation")
pD <- bars_by("department_grp", "D. Department")
binary_vars <- c(tb_exposure = "Self-reported TB exposure", agp = "Aerosol-generating procedures",
                 specimen = "Specimen handling", bcg = "BCG vaccinated", bcg_scar = "BCG scar present",
                 biosafety_training = "Biosafety training", mask_regular = "Regular mask use",
                 smoking = "Current or former smoker", covid_history = "History of COVID-19")
yes_level <- c(tb_exposure = "Yes", agp = "Yes", specimen = "Yes", bcg = "Yes", bcg_scar = "Present",
               biosafety_training = "Yes", mask_regular = "Regular", smoking = "Yes", covid_history = "Yes")
pE_df <- lapply(names(binary_vars), function(v) {
  a %>% filter(!is.na(.data[[v]])) %>% group_by(QFT) %>%
    summarise(pct = 100 * mean(.data[[v]] == yes_level[v]), .groups = "drop") %>%
    mutate(variable = binary_vars[v], miss = sum(is.na(a[[v]])))
}) %>% bind_rows() %>%
  mutate(variable = ifelse(miss > 0, sprintf("%s (%d missing)", variable, miss), variable))
order_E <- pE_df %>% group_by(variable) %>% summarise(m = mean(pct)) %>% arrange(m) %>% pull(variable)
pE <- pE_df %>% mutate(variable = factor(variable, order_E), QFT = factor(QFT, c("QFT+", "QFT-"))) %>%
  ggplot(aes(x = pct, y = variable, fill = QFT)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_dodge(width = 0.8), hjust = -0.15, size = 3) +
  scale_fill_manual(values = col_qft, breaks = c("QFT-", "QFT+"),
                    labels = c(sprintf("QFT- participants (n = %d)", n_neg), sprintf("QFT+ participants (n = %d)", n_pos))) +
  scale_x_continuous(limits = c(0, 108), breaks = seq(0, 100, 25), labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Participants answering 'yes' (%)", y = NULL, fill = NULL,
       title = "E. Exposure, prevention and medical history by QFT result",
       subtitle = "Each pair of bars compares QFT- and QFT+ participants; missing answers excluded from percentages") +
  th + theme(legend.position = "top", panel.grid.major.y = element_blank())
pF_df <- a %>% filter(!is.na(mask4)) %>% count(mask4, QFT) %>% group_by(mask4) %>%
  mutate(pct = 100 * n / sum(n), tot = sum(n)) %>% ungroup()
pF <- ggplot(pF_df, aes(mask4, pct, fill = QFT)) + geom_col(width = 0.7) +
  geom_text(data = filter(pF_df, QFT == "QFT+"), aes(y = pct / 2, label = sprintf("%.0f%%", pct)),
            colour = "white", size = 3) +
  geom_text(data = distinct(pF_df, mask4, tot), inherit.aes = FALSE,
            aes(mask4, 104, label = paste0("n=", tot)), size = 2.8) +
  scale_fill_manual(values = col_qft, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = c(0, 50, 100)) +
  labs(x = NULL, y = "Participants (%)", title = "F. Mask use",
       subtitle = sprintf("%d missing; white figure = %% QFT+", sum(is.na(a$mask4)))) + th
fig2 <- wrap_elements(pA | pB | pF) / wrap_elements(pC | pD) / wrap_elements(pE) +
  plot_layout(heights = c(1, 1.1, 1.5)) +
  plot_annotation(caption = sprintf(paste("n = %d participants. Blue = QFT-, red = QFT+.",
    "Percentages on available data; Wilson 95%% CI."), n_tot))
save_fig(fig2, "Figure2_descriptive_overview.png", 14, 15)

# Figure 3 – prevalence by department and occupation
prev_plot <- function(df, title) {
  ggplot(df, aes(reorder(group, prevalence), prevalence)) +
    geom_hline(yintercept = 100 * n_pos / n_tot, linetype = 2, colour = "grey40") +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high), colour = "#C44E52") +
    geom_text(aes(label = paste0(positive, "/", n)), nudge_x = 0.3, size = 3) +
    coord_flip() + scale_y_continuous(limits = c(0, 100)) +
    labs(x = NULL, y = "QFT+ prevalence (%) with 95% CI", title = title) + th
}
fig3 <- (prev_plot(prev_department, "A. By department") | prev_plot(prev_occupation, "B. By occupation")) +
  plot_annotation(caption = sprintf("Dashed line: overall prevalence (%.1f%%). Wilson 95%% confidence intervals.",
                                    100 * n_pos / n_tot))
save_fig(fig3, "Figure3_prevalence_department_occupation.png", 12, 4.5)

# Figure 4 – forest plot of adjusted odds ratios (three estimation methods)
method_labels <- c("Multiple imputation" = sprintf("Multiple imputation (primary, n = %d)", n_tot),
                   "Complete cases" = sprintf("Complete cases (n = %d)", nrow(cc)),
                   "Firth penalised" = sprintf("Firth penalised (n = %d)", nrow(cc)))
forest_labels <- c(age = "Age (per year)", sexMale = "Male vs female", tb_exposureYes = "Self-reported TB exposure",
                   agpYes = "Aerosol-generating procedures", mask_regularIrregular = "Irregular vs regular mask use",
                   bcgNo = "Not BCG vaccinated")
fig4 <- res_multi %>%
  mutate(term = factor(forest_labels[term], rev(forest_labels)),
         method = factor(method_labels[method], method_labels)) %>%
  ggplot(aes(OR, term, colour = method)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6)) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8)) +
  scale_colour_manual(values = c("#C44E52", "#4C72B0", "#55A868")) +
  labs(x = "Adjusted odds ratio (log scale) with 95% CI", y = NULL, colour = NULL) +
  th + theme(legend.position = "bottom", legend.direction = "vertical")
save_fig(fig4, "Figure4_forest_plot.png", 8, 5.5)

# Figure 5 – model performance: ROC curve and bootstrap calibration
p5a <- pROC::ggroc(roc_cc, legacy.axes = TRUE, colour = "#C44E52", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  annotate("text", x = 0.62, y = 0.15, size = 3.4,
           label = sprintf("AUC %.2f (95%% CI %.2f–%.2f)\nOptimism-corrected AUC %.2f",
                           roc_cc$auc, roc_cc$ci[1], roc_cc$ci[3], auc_corrected)) +
  coord_equal() + labs(x = "1 – specificity", y = "Sensitivity", title = "A. ROC curve") + th
cal_df <- data.frame(pred = cal[, "predy"], apparent = cal[, "calibrated.orig"],
                     corrected = cal[, "calibrated.corrected"]) %>%
  pivot_longer(-pred, names_to = "type", values_to = "obs") %>%
  mutate(type = factor(type, c("apparent", "corrected"), c("Apparent", "Bias-corrected (bootstrap)")))
p5b <- ggplot(cal_df, aes(pred, obs, linetype = type)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50") +
  geom_line(colour = "#4C72B0", linewidth = 0.9) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  annotate("text", x = 0.68, y = 0.1, size = 3.4,
           label = sprintf("Corrected calibration slope %.2f", slope_corrected)) +
  labs(x = "Predicted probability", y = "Observed proportion", linetype = NULL, title = "B. Calibration") +
  th + theme(legend.position = "bottom")
fig5 <- (p5a | p5b) + plot_annotation(caption = sprintf(paste("Complete-case model (n = %d); 500 bootstrap",
  "resamples. Grey diagonal = chance (A) or perfect calibration (B)."), nrow(cc)))
save_fig(fig5, "Figure5_model_performance.png", 11, 5.8)

# Supplementary figures: distributions, Q-Q plots, missing-data map
long_num <- a %>% select(QFT, age, experience) %>% pivot_longer(-QFT, names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>% mutate(variable = labels_en[variable])
sA <- ggplot(long_num, aes(value)) + geom_histogram(bins = 20, fill = "grey70", colour = "white") +
  facet_wrap(~ variable, scales = "free") + labs(x = NULL, y = "Participants", title = "Distributions") + th
sB <- ggplot(long_num, aes(sample = value)) + stat_qq() + stat_qq_line(colour = "#C44E52") +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "Theoretical quantiles", y = "Observed quantiles", title = "Normal Q-Q plots") + th
save_fig(sA / sB, "S_distributions_qqplots.png", 8, 7, supp_dir)
miss_map <- a %>% mutate(row = row_number()) %>% select(row, all_of(vars_missing)) %>%
  mutate(across(-row, is.na)) %>% pivot_longer(-row, names_to = "variable", values_to = "missing") %>%
  mutate(variable = factor(labels_en[variable], rev(labels_en[vars_missing])))
sC <- ggplot(miss_map, aes(row, variable, fill = missing)) + geom_tile() +
  scale_fill_manual(values = c(`FALSE` = "grey90", `TRUE` = "#C44E52"), labels = c("Recorded", "Missing")) +
  labs(x = "Participant", y = NULL, fill = NULL, title = "Missing-data map") +
  th + theme(legend.position = "top", panel.grid = element_blank())
save_fig(sC, "S_missing_data_map.png", 9, 4.5, supp_dir)


## 9. EXPORT OF TABLES (WORD AND EXCEL) ========================================
ft_basic <- function(df, caption, note, widths, landscape_cols = NULL) {
  ft <- flextable(df) %>% bold(part = "header") %>% set_caption(caption) %>%
    add_footer_lines(note) %>% fontsize(size = 9, part = "all") %>% fontsize(size = 8, part = "footer") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    padding(padding.top = 1, padding.bottom = 1, part = "all") %>%
    align(j = seq_len(ncol(df))[-1], align = "center", part = "all") %>% border_outer()
  for (j in seq_along(widths)) ft <- width(ft, j = j, width = widths[j])
  ft
}
t2 <- table2
ft2 <- ft_basic(t2[, 1:2] %>% rename(`n (%) or value` = Value),
  sprintf("Table 2. Characteristics of the study population (n = %d)", n_tot),
  paste("Values are n (%) unless otherwise stated; percentages are calculated on available data.",
        "IQR, interquartile range; SD, standard deviation; EMT, emergency medical technician;",
        "EMS, emergency medical services."), c(4.3, 2.0)) %>%
  bold(i = which(t2$type == "section")) %>% bg(i = which(t2$type == "section"), bg = "#EEF2F7") %>%
  italic(i = which(t2$type == "variable"))
ft3 <- ft_basic(table3,
  sprintf("Table 3. Characteristics of participants by QuantiFERON result, with crude odds ratios (n = %d)", n_tot),
  paste("Values are n (column %) calculated on available data, or median [IQR]. p-values: Pearson chi-square",
        "test when all expected counts were ≥ 5, otherwise Fisher's exact test (Monte Carlo p-value with 100,000",
        "replicates for variables with more than two categories); Mann–Whitney U test for continuous variables.",
        "For variables with more than two categories, the p-value is a global test. Crude OR: univariable logistic",
        "regression with profile-likelihood 95% CI; ref., reference category. Because p-values and confidence",
        "intervals are derived from different methods, they may disagree for estimates close to the 0.05 threshold."),
  c(2.2, 1.8, 1.1, 1.1, 1.4, 0.6, 1.0)) %>% bold(j = 1, part = "body")
ft4 <- ft_basic(table4, "Table 4. Factors associated with QFT positivity: multivariable logistic regression",
  paste("Values are adjusted odds ratio (95% CI); p-value. All variables are mutually adjusted. Multiple",
        "imputation by chained equations (30 imputed datasets, 20 iterations), pooled with Rubin's rules.",
        "Complete-case estimates use profile-likelihood confidence intervals; Firth estimates use bias-reduced",
        "logistic regression."), c(2.3, 1.45, 1.45, 1.45))
ft5 <- ft_basic(table5, sprintf("Table 5. Diagnostics and performance of the multivariable model (complete cases, n = %d; %d events)",
                                nrow(cc), sum(cc$qft_pos)),
  "AUC, area under the receiver operating characteristic curve.", c(2.8, 1.7, 2.1))
ftS1 <- ft_basic(table_S1, sprintf("Supplementary Table S1. Missing data by variable (n = %d)", n_tot),
  "Medical history and high-risk practices: a blank answer was coded as \"none\".", c(3.5, 1.2, 1.2))

landscape <- block_section(prop_section(page_size = page_size(orient = "landscape"), type = "nextPage"))
portrait  <- block_section(prop_section(page_size = page_size(orient = "portrait"),  type = "nextPage"))
doc <- read_docx() %>%
  body_add_flextable(ft2) %>% body_end_block_section(portrait) %>%
  body_add_flextable(ft3) %>% body_end_block_section(landscape) %>%
  body_add_flextable(ft4) %>% body_add_par("") %>% body_add_flextable(ft5) %>%
  body_add_par("") %>% body_add_flextable(ftS1) %>% body_end_block_section(portrait)
print(doc, target = file.path(tab_dir, "Thesis_tables.docx"))

openxlsx::write.xlsx(list("Table 2" = table2[, 1:2], "Table 3" = table3, "Table 4" = table4,
                          "Table 5" = table5, "Table S1" = table_S1,
                          "Prevalence by department" = prev_department,
                          "Prevalence by occupation" = prev_occupation,
                          "Adjusted OR (all methods)" = res_multi),
                     file.path(tab_dir, "Thesis_tables.xlsx"))


## 11. MISSING DATA ANALYSIS ===================================================
# Describes how much data are missing, for which variables and participants,
# and whether missingness is related to the QFT result or to other
# characteristics (information supporting the missing-at-random assumption).
all_vars <- c("age", "sex", "occupation_grp", "department_grp", "experience", "tb_exposure",
              "agp", "specimen", "mask4", "biosafety_training", "bcg", "bcg_scar", "smoking",
              "tb_history", "covid_history")

# 11.1 Table S1 – missing values by variable, overall and by QFT result
#      p-value: Fisher's exact test comparing the proportion of missing values
#      between QFT-negative and QFT-positive participants.
npct <- function(m, n) sprintf("%d (%s)", m, f1(100 * m / n))
table_S1_full <- lapply(all_vars, function(v) {
  miss <- is.na(a[[v]])
  m_neg <- sum(miss[a$QFT == "QFT-"]); m_pos <- sum(miss[a$QFT == "QFT+"])
  p <- if (sum(miss) > 0) fisher.test(table(miss, a$QFT))$p.value else NA
  data.frame(Variable = unname(labels_en[v]),
             `Missing, n (%)` = npct(sum(miss), n_tot),
             `Missing in QFT− (n = 100), n (%)` = npct(m_neg, n_neg),
             `Missing in QFT+ (n = 64), n (%)` = npct(m_pos, n_pos),
             p = fp(p), n_missing = sum(miss), check.names = FALSE)
}) %>% bind_rows() %>% arrange(desc(n_missing))
names(table_S1_full)[3:4] <- c(sprintf("Missing in QFT− (n = %d), n (%%)", n_neg),
                               sprintf("Missing in QFT+ (n = %d), n (%%)", n_pos))
table_S1_full$p[table_S1_full$n_missing == 0] <- "—"
table_S1_full <- table_S1_full %>% select(-n_missing)
print(table_S1_full, row.names = FALSE)

# 11.2 Number of missing values per participant
miss_vars <- all_vars[sapply(a[all_vars], function(x) any(is.na(x)))]
n_miss_person <- rowSums(is.na(a[miss_vars]))
table_S2 <- as.data.frame(table(n_miss_person)) %>%
  transmute(`Number of missing values` = as.character(n_miss_person), `Participants, n` = Freq,
            `Participants, %` = f1(100 * Freq / n_tot))
cat("\nParticipants with complete data on all variables:", sum(n_miss_person == 0),
    sprintf("(%.1f%%)", 100 * mean(n_miss_person == 0)), "\n")
cat("Total missing values:", sum(is.na(a[miss_vars])), "of", n_tot * length(all_vars),
    sprintf("(%.1f%% of all cells)", 100 * sum(is.na(a[miss_vars])) / (n_tot * length(all_vars))), "\n")

# 11.3 Table S3 – missing-data patterns (which variables are missing together)
short_names <- c(bcg_scar = "BCG scar", biosafety_training = "Biosafety training", mask4 = "Mask use",
                 experience = "Experience", bcg = "BCG vaccination", tb_exposure = "TB exposure",
                 smoking = "Smoking", occupation_grp = "Occupation", department_grp = "Department")
pattern <- apply(is.na(a[miss_vars]), 1, function(r) {
  if (!any(r)) "None (complete data)" else paste(short_names[miss_vars[r]], collapse = ", ")
})
table_S3 <- data.frame(pattern, n_miss = n_miss_person) %>%
  count(`Missing variables` = pattern, `Number missing` = n_miss, name = "Participants, n") %>%
  arrange(desc(`Participants, n`)) %>%
  mutate(`Participants, %` = f1(100 * `Participants, n` / n_tot))
print(table_S3, row.names = FALSE)

# 11.4 Table S4 – participants with complete vs incomplete data for the model
#      variables (complete-case analysis vs excluded participants)
a$complete_model <- factor(ifelse(complete.cases(a[model_vars]), "Complete", "Incomplete"),
                           levels = c("Complete", "Incomplete"))
cmp_cat <- function(var, label, level) {
  tab <- table(a[[var]], a$complete_model); col <- colSums(tab)
  p <- fisher.test(tab, simulate.p.value = nrow(tab) > 2, B = 1e5)$p.value
  data.frame(Characteristic = label,
             Complete = sprintf("%d (%s)", tab[level, "Complete"], f1(100 * tab[level, "Complete"] / col["Complete"])),
             Incomplete = sprintf("%d (%s)", tab[level, "Incomplete"], f1(100 * tab[level, "Incomplete"] / col["Incomplete"])),
             p = fp(p), Test = "Fisher")
}
cmp_num <- function(var, label) {
  x <- a[[var]]; g <- a$complete_model
  med <- function(z) sprintf("%s [%s–%s]", f1(median(z, na.rm = TRUE)),
                             f1(quantile(z, .25, na.rm = TRUE)), f1(quantile(z, .75, na.rm = TRUE)))
  data.frame(Characteristic = paste0(label, ", median [IQR]"), Complete = med(x[g == "Complete"]),
             Incomplete = med(x[g == "Incomplete"]),
             p = fp(wilcox.test(x ~ g, exact = FALSE)$p.value), Test = "Mann–Whitney")
}
set.seed(SEED)
table_S4 <- rbind(
  cmp_cat("QFT", "QFT positive", "QFT+"),
  cmp_num("age", "Age (years)"),
  cmp_cat("sex", "Female sex", "Female"),
  cmp_cat("occupation_grp", "Nurse", "Nurse"),
  cmp_cat("department_grp", "Intensive care department", "Intensive care"),
  cmp_cat("agp", "Aerosol-generating procedures", "Yes")
)
names(table_S4)[2:3] <- c(sprintf("Complete data (n = %d)", sum(a$complete_model == "Complete")),
                          sprintf("Incomplete data (n = %d)", sum(a$complete_model == "Incomplete")))
print(table_S4, row.names = FALSE)

# 11.5 Supplementary figures: missing values by variable and by QFT result
miss_plot_df <- lapply(miss_vars, function(v) a %>% group_by(QFT) %>%
  summarise(pct = 100 * mean(is.na(.data[[v]])), .groups = "drop") %>% mutate(variable = labels_en[v])) %>%
  bind_rows()
order_miss <- miss_plot_df %>% group_by(variable) %>% summarise(m = mean(pct)) %>% arrange(m) %>% pull(variable)
fig_miss <- miss_plot_df %>% mutate(variable = factor(variable, order_miss), QFT = factor(QFT, c("QFT+", "QFT-"))) %>%
  ggplot(aes(pct, variable, fill = QFT)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), position = position_dodge(width = 0.8), hjust = -0.15, size = 3) +
  scale_fill_manual(values = col_qft, breaks = c("QFT-", "QFT+")) +
  scale_x_continuous(limits = c(0, 35), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0))) +
  labs(x = "Missing values (%)", y = NULL, fill = NULL,
       title = "Missing values by variable and QuantiFERON result") +
  th + theme(legend.position = "top", panel.grid.major.y = element_blank())
save_fig(fig_miss, "S_missing_by_variable_and_QFT.png", 8, 5, supp_dir)

fig_pattern <- table_S3 %>% filter(`Missing variables` != "None (complete data)") %>%
  mutate(`Missing variables` = reorder(str_wrap(`Missing variables`, 45), `Participants, n`)) %>%
  ggplot(aes(`Participants, n`, `Missing variables`)) + geom_col(fill = "#C44E52", width = 0.7) +
  geom_text(aes(label = `Participants, n`), hjust = -0.3, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Participants", y = NULL, title = "Missing-data patterns",
       subtitle = sprintf("%d participants (%.1f%%) had complete data on all variables",
                          sum(n_miss_person == 0), 100 * mean(n_miss_person == 0))) + th
save_fig(fig_pattern, "S_missing_data_patterns.png", 9, 7.5, supp_dir)

# 11.6 Export (Word and Excel)
note_blank <- paste("Medical history, aerosol-generating procedures and specimen handling were free-text",
                    "fields in which a blank answer was coded as \"none\"; they therefore have no missing values.")
ftS1f <- ft_basic(table_S1_full, sprintf("Supplementary Table S1. Missing values by variable, overall and by QuantiFERON result (n = %d)", n_tot),
  paste("p-value: Fisher's exact test comparing the proportion of missing values between QFT-negative and",
        "QFT-positive participants.", note_blank), c(2.3, 1.1, 1.4, 1.4, 0.6))
ftS2 <- ft_basic(table_S2, "Supplementary Table S2. Number of missing values per participant",
  sprintf("Calculated over the %d variables with at least one missing value.", length(miss_vars)), c(2.0, 1.3, 1.3))
ftS3 <- ft_basic(table_S3, "Supplementary Table S3. Missing-data patterns",
  "Each row is a combination of variables missing together in the same participants.", c(3.6, 1.0, 1.0, 1.0)) %>%
  align(j = 1, align = "left", part = "all")
ftS4 <- ft_basic(table_S4, "Supplementary Table S4. Characteristics of participants with complete and incomplete data for the multivariable model",
  paste("Complete data: no missing value among the six model variables (age, sex, TB exposure, aerosol-generating",
        "procedures, mask use, BCG vaccination). Values are n (%) or median [IQR]. Fisher's exact test",
        "(Monte Carlo p-value for variables with more than two categories) or Mann–Whitney U test.",
        "For occupation and department, the p-value compares the whole distribution across all categories."),
  c(2.4, 1.4, 1.4, 0.6, 1.0))
doc_miss <- read_docx() %>% body_add_flextable(ftS1f) %>% body_add_par("") %>%
  body_add_flextable(ftS2) %>% body_add_par("") %>% body_add_flextable(ftS3) %>% body_add_par("") %>%
  body_add_flextable(ftS4) %>% body_add_par("") %>%
  body_add_img(file.path(supp_dir, "S_missing_by_variable_and_QFT.png"), width = 6.3, height = 3.9) %>%
  body_add_par("Supplementary Figure S1. Missing values by variable and QuantiFERON result.") %>%
  body_add_img(file.path(supp_dir, "S_missing_data_patterns.png"), width = 6.3, height = 5.25) %>%
  body_add_par("Supplementary Figure S2. Missing-data patterns.") %>%
  body_add_img(file.path(supp_dir, "S_missing_data_map.png"), width = 6.3, height = 3.15) %>%
  body_add_par("Supplementary Figure S3. Missing-data map (each column is a participant; red cells are missing values).")
print(doc_miss, target = file.path(tab_dir, "Missing_data_tables.docx"))
openxlsx::write.xlsx(list("S1 Missing by variable" = table_S1_full, "S2 Missing per participant" = table_S2,
                          "S3 Missing patterns" = table_S3, "S4 Complete vs incomplete" = table_S4),
                     file.path(tab_dir, "Missing_data_tables.xlsx"))
cat("\nMissing-data analysis exported to", tab_dir, "\n")


## 12. REPRODUCIBILITY ========================================================
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("\nAnalysis complete. Tables in", tab_dir, "| figures in", fig_dir, "\n")
