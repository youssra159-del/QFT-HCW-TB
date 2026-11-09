###############################################################################
#  Infection tuberculeuse latente (ITL) chez le personnel de santé
#  Analyse statistique complète – résultat QuantiFERON (QFT)
#
#  Questions de recherche
#   Q1. Quelle est la prévalence de l'ITL (QFT+) dans cette population ?
#   Q2. Quels facteurs sont associés à un QFT+ en analyse bivariée ?
#   Q3. Quels facteurs restent associés après ajustement (régression logistique
#       multivariable), et le modèle est-il fiable ?
#
###############################################################################

## 0. Packages et paramètres ---------------------------------------------------
pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "forcats", "ggplot2",
          "broom", "car", "pROC", "rms", "brglm2", "mice", "openxlsx")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants) > 0) install.packages(manquants)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))


fichier  <- "Data_QFT_Youssra.xlsx"   # adapter le chemin si besoin
dossier  <- "resultats"               # dossier de sortie (tableaux, figures)
dir.create(dossier, showWarnings = FALSE)
set.seed(2026)                        # reproductibilité (bootstrap, imputation)
options(digits = 4)

## 1. Import -------------------------------------------------------------------
brut <- read_excel(fichier, na = c("", "NA"))   # "NA" saisi en texte = manquant
names(brut) <- c("id", "age", "sexe", "profession", "service", "exposition",
                 "experience", "masque", "formation", "bcg", "cicatrice",
                 "tabac", "qft", "antecedents", "pratiques")
table(brut$qft)
cat("Observations importées :", nrow(brut), "\n")

## 2. Nettoyage et recodage ----------------------------------------------------

txt <- function(x) str_squish(as.character(x))
oui_non <- function(x) {
  x <- str_to_lower(txt(x))
  case_when(x %in% c("yes", "oui", "rarement") ~ "Oui",   # "rarement" = exposé
            x %in% c("no", "n", "non")         ~ "Non",
            TRUE ~ NA_character_)                          # "RAS" = non interprétable
}

d <- brut %>%
  mutate(
    across(everything(), ~ if (is.character(.x)) txt(.x) else .x),
    age        = as.numeric(age),
    experience = suppressWarnings(as.numeric(experience)),
    # Issue principale : QFT+ (1) / QFT- (0) ; indéterminé = manquant
    qft_pos = case_when(qft == "QFT+" ~ 1L, qft == "QFT-" ~ 0L, TRUE ~ NA_integer_),
    sexe = factor(sexe, levels = c("Female", "Male"), labels = c("Femme", "Homme")),

    # Profession regroupée en catégories cliniquement homogènes
    profession_g = case_when(
      str_detect(profession, regex("doctor|médcin|docteur|intern", TRUE)) ~ "Médecin / interne",
      str_detect(profession, regex("^nurse", TRUE))                       ~ "Infirmier",
      str_detect(profession, regex("nursing assistant|ambulancier", TRUE)) ~ "Aide-soignant / ambulancier",
      str_detect(profession, regex("laboratory|scientifique|technical|pharmacien", TRUE)) ~ "Laboratoire / pharmacie",
      !is.na(profession) ~ "Personnel non soignant",
      TRUE ~ NA_character_),
    profession_g = factor(profession_g, levels = c("Infirmier", "Médecin / interne",
                          "Aide-soignant / ambulancier", "Laboratoire / pharmacie",
                          "Personnel non soignant")),

    # Service regroupé
    service_g = case_when(
      str_detect(service, regex("intensive care", TRUE))               ~ "Réanimation",
      str_detect(service, regex("pulmonology", TRUE))                  ~ "Pneumologie",
      str_detect(service, regex("emergency|resuscitation", TRUE))      ~ "Urgences / SAMU",
      str_detect(service, regex("laboratory|genetics", TRUE))          ~ "Laboratoire",
      str_detect(service, regex("consultation|internal|ophthalmo", TRUE)) ~ "Médecine / consultations",
      !is.na(service) ~ "Administration / sécurité",
      TRUE ~ NA_character_),
    service_g = factor(service_g, levels = c("Réanimation", "Pneumologie", "Urgences / SAMU",
                       "Laboratoire", "Médecine / consultations", "Administration / sécurité")),

    exposition = factor(oui_non(exposition), levels = c("Non", "Oui")),
    formation  = factor(oui_non(formation),  levels = c("Oui", "Non")),
    bcg        = factor(oui_non(bcg),        levels = c("Oui", "Non")),
    cicatrice  = factor(oui_non(cicatrice),  levels = c("Oui", "Non")),

    # Port du masque : variable ordinale ; "NeveRarely" = Rarement ; "P" = manquant
    masque = case_when(masque %in% c("Never")                 ~ "Jamais / rarement",
                       masque %in% c("Rarely", "NeveRarely")  ~ "Jamais / rarement",
                       masque == "Sometimes" ~ "Parfois",
                       masque == "Often"     ~ "Souvent",
                       masque == "Always"    ~ "Toujours",
                       TRUE ~ NA_character_),
    masque = factor(masque, levels = c("Toujours", "Souvent", "Parfois", "Jamais / rarement")),
    masque_reg = factor(ifelse(is.na(masque), NA,
                        ifelse(masque %in% c("Toujours", "Souvent"), "Régulier", "Irrégulier")),
                        levels = c("Régulier", "Irrégulier")),

    # Tabac : fumeur actuel ou ancien = "Oui"
    tabac = case_when(str_to_lower(tabac) == "no" ~ "Non",
                      is.na(tabac) ~ NA_character_,
                      TRUE ~ "Oui"),
    tabac = factor(tabac, levels = c("Non", "Oui")),

    # Antécédents (champ libre ; vide = aucun antécédent déclaré)
    antecedents = ifelse(is.na(antecedents) | str_to_upper(antecedents) == "NA", "", antecedents),
    atcd_tb    = factor(ifelse(str_detect(antecedents, regex("tubercul", TRUE)), "Oui", "Non"), c("Non", "Oui")),
    atcd_covid = factor(ifelse(str_detect(antecedents, regex("covid", TRUE)), "Oui", "Non"), c("Non", "Oui")),

    # Pratiques à risque (champ libre ; vide = aucune pratique déclarée)
    pratiques = ifelse(is.na(pratiques), "", pratiques),
    pgae = factor(ifelse(str_detect(pratiques,
                  regex("intubation|suction|aspiration|bronchoscop|sputum|expector|exprector|r[ée]animation", TRUE)),
                  "Oui", "Non"), c("Non", "Oui")),   # procédures générant des aérosols
    prelevements = factor(ifelse(str_detect(pratiques, regex("sample|sampling|prél", TRUE)), "Oui", "Non"), c("Non", "Oui")),

    # Groupes d'âge et d'ancienneté (pour la description uniquement)
    age_cl = cut(age, c(-Inf, 24, 29, Inf), labels = c("< 25 ans", "25–29 ans", "≥ 30 ans")),
    exp_cl = cut(experience, c(-Inf, 1, 4, Inf), labels = c("≤ 1 an", "> 1 à 4 ans", "> 4 ans"))
  )


## 3. Contrôle de qualité --------------------------------------------------------
cat("\nDoublons d'identifiant :", sum(duplicated(d$id)), "\n")
cat("QFT indéterminé (exclu) :", sum(is.na(d$qft_pos)), "\n")
# Cohérence entre le préfixe de l'identifiant (ITL = cas, TS = témoin sain) et le QFT
d$prefixe <- str_extract(d$id, "(?<=-)[A-Z]+")
print(table(Prefixe = d$prefixe, QFT = d$qft, useNA = "ifany"))
print(d %>% filter((prefixe %in% c("ITL", "ILT") & qft != "QFT+") |
                   (prefixe == "TS" & qft != "QFT-")) %>% select(id, qft))
# Valeurs plausibles
print(summary(d[, c("age", "experience")]))
cat("Ancienneté > âge - 16 ans :", sum(d$experience > d$age - 16, na.rm = TRUE), "\n")

# Population d'analyse : QFT interprétable
a <- d %>% filter(!is.na(qft_pos)) %>%
  mutate(QFT = factor(qft_pos, 0:1, c("QFT-", "QFT+")))
cat("Population d'analyse : n =", nrow(a), "\n")

vars_cat <- c("sexe", "age_cl", "profession_g", "service_g", "exposition", "exp_cl",
              "masque", "masque_reg", "formation", "bcg", "cicatrice", "tabac",
              "atcd_tb", "atcd_covid", "pgae", "prelevements")
vars_num <- c("age", "experience")
libelles <- c(sexe = "Sexe", age_cl = "Classe d'âge", profession_g = "Profession",
  service_g = "Service", exposition = "Exposition déclarée à la tuberculose",
  exp_cl = "Ancienneté professionnelle", masque = "Port du masque",
  masque_reg = "Port du masque (2 classes)", formation = "Formation en biosécurité",
  bcg = "Vaccination BCG", cicatrice = "Cicatrice BCG", tabac = "Tabagisme (actuel ou ancien)",
  atcd_tb = "Antécédent de tuberculose", atcd_covid = "Antécédent de COVID-19",
  pgae = "Procédures générant des aérosols", prelevements = "Manipulation de prélèvements",
  age = "Âge (ans)", experience = "Ancienneté (ans)")

# Données manquantes par variable
manq <- data.frame(Variable = unname(libelles[c(vars_num, vars_cat)]),
                   n_manquant = sapply(a[c(vars_num, vars_cat)], function(x) sum(is.na(x)))) %>%
  mutate(pourcentage = round(100 * n_manquant / nrow(a), 1))
print(manq, row.names = FALSE)

## 4. Q1 – Prévalence de l'ITL ---------------------------------------------------
n_pos <- sum(a$qft_pos); n_tot <- nrow(a)
wilson <- prop.test(n_pos, n_tot, correct = FALSE)$conf.int
exact  <- binom.test(n_pos, n_tot)$conf.int
prevalence <- data.frame(
  n_positifs = n_pos, n_total = n_tot, prevalence_pct = round(100 * n_pos / n_tot, 1),
  IC95_Wilson = sprintf("%.1f–%.1f", 100 * wilson[1], 100 * wilson[2]),
  IC95_exact  = sprintf("%.1f–%.1f", 100 * exact[1], 100 * exact[2]))
print(prevalence)

## 5. Description et analyse bivariée (Q2) ------------------------------------------
# Variables quantitatives : normalité (Shapiro-Wilk) puis test adapté
desc_num <- lapply(vars_num, function(v) {
  x <- a[[v]]; g <- a$QFT
  sw <- by(x, g, function(z) shapiro.test(z[!is.na(z)])$p.value)
  normal <- all(unlist(sw) > 0.05)
  test <- if (normal) t.test(x ~ g) else wilcox.test(x ~ g, exact = FALSE)
  med <- function(z) sprintf("%.1f [%.1f–%.1f]", median(z, na.rm = TRUE),
                             quantile(z, .25, na.rm = TRUE), quantile(z, .75, na.rm = TRUE))
  data.frame(Variable = unname(libelles[v]),
             Total = med(x), `QFT-` = med(x[g == "QFT-"]), `QFT+` = med(x[g == "QFT+"]),
             Test = ifelse(normal, "t de Student", "Mann-Whitney"),
             p = test$p.value, n_manquant = sum(is.na(x)), check.names = FALSE)
}) %>% bind_rows()

# Variables qualitatives : chi² si tous les effectifs attendus ≥ 5, sinon Fisher
desc_cat <- lapply(vars_cat, function(v) {
  tab <- table(a[[v]], a$QFT)
  attendu <- suppressWarnings(chisq.test(tab)$expected)
  use_fisher <- any(attendu < 5)
  p <- if (use_fisher) fisher.test(tab, simulate.p.value = nrow(tab) > 2, B = 1e5)$p.value
       else chisq.test(tab, correct = FALSE)$p.value
  tot <- rowSums(tab)
  data.frame(Variable = unname(libelles[v]), Modalite = rownames(tab),
             Total = sprintf("%d (%.1f)", tot, 100 * tot / sum(tot)),
             `QFT-` = sprintf("%d (%.1f)", tab[, 1], 100 * tab[, 1] / sum(tab[, 1])),
             `QFT+` = sprintf("%d (%.1f)", tab[, 2], 100 * tab[, 2] / sum(tab[, 2])),
             `% QFT+ dans la ligne` = sprintf("%.1f", 100 * tab[, 2] / tot),
             Test = ifelse(use_fisher, "Fisher exact", "Chi²"),
             p = c(p, rep(NA, nrow(tab) - 1)), n_manquant = c(sum(is.na(a[[v]])), rep(NA, nrow(tab) - 1)),
             check.names = FALSE)
}) %>% bind_rows()

# Tendance pour le port du masque (variable ordinale) : test de Cochran-Armitage
tm <- table(a$masque, a$QFT)
p_tendance_masque <- prop.trend.test(tm[, 2], rowSums(tm))$p.value
cat("Test de tendance (port du masque) : p =", signif(p_tendance_masque, 3), "\n")

# OR bruts (régression logistique simple, IC par vraisemblance profilée)
or_brut <- function(v) {
  m <- glm(reformulate(v, "qft_pos"), data = a, family = binomial)
  ci <- suppressMessages(exp(confint(m)))
  data.frame(Variable = unname(libelles[v]), Terme = names(coef(m))[-1],
             OR = exp(coef(m))[-1], IC_bas = ci[-1, 1], IC_haut = ci[-1, 2],
             p_Wald = summary(m)$coefficients[-1, 4],
             p_global_LRT = c(drop1(m, test = "LRT")[2, "Pr(>Chi)"], rep(NA, length(coef(m)) - 2)),
             n = nobs(m))
}
vars_or <- c("age", "experience", "sexe", "profession_g", "service_g", "exposition",
             "masque_reg", "formation", "bcg", "cicatrice", "tabac", "atcd_covid",
             "pgae", "prelevements")
tab_or_brut <- bind_rows(lapply(vars_or, or_brut))
print(tab_or_brut %>% mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

## 6. Q3 – Régression logistique multivariable ---------------------------------
# 6.1 Choix des variables A PRIORI (connaissances, pas p-values)
#   Âge (exposition cumulée), sexe, exposition déclarée à la tuberculose,
#   procédures générant des aérosols, port régulier du masque, vaccination BCG.
#   Le QFT n'est pas faussé par le BCG (antigènes ESAT-6/CFP-10 absents du BCG) ;
#   le BCG est inclus comme facteur potentiellement protecteur contre l'infection.
#   L'ancienneté n'est pas incluse : fortement liée à l'âge et 16 % manquante.
vars_mod <- c("age", "sexe", "exposition", "pgae", "masque_reg", "bcg")
f_mod <- reformulate(vars_mod, "qft_pos")

# 6.2 Nombre d'événements par paramètre (EPV)
n_param <- 6
cat("\nÉvénements (QFT+) :", n_pos, "| paramètres :", n_param,
    "| EPV =", round(n_pos / n_param, 1), "\n")

# 6.3 Analyse en cas complets
cc <- a %>% select(qft_pos, all_of(vars_mod)) %>% na.omit()
cat("Cas complets : n =", nrow(cc), "| événements =", sum(cc$qft_pos),
    "| EPV =", round(sum(cc$qft_pos) / n_param, 1), "\n")
m_cc <- glm(f_mod, data = cc, family = binomial)

tidy_or <- function(m, methode) {
  t <- broom::tidy(m, conf.int = TRUE, exponentiate = TRUE, conf.method = "profile")
  t %>% filter(term != "(Intercept)") %>%
    transmute(Methode = methode, Terme = term, OR = estimate,
              IC_bas = conf.low, IC_haut = conf.high, p = p.value)
}
res_cc <- tidy_or(m_cc, "Cas complets (MV)")

# 6.4 Diagnostics du modèle
# Colinéarité
print(car::vif(m_cc))
# Linéarité de l'âge sur le logit : spline cubique restreinte (3 nœuds) vs linéaire
dd <- datadist(cc); options(datadist = "dd")
m_lin <- lrm(f_mod, data = cc, x = TRUE, y = TRUE)
m_rcs <- lrm(qft_pos ~ rcs(age, 3) + sexe + exposition + pgae + masque_reg + bcg, data = cc)
p_lin <- 1 - pchisq(m_rcs$stats["Model L.R."] - m_lin$stats["Model L.R."], 1)
cat("Test de non-linéarité de l'âge : p =", signif(p_lin, 3), "\n")
# Observations influentes
infl <- which(cooks.distance(m_cc) > 4 / nrow(cc))
cat("Observations avec distance de Cook > 4/n :", length(infl), "\n")

# 6.5 Performance : discrimination, calibration, validation interne
roc_cc <- roc(cc$qft_pos, fitted(m_cc), quiet = TRUE, ci = TRUE)
cat("AUC apparente =", sprintf("%.3f (IC95 %.3f–%.3f)", roc_cc$auc,
                                roc_cc$ci[1], roc_cc$ci[3]), "\n")
# Test global d'ajustement de le Cessie–van Houwelingen (alternative à Hosmer-Lemeshow)
gof <- residuals(m_lin, type = "gof"); print(gof)
# Hosmer–Lemeshow (10 groupes), codé directement
hosmer_lemeshow <- function(y, p, g = 10) {
  grp <- cut(p, unique(quantile(p, seq(0, 1, length.out = g + 1))), include.lowest = TRUE)
  obs <- tapply(y, grp, sum); att <- tapply(p, grp, sum); n <- tapply(y, grp, length)
  stat <- sum((obs - att)^2 / (att * (1 - att / n)))
  c(chi2 = stat, ddl = nlevels(grp) - 2, p = 1 - pchisq(stat, nlevels(grp) - 2))
}
print(hosmer_lemeshow(cc$qft_pos, fitted(m_cc)))
# Validation interne par bootstrap (B = 500) : optimisme de l'AUC et pente de calibration
val <- validate(m_lin, method = "boot", B = 500)
auc_corr <- 0.5 + val["Dxy", "index.corrected"] / 2
cat("AUC corrigée de l'optimisme =", round(auc_corr, 3),
    "| pente de calibration corrigée =", round(val["Slope", "index.corrected"], 3), "\n")
cal <- calibrate(m_lin, method = "boot", B = 500)
brier <- mean((fitted(m_cc) - cc$qft_pos)^2)
cat("Score de Brier =", round(brier, 3), "\n")

# 6.6 Régression de Firth (biais des petits effectifs) – sensibilité
m_firth <- glm(f_mod, data = cc, family = binomial, method = "brglmFit", type = "AS_mean")
res_firth <- suppressWarnings(broom::tidy(m_firth)) %>% filter(term != "(Intercept)") %>%
  transmute(Methode = "Firth (cas complets)", Terme = term, OR = exp(estimate),
            IC_bas = exp(estimate - 1.96 * std.error), IC_haut = exp(estimate + 1.96 * std.error),
            p = p.value)

# 6.7 Imputation multiple (analyse principale ; hypothèse MAR)
vars_imp <- c("qft_pos", vars_mod, "experience", "profession_g", "service_g", "formation", "tabac")
imp_data <- a %>% select(all_of(vars_imp))
imp <- mice(imp_data, m = 30, maxit = 20, printFlag = FALSE, seed = 2026)
fit_mi <- with(imp, glm(qft_pos ~ age + sexe + exposition + pgae + masque_reg + bcg, family = binomial))
pool_mi <- as.data.frame(summary(pool(fit_mi), conf.int = TRUE))
res_mi <- pool_mi %>% filter(term != "(Intercept)") %>%
  transmute(Methode = "Imputation multiple (principale)", Terme = as.character(term),
            OR = exp(estimate), IC_bas = exp(`2.5 %`), IC_haut = exp(`97.5 %`), p = p.value)

res_multi <- bind_rows(res_mi, res_cc, res_firth)
print(as.data.frame(res_multi) %>% mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)

# 6.8 Sensibilité : ajout du service (modèle exploratoire, mêmes patients)
d_s <- a %>% select(qft_pos, all_of(vars_mod), service_g) %>% na.omit()
m_sans <- glm(f_mod, data = d_s, family = binomial)
m_avec <- update(m_sans, . ~ . + service_g)
p_service <- anova(m_sans, m_avec, test = "LRT")[2, "Pr(>Chi)"]
cat("Ajout du service au modèle (LRT, 5 ddl) : p =", signif(p_service, 3), "\n")

## 7. Figures -------------------------------------------------------------------
# 7.1 Prévalence du QFT+ par service (IC95 de Wilson)
prev_serv <- a %>% filter(!is.na(service_g)) %>% group_by(service_g) %>%
  summarise(n = n(), pos = sum(qft_pos), .groups = "drop") %>%
  rowwise() %>%
  mutate(prev = pos / n,   # IC de Wilson (avertissement normal pour les petits effectifs)
         bas  = suppressWarnings(prop.test(pos, n, correct = FALSE)$conf.int[1]),
         haut = suppressWarnings(prop.test(pos, n, correct = FALSE)$conf.int[2])) %>% ungroup()
g1 <- ggplot(prev_serv, aes(x = reorder(service_g, prev), y = 100 * prev)) +
  geom_pointrange(aes(ymin = 100 * bas, ymax = 100 * haut)) +
  geom_hline(yintercept = 100 * n_pos / n_tot, linetype = 2) +
  geom_text(aes(label = paste0(pos, "/", n)), nudge_x = 0.3, size = 3) +
  coord_flip() + labs(x = NULL, y = "Prévalence QFT+ (%) et IC à 95 %",
  title = "Prévalence de l'ITL par service", subtitle = "Ligne pointillée : prévalence globale") +
  theme_minimal(base_size = 11)
ggsave(file.path(dossier, "fig1_prevalence_service.png"), g1, width = 7, height = 4.5, dpi = 300, bg = "white")

# 7.2 Forest plot des OR ajustés (imputation multiple)
noms_termes <- c(age = "Âge (par année)", sexeHomme = "Homme vs femme",
  expositionOui = "Exposition déclarée à la TB", pgaeOui = "Procédures générant des aérosols",
  "masque_regIrrégulier" = "Port du masque irrégulier", bcgNon = "Non vacciné BCG")
g2 <- res_multi %>% mutate(Terme = noms_termes[Terme]) %>%
  ggplot(aes(x = OR, y = Terme, colour = Methode)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = IC_bas, xmax = IC_haut), position = position_dodge(width = 0.6)) +
  scale_x_log10() + labs(x = "Odds ratio ajusté (échelle log) et IC à 95 %", y = NULL, colour = NULL,
  title = "Facteurs associés au QFT+ : modèle multivariable") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(dossier, "fig2_forest_plot.png"), g2, width = 7.5, height = 5.5, dpi = 300, bg = "white")

# 7.3 Courbe ROC et 7.4 courbe de calibration (bootstrap)
png(file.path(dossier, "fig3_roc.png"), width = 1600, height = 1600, res = 300, bg = "white")
plot(roc_cc, print.auc = TRUE, legacy.axes = TRUE, main = "Courbe ROC (cas complets)")
dev.off()
png(file.path(dossier, "fig4_calibration.png"), width = 1800, height = 1600, res = 300, bg = "white")
plot(cal, main = "Calibration (bootstrap, B = 500)", xlab = "Probabilité prédite", ylab = "Proportion observée")
dev.off()

## 8. Export des tableaux --------------------------------------------------------
arrondir <- function(df) df %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
perf <- data.frame(
  Indicateur = c("n (cas complets)", "Événements", "EPV", "AUC apparente", "IC95 AUC",
                 "AUC corrigée (bootstrap)", "Pente de calibration corrigée",
                 "Hosmer-Lemeshow p", "le Cessie-van Houwelingen p", "Score de Brier",
                 "Non-linéarité de l'âge p", "VIF maximal", "Ajout du service (LRT) p"),
  Valeur = c(nrow(cc), sum(cc$qft_pos), round(sum(cc$qft_pos) / n_param, 1), round(roc_cc$auc, 3),
             sprintf("%.3f–%.3f", roc_cc$ci[1], roc_cc$ci[3]), round(auc_corr, 3),
             round(val["Slope", "index.corrected"], 3),
             round(hosmer_lemeshow(cc$qft_pos, fitted(m_cc))["p"], 3), round(gof["P"], 3),
             round(brier, 3), round(p_lin, 3), round(max(car::vif(m_cc)), 2), round(p_service, 3)))
wb <- list("Données manquantes" = manq, "Prévalence" = prevalence,
           "Bivarié quantitatif" = arrondir(desc_num), "Bivarié qualitatif" = arrondir(desc_cat),
           "OR bruts" = arrondir(tab_or_brut), "OR ajustés" = arrondir(as.data.frame(res_multi)),
           "Performance modèle" = perf)
openxlsx::write.xlsx(wb, file.path(dossier, "tableaux_resultats.xlsx"))

## 9. Reproductibilité -----------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(dossier, "sessionInfo.txt"))
cat("\nAnalyse terminée. Résultats dans le dossier :", dossier, "\n")

## 10. Graphics
dir_fig <- file.path(dossier, "figures_descriptives")
dir.create(dir_fig, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 11))
sauver <- function(g, nom, w = 7, h = 4.5) {
  ggsave(file.path(dir_fig, nom), g, width = w, height = h, dpi = 300, bg = "white")
}
couleurs_qft <- c("QFT-" = "#4C72B0", "QFT+" = "#C44E52")

## 10.1 Tableau descriptif global des variables quantitatives -------------------
desc_quanti <- lapply(vars_num, function(v) {
  x <- a[[v]]; xn <- x[!is.na(x)]
  data.frame(
    Variable   = unname(libelles[v]),
    n          = length(xn),
    Manquants  = sum(is.na(x)),
    Moyenne    = mean(xn),
    Ecart_type = sd(xn),
    Mediane    = median(xn),
    Q1         = unname(quantile(xn, 0.25)),
    Q3         = unname(quantile(xn, 0.75)),
    Minimum    = min(xn),
    Maximum    = max(xn),
    Asymetrie  = mean((xn - mean(xn))^3) / sd(xn)^3,     # > 1 : forte asymétrie à droite
    Shapiro_p  = shapiro.test(xn)$p.value)               # < 0,05 : distribution non normale
}) %>% bind_rows()
print(desc_quanti %>% mutate(across(where(is.numeric), ~ round(.x, 2))), row.names = FALSE)

## 10.2 Tableau descriptif global des variables qualitatives --------------------
desc_quali <- lapply(vars_cat, function(v) {
  x <- a[[v]]; tab <- table(x)
  data.frame(Variable = unname(libelles[v]), Modalite = names(tab),
             n = as.integer(tab),
             Pourcentage = round(100 * as.integer(tab) / sum(tab), 1),   # sur les répondants
             Manquants = c(sum(is.na(x)), rep(NA, length(tab) - 1)))
}) %>% bind_rows()
print(desc_quali, row.names = FALSE)

## 10.3 Distribution des variables quantitatives : histogramme + densité --------
#      (la ligne rouge = moyenne, la ligne bleue = médiane)
for (v in vars_num) {
  df <- a %>% filter(!is.na(.data[[v]]))
  g <- ggplot(df, aes(x = .data[[v]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "grey75", colour = "white") +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = mean(df[[v]]), colour = "#C44E52", linetype = 2) +
    geom_vline(xintercept = median(df[[v]]), colour = "#4C72B0", linetype = 2) +
    labs(x = libelles[v], y = "Densité",
         title = paste("Distribution :", libelles[v]),
         subtitle = sprintf("n = %d ; moyenne (rouge) = %.1f ; médiane (bleu) = %.1f",
                            nrow(df), mean(df[[v]]), median(df[[v]])))
  sauver(g, paste0("histogramme_", v, ".png"))
}

## 10.4 Q-Q plots : vérification visuelle de la normalité -----------------------
qq <- a %>% select(all_of(vars_num)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "valeur") %>%
  filter(!is.na(valeur)) %>% mutate(variable = libelles[variable])
g_qq <- ggplot(qq, aes(sample = valeur)) + stat_qq() + stat_qq_line(colour = "#C44E52") +
  facet_wrap(~ variable, scales = "free") +
  labs(x = "Quantiles théoriques (loi normale)", y = "Quantiles observés",
       title = "Q-Q plots : les points doivent suivre la droite si la distribution est normale")
sauver(g_qq, "qqplot_quantitatives.png", w = 8, h = 4)

## 10.5 Boîtes à moustaches selon le résultat du QFT -----------------------------
bx <- a %>% select(QFT, all_of(vars_num)) %>%
  pivot_longer(-QFT, names_to = "variable", values_to = "valeur") %>%
  filter(!is.na(valeur)) %>% mutate(variable = libelles[variable])
g_box <- ggplot(bx, aes(x = QFT, y = valeur, fill = QFT)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values = couleurs_qft, guide = "none") +
  labs(x = NULL, y = NULL, title = "Âge et ancienneté selon le résultat du QFT",
       subtitle = "Boîte = médiane et quartiles ; points = participants")
sauver(g_box, "boxplots_selon_QFT.png", w = 8, h = 4.5)

## 10.6 Répartition des variables qualitatives (diagrammes en barres) -----------
vars_barres <- c("sexe", "age_cl", "profession_g", "service_g", "exposition", "masque",
                 "formation", "bcg", "cicatrice", "tabac", "pgae", "atcd_covid")
for (v in vars_barres) {
  df <- a %>% count(modalite = .data[[v]]) %>%
    mutate(modalite = fct_na_value_to_level(factor(modalite), "Manquant"),
           pct = 100 * n / sum(n))
  g <- ggplot(df, aes(x = modalite, y = n)) +
    geom_col(fill = ifelse(df$modalite == "Manquant", "grey80", "#4C72B0")) +
    geom_text(aes(label = sprintf("%d (%.0f %%)", n, pct)), vjust = -0.4, size = 3.2) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Effectif", title = libelles[v],
         subtitle = paste("n =", nrow(a), "; les données manquantes sont en gris")) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  sauver(g, paste0("barres_", v, ".png"))
}

## 10.7 Proportion de QFT+ selon chaque variable qualitative (barres empilées à 100 %)
for (v in vars_barres) {
  df <- a %>% filter(!is.na(.data[[v]])) %>% count(modalite = .data[[v]], QFT) %>%
    group_by(modalite) %>% mutate(pct = 100 * n / sum(n), total = sum(n)) %>% ungroup()
  g <- ggplot(df, aes(x = modalite, y = pct, fill = QFT)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.0f %%", pct)), position = position_stack(vjust = 0.5),
              colour = "white", size = 3.2) +
    geom_text(data = distinct(df, modalite, total), inherit.aes = FALSE,
              aes(x = modalite, y = 103, label = paste0("n=", total)), size = 3) +
    scale_fill_manual(values = couleurs_qft) +
    labs(x = NULL, y = "% des participants", fill = NULL,
         title = paste("Résultat du QFT selon :", libelles[v])) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")
  sauver(g, paste0("QFT_selon_", v, ".png"))
}

## 10.8 Répartition du résultat QFT dans toute la population ---------------------
df_qft <- d %>% mutate(resultat = case_when(qft == "QFT+" ~ "QFT+", qft == "QFT-" ~ "QFT-",
                                            TRUE ~ "Indéterminé")) %>%
  count(resultat) %>% mutate(pct = 100 * n / sum(n),
                             resultat = factor(resultat, c("QFT-", "QFT+", "Indéterminé")))
g_qft <- ggplot(df_qft, aes(x = resultat, y = n, fill = resultat)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d (%.1f %%)", n, pct)), vjust = -0.4) +
  scale_fill_manual(values = c(couleurs_qft, "Indéterminé" = "grey70"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Effectif", title = "Résultat du QuantiFERON (n = 165)")
sauver(g_qft, "repartition_QFT.png", w = 5.5, h = 4.5)

## 10.9 Pyramide des âges par sexe ------------------------------------------------
pyr <- a %>% mutate(classe = cut(age, c(15, 20, 25, 30, 35, 40, 45, 60), right = FALSE,
                                 labels = c("15–19", "20–24", "25–29", "30–34", "35–39", "40–44", "≥ 45"))) %>%
  count(classe, sexe) %>% mutate(n = ifelse(sexe == "Homme", -n, n))
g_pyr <- ggplot(pyr, aes(x = classe, y = n, fill = sexe)) + geom_col(width = 0.8) +
  coord_flip() + scale_y_continuous(labels = abs) +
  scale_fill_manual(values = c(Femme = "#DD8452", Homme = "#55A868")) +
  labs(x = "Âge (ans)", y = "Effectif", fill = NULL, title = "Pyramide des âges par sexe") +
  theme(legend.position = "top")
sauver(g_pyr, "pyramide_ages.png", w = 6, h = 4.5)

## 10.10 Prévalence du QFT+ par profession avec IC à 95 % (Wilson) ---------------
prev_prof <- a %>% filter(!is.na(profession_g)) %>% group_by(profession_g) %>%
  summarise(n = n(), pos = sum(qft_pos), .groups = "drop") %>% rowwise() %>%
  mutate(prev = pos / n,
         bas  = suppressWarnings(prop.test(pos, n, correct = FALSE)$conf.int[1]),
         haut = suppressWarnings(prop.test(pos, n, correct = FALSE)$conf.int[2])) %>% ungroup()
g_prof <- ggplot(prev_prof, aes(x = reorder(profession_g, prev), y = 100 * prev)) +
  geom_pointrange(aes(ymin = 100 * bas, ymax = 100 * haut)) +
  geom_hline(yintercept = 100 * mean(a$qft_pos), linetype = 2) +
  geom_text(aes(label = paste0(pos, "/", n)), nudge_x = 0.3, size = 3) +
  coord_flip() + labs(x = NULL, y = "Prévalence QFT+ (%) et IC à 95 %",
                      title = "Prévalence de l'ITL par profession",
                      subtitle = "Ligne pointillée : prévalence globale")
sauver(g_prof, "prevalence_profession.png")

## 10.11 Carte des données manquantes ---------------------------------------------
vars_manq <- c(vars_num, "sexe", "profession_g", "service_g", "exposition", "masque",
               "formation", "bcg", "cicatrice", "tabac")
carte <- a %>% mutate(ligne = row_number()) %>%
  select(ligne, all_of(vars_manq)) %>%
  mutate(across(-ligne, ~ is.na(.x))) %>%
  pivot_longer(-ligne, names_to = "variable", values_to = "manquant") %>%
  mutate(variable = factor(libelles[variable], levels = rev(libelles[vars_manq])))
g_carte <- ggplot(carte, aes(x = ligne, y = variable, fill = manquant)) + geom_tile() +
  scale_fill_manual(values = c(`FALSE` = "grey90", `TRUE` = "#C44E52"),
                    labels = c("Renseigné", "Manquant")) +
  labs(x = "Participant", y = NULL, fill = NULL, title = "Carte des données manquantes") +
  theme(legend.position = "top", panel.grid = element_blank())
sauver(g_carte, "carte_donnees_manquantes.png", w = 9, h = 4.5)

g_manq <- manq %>% filter(n_manquant > 0) %>%
  ggplot(aes(x = reorder(Variable, pourcentage), y = pourcentage)) +
  geom_col(fill = "#C44E52") + geom_text(aes(label = paste0(pourcentage, " %")), hjust = -0.2, size = 3.2) +
  coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% de données manquantes", title = "Données manquantes par variable")
sauver(g_manq, "pourcentage_manquants.png")

## 10.12 Export Excel des tableaux descriptifs -----------------------------------
openxlsx::write.xlsx(list("Quantitatives" = desc_quanti %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
                          "Qualitatives"  = desc_quali),
                     file.path(dossier, "tableaux_descriptifs.xlsx"))
cat("\nAnalyse descriptive terminée :", length(list.files(dir_fig)), "figures dans", dir_fig, "\n")


