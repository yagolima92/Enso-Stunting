# ============================================================
# ENSO & Childhood Stunting in Brazil — Full Analysis Script
# Revised version — all four ENSO-FOOD CPI specifications
# ============================================================
# PART I  — National mechanisms (2000–2025)
# PART II — ENSO → FOOD CPI: four model specifications
# PART III— Climatic pathways (2010–2024)
# PART IV — State-level panel (2010–2024)
# PART V  — Regional summary
# ============================================================

library(tidyverse)
library(lmtest)
library(sandwich)
library(plm)
library(mgcv)

# ============================================================
# DATA INPUT
# ============================================================

national <- data.frame(
  year       = 2000:2024,
  pou        = c(12.12,10.70,9.30,7.90,6.80,6.50,6.20,5.60,5.10,4.30,
                 3.70,3.20,3.00,2.80,2.60,2.50,2.50,2.50,2.50,2.60,
                 4.10,3.40,3.20,2.5,2.5),
  stunting   = c(9.9,9.2,8.6,8.1,7.7,7.4,7.2,7.0,6.9,6.8,
                 6.7,6.6,6.6,6.7,6.8,6.9,7.1,7.3,7.5,7.7,
                 8.0,8.2,8.5,8.7,8.9),
  anemia_mom = c(30.9,30.0,29.2,28.5,27.7,26.9,26.1,25.3,24.4,23.5,
                 22.7,22.0,21.5,21.1,20.8,20.6,20.5,20.5,20.4,20.4,
                 20.4,20.6,20.9,21.3,20.7),
  meat       = c(20.87,20.68,20.54,20.81,23.26,25.91,26.47,24.71,23.46,24.95,
                 25.24,24.99,25.70,29.70,27.68,29.72,32.08,30.76,26.99,25.67,
                 26.32,18.00,17.00,19.00,18.00),
  enso       = c(-0.65,-0.22,0.96,0.12,0.53,-0.20,0.46,-0.97,-0.48,0.77,
                 -1.24,-0.76,0.12,-0.29,0.33,1.92,-0.42,-0.27,0.40,0.37,
                 -0.77,-0.64,-0.92,1.38,-0.13),
  ipca       = c(32.18,35.27,42.14,45.29,47.04,47.97,48.57,53.80,59.78,61.67,
                 68.08,72.98,80.18,86.97,93.96,105.24,114.30,112.17,116.69,124.14,
                 141.63,152.89,170.70,172.40,185.70)
)

# Lagged ENSO (t-1)
national$enso_lag1 <- c(NA, national$enso[-nrow(national)])

# Detrending function
detrend <- function(x, year) resid(lm(x ~ year))

national_dt <- national %>%
  mutate(
    stunting_dt = detrend(stunting,   year),
    pou_dt      = detrend(pou,        year),
    anemia_dt   = detrend(anemia_mom, year),
    meat_dt     = detrend(meat,       year),
    ipca_dt     = detrend(ipca,       year),
    enso_dt     = detrend(enso,       year)
  )
national_dt$enso_lag1_dt <- c(NA, national_dt$enso_dt[-nrow(national_dt)])

cat("=== Data loaded:", nrow(national), "observations ===\n\n")

# ============================================================
# PART I — NATIONAL STRUCTURAL MECHANISMS (detrended)
# ============================================================
cat("============================================================\n")
cat("PART I — NATIONAL STRUCTURAL MECHANISMS (2000–2025)\n")
cat("============================================================\n\n")

nat_dt <- national_dt %>% filter(!is.na(enso_lag1_dt))

# --- Correlations (detrended) --------------------------------
cat("--- PEARSON CORRELATIONS (detrended series) ---\n")
cor_pou_stunt  <- cor.test(national_dt$pou_dt,       national_dt$stunting_dt)
cor_anem_stunt <- cor.test(national_dt$anemia_dt,    national_dt$stunting_dt)
cor_meat_stunt <- cor.test(national_dt$meat_dt,      national_dt$stunting_dt)
cor_ipca_meat  <- cor.test(national_dt$ipca_dt,      national_dt$meat_dt)

cat(sprintf("PoU    ↔ Stunting : r = %.3f (p = %.4f)\n",
            cor_pou_stunt$estimate,  cor_pou_stunt$p.value))
cat(sprintf("Anemia ↔ Stunting : r = %.3f (p = %.4f)\n",
            cor_anem_stunt$estimate, cor_anem_stunt$p.value))
cat(sprintf("Meat   ↔ Stunting : r = %.3f (p = %.4f)\n",
            cor_meat_stunt$estimate, cor_meat_stunt$p.value))
cat(sprintf("IPCA   → Meat     : r = %.3f (p = %.4f)\n\n",
            cor_ipca_meat$estimate,  cor_ipca_meat$p.value))

# --- IPCA → Meat (detrended) ---------------------------------
cat("--- MODEL: FOOD CPI → Meat consumption (detrended) ---\n")
m_ipca_meat <- lm(meat_dt ~ ipca_dt, data = national_dt)
print(summary(m_ipca_meat))

# --- Meat → Stunting -----------------------------------------
cat("\n--- MODEL: Meat → Stunting ---\n")
m_meat_stunt <- lm(stunting ~ meat, data = national)
print(summary(m_meat_stunt))

# --- ENSO → Stunting (direct, detrended) --------------------
cat("\n--- MODEL: ENSO(t-1) → Stunting (direct, detrended) ---\n")
m_enso_stunt <- lm(stunting_dt ~ enso_lag1_dt, data = nat_dt)
print(summary(m_enso_stunt))

# --- Final multivariate model (detrended) --------------------
cat("\n--- FINAL MULTIVARIATE MODEL (detrended) ---\n")
cat("Stunting ~ PoU + Anemia + Meat + IPCA + ENSO(t-1)\n\n")
m_final <- lm(stunting_dt ~ pou_dt + anemia_dt + meat_dt + ipca_dt + enso_lag1_dt,
              data = nat_dt)
print(summary(m_final))

# ============================================================
# PART II — ENSO → FOOD CPI: FOUR MODEL SPECIFICATIONS
# ============================================================
cat("\n============================================================\n")
cat("PART II — ENSO → FOOD CPI: FOUR MODEL SPECIFICATIONS\n")
cat("(Following Hsiang & Meng, 2015)\n")
cat("============================================================\n\n")

df <- national %>%
  filter(!is.na(enso_lag1)) %>%
  mutate(
    enso2          = enso^2,
    enso3          = enso^3,
    enso_lag2      = enso_lag1^2,
    enso_lag3      = enso_lag1^3,
    d_strong_nina  = as.integer(enso <= -1.0),
    d_weak_nina    = as.integer(enso > -1.0 & enso <= -0.5),
    d_neutral      = as.integer(enso > -0.5 & enso < 0.5),
    d_weak_nino    = as.integer(enso >= 0.5 & enso < 1.0),
    d_strong_nino  = as.integer(enso >= 1.0),
    dl_strong_nina = as.integer(enso_lag1 <= -1.0),
    dl_weak_nina   = as.integer(enso_lag1 > -1.0 & enso_lag1 <= -0.5),
    dl_neutral     = as.integer(enso_lag1 > -0.5 & enso_lag1 < 0.5),
    dl_weak_nino   = as.integer(enso_lag1 >= 0.5 & enso_lag1 < 1.0),
    dl_strong_nino = as.integer(enso_lag1 >= 1.0),
    t              = year - min(year)
  )

# --- Model 1: Linear -----------------------------------------
cat("--- MODEL 1: LINEAR (baseline Hsiang & Meng) ---\n")
m_lin <- lm(ipca ~ enso + enso_lag1 + t, data = df)
print(summary(m_lin))

# --- Model 2: Cubic ------------------------------------------
cat("\n--- MODEL 2: CUBIC ---\n")
m_cub <- lm(ipca ~ enso + enso2 + enso3 +
                   enso_lag1 + enso_lag2 + enso_lag3 + t, data = df)
print(summary(m_cub))

# --- Model 3: GAM --------------------------------------------
cat("\n--- MODEL 3: GAM ---\n")
m_gam <- gam(ipca ~ s(enso) + t, data = df)
print(summary(m_gam))
cat("\n--- ANOVA: GAM vs Linear ---\n")
m_gam_lin <- gam(ipca ~ enso + t, data = df)
print(anova(m_gam_lin, m_gam, test = "F"))

# --- Model 4: Dummy (preferred) ------------------------------
cat("\n--- MODEL 4: DUMMY (preferred specification) ---\n")
m_dum <- lm(ipca ~ d_weak_nina + d_neutral + d_weak_nino + d_strong_nino +
                   dl_weak_nina + dl_neutral + dl_weak_nino + dl_strong_nino + t,
            data = df)
print(summary(m_dum))

cat("\n--- SUMMARY: Dummy model significant coefficients ---\n")
coefs <- summary(m_dum)$coefficients
sig   <- coefs[coefs[,4] < 0.10, ]
print(round(sig, 4))

# ============================================================
# PART III — CLIMATIC PATHWAYS (2010–2024)
# ============================================================
cat("\n============================================================\n")
cat("PART III — CLIMATIC PATHWAYS: ENSO → TEMP / PREC (2010–2024)\n")
cat("============================================================\n\n")

states <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA",
            "MG","MS","MT","PA","PB","PE","PI","PR","RJ","RN",
            "RO","RR","RS","SC","SE","SP","TO")

region_map <- c(
  AC="Norte", AL="Nordeste", AM="Norte",  AP="Norte",  BA="Nordeste",
  CE="Nordeste", DF="Centro-Oeste", ES="Sudeste", GO="Centro-Oeste",
  MA="Nordeste", MG="Sudeste", MS="Centro-Oeste", MT="Centro-Oeste",
  PA="Norte", PB="Nordeste", PE="Nordeste", PI="Nordeste",
  PR="Sul", RJ="Sudeste", RN="Nordeste", RO="Norte",
  RR="Norte", RS="Sul", SC="Sul", SE="Nordeste", SP="Sudeste", TO="Norte"
)

stunting_state <- matrix(c(
  23.45,22.81,22.15,21.29,22.92,18.95,19.74,17.86,18.48,19.63,18.27,18.02,17.27,17.57,16.1,
  16.68,15.01,13.79,13.43,14.35,13.41,13.85,13.77,14.07,14.79,15.31,13.53,12.65,12.1,12.06,
  25.85,25.53,22.62,23.59,22.69,19.79,20.05,21.04,20.13,21.78,19.56,17.35,16.77,17.49,17.9,
  22.33,21.88,20.05,22.24,25.17,21.29,22.34,22.68,22.97,22.21,19.44,18.38,17.5,17.48,16.97,
  13.33,12.48,11.24,11.55,11.8,11.07,11.84,11.62,11.77,11.98,12.53,11.3,11.31,11.37,11.44,
  15.61,15.2,13.34,13.82,13.93,12.87,12.87,13.36,13.24,13.2,12.77,11.29,12.17,11.43,11.33,
  7.78,8.77,9.68,10.12,10.32,10.95,10.39,10.42,10,9.91,7.34,7.47,8.52,9.29,10.31,
  8.85,8.08,7.86,7.93,8.75,7.91,8.75,8.74,8.32,8.84,7.99,6.98,7.45,8.23,8.66,
  12.6,12.26,11,11.35,12.82,11.66,11.61,11.79,12.41,13.54,12.55,9.59,10.35,11.39,11.72,
  24.31,23.3,20.14,19.88,20.38,19,19.26,18.78,18.46,18.97,18.54,17.33,17.24,16.3,16.28,
  12.13,11.26,10.78,11.46,11.5,10.69,10.82,10.88,10.52,10.35,11.13,10.01,10.22,10.41,10.26,
  14.23,11.96,11.91,10.39,12.34,10.43,11.83,12.89,12.9,13.63,14.34,12.57,13.85,11.02,10.93,
  12.45,11.9,12.22,11.67,12.83,11.67,11.86,11.79,11.72,12.85,11.38,12.24,10.01,11.17,10.87,
  23.76,22.83,20.57,23.34,20.92,18.31,18.71,19.37,19.22,19.33,17.15,15.68,15.69,16.01,15.58,
  12.19,11.17,10.28,10.5,10.93,10.43,10.86,11.07,11.06,11.44,12.58,12.03,11.82,11.66,11.26,
  14.72,14.24,13.31,13.32,13.28,12.75,13.08,13.39,13.87,14.6,14.83,13.41,13.56,13.17,12.81,
  16.21,15.02,14.97,13.39,13.57,12.59,13.33,12.69,12.64,13.04,11.64,11.31,11.16,11,11.16,
  9.35,8.67,8.13,8.68,9.01,8.16,9.78,9.7,9.59,8.61,11.02,10.89,10.6,10.15,10.29,
  11.23,11.69,11.09,12.98,13.97,13.38,16.69,16.6,19.21,16.87,14.55,13.62,15.63,16.59,14.99,
  15,14.11,13.07,13.69,13.73,12.46,12.95,13.06,12.46,13.12,11.7,10.76,11.73,10.98,10.85,
  11.29,11.08,10.64,11.3,11.18,12.11,9.48,10.48,11.39,12.31,9.39,9.26,9.4,9.42,9.37,
  19.46,19.66,16.81,27.53,26.98,16.88,16.53,17.75,17.83,16.09,13.25,12.65,12.51,14.18,13.61,
  10.01,9.28,9.04,9.55,9.79,10.6,8.76,8.78,8.98,9.67,10.19,9.17,9.63,9.94,10.06,
  10.75,9.5,9.5,10.04,10.03,9.46,8.74,8.67,8.46,8.82,9.56,9.37,8.61,9.01,9.37,
  15.64,14.47,12.55,12.98,13.17,12.32,12.61,13.95,14.62,13.83,14.38,13.61,14.26,14.41,13.91,
  9.08,8.53,8.2,8.98,9.39,9.92,10.2,10.25,9.63,10.57,10.33,7.65,7.93,8.42,8.44,
  16.4,15.84,14.52,14.64,14.97,14.03,14.1,13.25,14.36,13.8,13.72,12.89,12.88,13.03,13.14
), nrow=27, ncol=15, byrow=TRUE)

temp_state <- matrix(c(
  26,25,25,25,25,25,25,25,25,25,26,25,25,26,26,
  26,26,26,27,26,27,27,26,27,27,27,27,26,27,27,
  26,26,26,25,25,26,26,26,25,26,26,26,26,27,27,
  26,26,26,26,26,26,26,25,25,26,26,25,25,26,27,
  26,25,26,26,25,27,26,26,26,27,25,26,25,27,26,
  29,28,30,30,29,29,30,29,29,29,29,30,28,29,30,
  23,22,23,23,23,23,24,23,23,24,23,23,23,24,23,
  23,23,23,23,23,25,24,24,23,24,23,23,23,24,24,
  25,25,25,25,25,26,26,26,25,26,25,25,25,26,26,
  28,27,28,27,27,28,28,28,27,28,28,27,27,28,28,
  23,22,23,23,23,24,24,23,23,24,23,23,23,24,24,
  25,25,25,25,25,26,25,25,25,26,26,26,25,26,27,
  26,26,26,26,26,27,26,26,26,27,27,26,26,27,27,
  27,26,26,26,26,27,27,27,27,26,27,26,26,27,27,
  28,27,29,29,28,29,29,28,28,28,28,29,28,29,29,
  27,27,28,28,27,28,28,28,27,28,27,27,27,28,28,
  29,28,29,29,29,30,30,29,29,29,28,29,28,29,29,
  20,20,21,20,20,21,20,20,20,21,21,20,20,21,22,
  22,22,22,22,23,23,23,23,23,24,23,23,22,23,23,
  30,28,30,30,29,30,30,29,29,29,29,30,29,30,30,
  26,26,26,25,25,26,26,26,25,26,26,26,26,27,27,
  26,25,25,25,25,26,26,25,25,26,26,25,25,27,28,
  19,18,20,18,19,19,18,19,19,19,19,19,19,20,19,
  18,17,18,17,18,18,17,18,18,19,18,18,17,18,19,
  26,26,26,27,26,27,27,26,27,27,27,27,26,27,27,
  22,22,22,22,23,23,22,23,23,24,22,22,22,23,23,
  28,27,28,28,28,28,28,28,28,28,28,28,27,29,28
), nrow=27, ncol=15, byrow=TRUE)

prec_state <- matrix(c(
  1601,2024,2198,2139,2276,2137,1727,1707,2072,2122,2164,1826,1612,1707,1536,
  1153,997,477,943,1058,757,626,1159,712,771,967,973,1427,837,936,
  2388,2661,2681,2738,2731,2316,2322,2531,2620,2729,2411,2763,2264,2126,1981,
  2752,2628,2330,2578,2467,2492,2336,2507,2732,2427,2164,2829,2815,1973,2119,
  897,872,534,889,762,619,798,675,821,663,1079,1023,940,533,919,
  706,1159,447,628,703,626,641,744,840,930,976,678,833,667,757,
  1398,1494,1350,1516,1388,1191,1201,1171,1433,1134,1536,1437,1064,1044,1252,
  1261,1335,1077,1540,919,730,949,1034,1460,1114,1384,1427,1301,864,1189,
  1373,1595,1373,1615,1464,1350,1360,1358,1511,1227,1421,1471,1263,1190,1345,
  1344,1749,1148,1386,1469,1123,1216,1507,1624,1541,1750,1561,1540,1147,1334,
  1281,1428,1097,1431,919,1047,1227,1038,1333,1044,1472,1320,1523,1080,1384,
  1326,1475,1350,1466,1577,1366,1605,1503,1506,1222,1152,1321,1293,1346,1063,
  1629,1796,1667,1860,1851,1460,1644,1715,1796,1583,1465,1302,1502,1215,1345,
  2100,2393,2128,2326,2301,1804,2155,2199,2270,2310,2226,2248,2157,1627,1748,
  846,1257,422,740,838,624,647,745,756,819,938,565,975,654,716,
  884,974,343,817,779,558,589,719,673,659,812,812,1040,582,742,
  852,1134,670,984,899,712,864,897,1042,978,1115,983,963,682,805,
  1841,1850,1532,1888,2027,2411,1953,1993,1635,1587,1589,1400,1773,1896,1469,
  1574,1534,1344,1566,904,1113,949,1042,1383,1238,1589,1306,1411,1144,1215,
  692,1267,426,695,779,607,603,726,849,933,892,541,1074,681,805,
  1879,1839,1624,1600,1658,1366,1439,1692,1782,1722,1685,1349,1309,1168,1127,
  2361,2261,2137,2013,1666,2316,1901,1770,2145,2215,2054,2456,2531,1725,1932,
  1709,1717,1642,1761,2227,2350,1750,2070,1856,1912,1589,1444,1509,2187,2071,
  2175,2403,1566,2042,2351,2457,1781,1759,2732,2429,1578,1392,1744,2209,1879,
  1180,883,403,781,884,781,949,1016,658,850,1117,966,1119,809,823,
  1484,1578,1548,1522,1238,1688,1569,1578,1330,1266,1685,1063,1324,1377,1307,
  1417,1807,1479,1850,1767,1316,1430,1510,1705,1409,1690,1769,1652,1190,1470
), nrow=27, ncol=15, byrow=TRUE)

years_state <- 2010:2024

panel_list <- list()
for (i in seq_along(states)) {
  panel_list[[i]] <- data.frame(
    state    = states[i],
    region   = region_map[states[i]],
    year     = years_state,
    stunting = stunting_state[i, ],
    temp     = temp_state[i, ],
    prec     = prec_state[i, ]
  )
}
panel <- bind_rows(panel_list) %>%
  group_by(state) %>%
  mutate(stunting_z = scale(stunting)[,1]) %>%
  ungroup() %>%
  left_join(national %>% select(year, enso, ipca), by = "year")

climate_nat <- panel %>%
  group_by(year) %>%
  summarise(temp_mean     = mean(temp, na.rm=TRUE),
            prec_mean     = mean(prec, na.rm=TRUE),
            stunting_mean = mean(stunting, na.rm=TRUE),
            .groups="drop") %>%
  left_join(national %>% select(year, enso, enso_lag1), by="year")

cat("--- ENSO(t-1) → Temperature ---\n")
m_enso_temp <- lm(temp_mean ~ enso_lag1, data = climate_nat)
print(summary(m_enso_temp))

cat("\n--- ENSO(t-1) → Precipitation ---\n")
m_enso_prec <- lm(prec_mean ~ enso_lag1, data = climate_nat)
print(summary(m_enso_prec))

cat("\n--- Precipitation → Stunting ---\n")
m_prec_stunt <- lm(stunting_mean ~ prec_mean, data = climate_nat)
print(summary(m_prec_stunt))

cat("\n--- Temperature + Precipitation → Stunting ---\n")
m_clim_stunt <- lm(stunting_mean ~ temp_mean + prec_mean, data = climate_nat)
print(summary(m_clim_stunt))

# ============================================================
# PART IV — STATE PANEL: FIXED EFFECTS
# ============================================================
cat("\n============================================================\n")
cat("PART IV — STATE PANEL: FIXED EFFECTS (2010–2024)\n")
cat("============================================================\n\n")

panel_plm <- pdata.frame(panel, index = c("state", "year"))

cat("--- FE Model: Stunting ~ Temp + Prec ---\n")
fe1 <- plm(stunting ~ temp + prec,
           data = panel_plm, model = "within", effect = "individual")
print(summary(fe1))

cat("\n--- FE Model (z-score): Stunting_z ~ Temp + Prec ---\n")
fe3 <- plm(stunting_z ~ temp + prec,
           data = panel_plm, model = "within", effect = "individual")
print(summary(fe3))

cat("\n--- Hausman Test (FE vs RE) ---\n")
fe2 <- plm(stunting ~ temp + prec + enso,
           data = panel_plm, model = "within", effect = "individual")
re1 <- plm(stunting ~ temp + prec + enso,
           data = panel_plm, model = "random")
print(phtest(fe2, re1))

# ============================================================
# PART V — REGIONAL SUMMARY
# ============================================================
cat("\n============================================================\n")
cat("PART V — REGIONAL SUMMARY\n")
cat("============================================================\n\n")

regional_summary <- panel %>%
  group_by(region) %>%
  summarise(
    stunting_mean  = round(mean(stunting, na.rm=TRUE), 2),
    stunting_trend = round(coef(lm(stunting ~ year))[2], 3),
    temp_mean      = round(mean(temp, na.rm=TRUE), 1),
    prec_mean      = round(mean(prec, na.rm=TRUE), 0),
    .groups = "drop"
  ) %>%
  arrange(desc(stunting_mean))

cat("Regional averages (2010–2024):\n")
print(regional_summary)

cat("\n=== Analysis complete ===\n")

#Figure 1 — ENSO variability and annual Food CPI variation (%), Brazil (2001–2024)#
library(tidyverse)

# Define ipca_raw explicitamente
ipca_raw <- c(32.18,35.27,42.14,45.29,47.04,47.97,48.57,53.80,59.78,61.67,
              68.08,72.98,80.18,86.97,93.96,105.24,114.30,112.17,116.69,124.14,
              141.63,152.89,170.70,172.40,185.70)

# Calculate annual variation inside national dataframe
national <- national %>%
  mutate(foodcpi_var = c(NA, diff(ipca) / ipca[-25] * 100))

# Create graphic dataframe
df_fig1 <- national %>%
  filter(!is.na(foodcpi_var)) %>%
  mutate(category = case_when(
    enso >= 1.0  ~ "Strong El Niño",
    enso >= 0.5  ~ "Weak El Niño",
    enso <= -1.0 ~ "Strong La Niña",
    enso <= -0.5 ~ "Weak La Niña",
    TRUE ~ "Neutral"
  ),
  category = factor(category, levels = c(
    "Strong El Niño","Weak El Niño","Neutral","Weak La Niña","Strong La Niña"
  )))

# Gráfico
ggplot(df_fig1, aes(x = year)) +
  geom_rect(aes(xmin=year-0.5, xmax=year+0.5,
                ymin=-Inf, ymax=Inf, fill=category), alpha=0.3) +
  geom_line(aes(y = enso * 3), color="#378ADD", linewidth=1, linetype="dashed") +
  geom_point(aes(y = enso * 3), color="#378ADD", size=2) +
  geom_line(aes(y = foodcpi_var), color="#C00000", linewidth=1.3) +
  geom_point(aes(y = foodcpi_var), color="#C00000", size=2.5) +
  scale_fill_manual(values = c(
    "Strong El Niño" = "#B41E1E",
    "Weak El Niño"   = "#EF8C14",
    "Strong La Niña" = "#0F3C82",
    "Weak La Niña"   = "#378ADD",
    "Neutral"        = "transparent"
  )) +
  scale_y_continuous(
    name = "Food CPI annual variation (%)",
    sec.axis = sec_axis(~./3, name = "ENSO (Niño 3.4 index)")
  ) +
  labs(
    title = "ENSO variability and Food CPI annual variation, Brazil (2001–2024)",
    subtitle = "Background = ENSO episode. Red line = Food CPI variation (left). Blue dashed = ENSO (right).",
    x = NULL, fill = "ENSO category"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("fig1_enso_foodcpi_v2.pdf", width = 10, height = 5, dpi = 300)

#Figure 2 — ENSO episode effects on Food CPI — dummy model#
library(tidyverse)

# Create coef_df
coef_df <- data.frame(
  label = c("Weak La Niña (t)","Neutral (t)","Weak El Niño (t)","Strong El Niño (t)",
            "Weak La Niña (t−1)","Neutral (t−1)","Weak El Niño (t−1)","Strong El Niño (t−1)"),
  beta  = c(25.33, 15.32, 28.00, 25.21, 22.47, 11.41, 26.56, 31.87),
  se    = c(11.88, 10.88, 13.21, 12.66, 10.25, 10.33, 11.89, 13.11),
  pval  = c(0.051, 0.181, 0.052, 0.066, 0.046, 0.288, 0.042, 0.029)
) %>%
  mutate(
    ci_lo  = beta - 1.96 * se,
    ci_hi  = beta + 1.96 * se,
    sig    = case_when(pval < 0.05 ~ "p < 0.05", pval < 0.10 ~ "p < 0.10", TRUE ~ "n.s."),
    mark   = case_when(pval < 0.05 ~ "*", pval < 0.10 ~ "⚠", TRUE ~ ""),
    label2 = sprintf("β=%.2f%s\n[%.1f,%.1f]\np=%.3f", beta, mark, ci_lo, ci_hi, pval),
    label  = factor(label, levels = rev(label))
  )

# Graphic
ggplot(coef_df, aes(x = label, y = beta, fill = sig, color = sig)) +
  geom_col(alpha = 0.2, width = 0.55) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(y = ci_hi + 2, label = label2),
            vjust = 0, size = 2.6, lineheight = 1.3, fontface = "plain") +
  scale_fill_manual(values = c("p < 0.05"="#E24B4A","p < 0.10"="#EF9F27","n.s."="#B4B2A9")) +
  scale_color_manual(values = c("p < 0.05"="#E24B4A","p < 0.10"="#EF9F27","n.s."="#B4B2A9")) +
  scale_y_continuous(limits = c(-15, 75)) +
  labs(
    title    = "ENSO episode effects on Food CPI — dummy model (Eq. 14)",
    subtitle = "Coefficients vs. strong La Niña reference. Error bars = 95% CI. Adj. R² = 0.964.",
    x = NULL, y = "β (Food CPI points)", color = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9, lineheight = 1.2))

ggsave("fig2_dummy_foodcpi_v3.pdf", width = 10, height = 6, dpi = 300)


#Figure 4 — Mean childhood stunting prevalence by state, Brazil (2010–2024)#
install.packages("geobr")
library(tidyverse)
library(geobr)      # install.packages("geobr") 
library(ggplot2)

states_geo <- geobr::read_state(year = 2020, showProgress = FALSE)

stunting_state_means <- data.frame(
  abbrev_state  = c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA",
                    "MG","MS","MT","PA","PB","PE","PI","PR","RJ","RN",
                    "RO","RR","RS","SC","SE","SP","TO"),
  stunting_mean = c(19.6,13.7,20.5,20.2,11.7,13.0,9.5,8.3,11.7,18.9,
                    10.8,12.4,11.7,18.4,11.4,13.6,12.8,9.7,14.6,12.6,
                    10.5,17.4,9.7,9.4,13.7,9.3,14.0)
)

states_geo_data <- states_geo %>%
  left_join(stunting_state_means, by = "abbrev_state") %>%
  mutate(label_pct = paste0(round(stunting_mean, 1), "%"))

centroids <- states_geo_data %>%
  sf::st_centroid() %>%
  mutate(
    lon = sf::st_coordinates(.)[,1],
    lat = sf::st_coordinates(.)[,2]
  ) %>%
  sf::st_drop_geometry()

ggplot(states_geo_data) +
  
  geom_sf(aes(fill = stunting_mean), color = "white", linewidth = 0.3) +
  
  geom_sf(data = states_geo_data %>% summarise(),
          fill = NA, color = "white", linewidth = 0.9) +
  
  geom_text(data     = centroids,
            aes(x = lon, y = lat + 0.25, label = abbrev_state),
            size     = 1.9,
            color    = "gray20",
            fontface = "plain",
            alpha    = 0.75) +
  
  geom_text(data     = centroids,
            aes(x = lon, y = lat - 0.25, label = label_pct),
            size     = 2.1,
            color    = "black",
            fontface = "bold") +
  
  # Gradient color
  scale_fill_gradient(
    low    = "#FCBBA1",   # vermelho muito claro
    high   = "#67000D",   # vermelho muito escuro
    name   = "Stunting (%)"
  ) +
  
  labs(
    title    = "Mean childhood stunting prevalence by state, Brazil (2010–2024)",
    subtitle = "% children under 5 with HAZ < −2 SD. Source: SISVAN.",
    caption  = "Light red = low prevalence; Dark red = high prevalence. State boundaries: IBGE/geobr (2020)."
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(size = 12, face = "bold"),
    plot.subtitle    = element_text(size = 10),
    plot.caption     = element_text(size = 8, color = "gray50"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("fig4_stunting_map_v5.png",
       width  = 8,
       height = 7,
       dpi    = 300,
       bg     = "white")



#Figure 3 - Local polynomial regression

install.packages("patchwork")
library(tidyverse)
library(patchwork)
library(gridExtra)

# Translate region for english - panel_data
panel_data <- panel_data %>%
  mutate(regiao = recode(regiao,
                         'Norte'        = 'North',
                         'Nordeste'     = 'Northeast',
                         'Centro-Oeste' = 'Central-West',
                         'Sudeste'      = 'Southeast',
                         'Sul'          = 'South'),
         regiao = factor(regiao,
                         levels = c('North','Northeast','Central-West',
                                    'Southeast','South')))

# Region color
cores_regiao <- c(
  "North"        = "#B41E1E",
  "Northeast"    = "#EF8C14",
  "Central-West" = "#2E7D32",
  "Southeast"    = "#1565C0",
  "South"        = "#6A1B9A"
)

# ── FIGURE 3A: Stunting × Temperature ─────────────────────────
fig_temp <- ggplot(panel_data,
                   aes(x = temp, y = stunting,
                       color = regiao, fill = regiao)) +
  geom_point(alpha = 0.20, size = 1.5) +
  geom_smooth(method    = "loess",
              span      = 0.85,
              se        = TRUE,
              alpha     = 0.15,
              linewidth = 1.2) +
  scale_color_manual(values = cores_regiao) +
  scale_fill_manual(values  = cores_regiao) +
  scale_y_continuous(limits = c(5, 30)) +
  labs(
    title = "A — Stunting prevalence by mean temperature",
    x     = "Mean annual temperature (°C)",
    y     = "Child stunting prevalence (%)",
    color = "Region",
    fill  = "Region"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 11, face = "bold")
  )

# ── FIGURE 3B: Stunting × Precipitation ───────────────────────
fig_prec <- ggplot(panel_data,
                   aes(x = prec, y = stunting,
                       color = regiao, fill = regiao)) +
  geom_point(alpha = 0.20, size = 1.5) +
  geom_smooth(method    = "loess",
              span      = 0.85,
              se        = TRUE,
              alpha     = 0.15,
              linewidth = 1.2) +
  scale_color_manual(values = cores_regiao) +
  scale_fill_manual(values  = cores_regiao) +
  scale_y_continuous(limits = c(5, 30)) +
  labs(
    title = "B — Stunting prevalence by annual precipitation",
    x     = "Mean annual precipitation (mm)",
    y     = "Child stunting prevalence (%)",
    color = "Region",
    fill  = "Region"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(size = 11, face = "bold")
  )

# ── SAVE FIGURE 3 ──────────────────────────────────
ggsave("fig3A_lpoly_stunting_temp.pdf",
       plot = fig_temp, width = 7, height = 5, dpi = 300)

ggsave("fig3B_lpoly_stunting_prec.pdf",
       plot = fig_prec, width = 7, height = 5, dpi = 300)

# ── SAVE FIGURE - MERGED ───────────────────────────────────
pdf("fig3_stunting_climate_combined.pdf", width = 14, height = 6)
grid.arrange(fig_temp, fig_prec, ncol = 2)
dev.off()

png("fig3_stunting_climate_combined.png", width = 14, height = 6,
    units = "in", res = 300)
grid.arrange(fig_temp, fig_prec, ncol = 2)
dev.off()

cat("Figuras salvas!\n")
