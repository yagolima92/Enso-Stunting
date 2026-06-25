# ENSO and Childhood Stunting in Brazil: Climate Effects Mediated by Economic and Nutritional Pathways

This repository contains the R code used for all statistical analyses reported in the study:

> Lima, A. et al. *ENSO and Childhood Stunting in Brazil: Climate Effects Mediated by Economic and Nutritional Pathways*. IDEAS 2026 — International Database Engineered Applications Symposium. Springer LNCS (2026).

---

## Repository Contents

- `CODE.md`: Complete R code for all analyses reported in the article, organized by analytical part.
- `README.md`: This file.

---

## Study Overview

This study investigates the relationship between El Niño–Southern Oscillation (ENSO) variability and childhood stunting in Brazil (2000–2025), using national time-series and state-level panel data for 27 states (2010–2024). Four model specifications following Hsiang and Meng (2015) reveal that lagged ENSO episodes are significantly associated with higher Food CPI levels, which in turn reduce animal protein consumption and increase childhood stunting prevalence. The study identifies an indirect climate–nutrition pathway: ENSO instability → inflation → reduced dietary quality → food insecurity → increased child stunting.

---

## Analytical Structure

The R code is organized into five parts:

### Part I — Structural National Mechanisms (2000–2025)
- Pearson correlation coefficients between detrended series (PoU, maternal anemia, meat consumption, Food CPI, stunting)
- Bivariate OLS models (Equations 6–8)
- Multivariate OLS model (Equation 9)

### Part II — ENSO → Food CPI: Four Model Specifications
- Specification 1: Linear model (Equation 10)
- Specification 2: Cubic polynomial model (Equation 11)
- Specification 3: Generalized Additive Model — GAM (Equation 12)
- Specification 4: Dummy model — preferred (Equation 13)

### Part III — Climatic Pathway Models (2010–2024)
- ENSO → Temperature (Equation 14)
- ENSO → Precipitation (Equation 15)
- Stunting ~ Temperature + Precipitation (Equation 16)

### Part IV — State-Level Fixed Effects Panel (2010–2024)
- Fixed effects model (Equation 17)
- Hausman test (Equation 18)
- Regional OLS models with HC3 robust standard errors (Equation 19)
- F-test for regional coefficient heterogeneity (Equation 20)

### Part V — Figures and Visualization
- Figure 1: ENSO variability and annual Food CPI variation (2001–2024)
- Figure 2: ENSO episode effects on Food CPI — dummy model coefficients
- Figure 3: LOESS regression — stunting by temperature and precipitation across Brazilian regions
- Figure 4: Choropleth map — mean childhood stunting prevalence by state (2010–2024)

---

## Data Sources

All data used in this study are publicly available. No ethical approval was required.

| Variable | Source | Coverage | Scale |
|---|---|---|---|
| Childhood stunting prevalence | SISVAN / FAOSTAT | 2000–2025 | National / State |
| Prevalence of Undernourishment (PoU) | FAOSTAT | 2000–2025 | National |
| Maternal anemia prevalence | FAOSTAT | 2000–2025 | National |
| Animal protein consumption | FAOSTAT | 2000–2025 | National |
| ENSO — Niño 3.4 index | NOAA | 2000–2025 | Global |
| Food CPI | Banco Central do Brasil | 2000–2025 | National |
| Mean surface temperature (°C) | MapBiomas | 2010–2024 | State |
| Annual precipitation (mm) | MapBiomas | 2010–2024 | State |

**Data access links:**
- SISVAN: https://sisvan.saude.gov.br
- FAOSTAT: https://www.fao.org/faostat
- NOAA Niño 3.4 index: https://www.cpc.ncep.noaa.gov
- Banco Central do Brasil: https://www.bcb.gov.br
- MapBiomas: https://mapbiomas.org

---

## Required R Packages

```r
install.packages(c(
  "tidyverse",   # data manipulation and visualization
  "ggplot2",     # figures
  "plm",         # panel data models and Hausman test
  "lmtest",      # coefficient tests
  "sandwich",    # HC3 robust standard errors
  "mgcv",        # generalized additive models (GAM)
  "geobr",       # Brazilian state boundaries (Fig. 4)
  "sf",          # spatial features
  "gridExtra",   # panel figure layout
  "patchwork"    # optional: alternative panel layout
))
```

---

## How to Run

1. Clone or download this repository.
2. Obtain the data from the sources listed above and organize them as described in the article (Table 1).
3. Open `CODE.md` and run each section sequentially in R or RStudio.
4. All figures will be saved as `.pdf` and `.png` files in your working directory.

> **Note:** The code uses simulated panel data based on reported state-level means for the regional and LOESS analyses (Parts IV and V), since the full anonymized state × year dataset is not publicly available in a single file. The national time-series data (Parts I–III) can be fully reproduced from the public sources listed above.

---

## Software

All analyses were performed in **R (version 4.x)**. The code was developed and tested in RStudio.

```
R version 4.x.x
Platform: x86_64-pc-linux-gnu
```

---

## Citation

If you use this code, please cite the original article:

> Lima, A. et al. ENSO and Childhood Stunting in Brazil: Climate Effects Mediated by Economic and Nutritional Pathways. In: *Proceedings of IDEAS 2026*. Lecture Notes in Computer Science. Springer, Cham (2026).


---

## Contact

For questions or issues regarding this code, please contact:

**Alan Lima**
Email: alanlima.agro@hotmail.com
