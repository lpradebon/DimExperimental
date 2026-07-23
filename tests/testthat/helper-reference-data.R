## Reference data shared by the regression tests.
## Files named helper-*.R are sourced by testthat before the tests run.

## ---------------------------------------------------------------------------
## Chickpea uniformity trials, Cargnelutti Filho et al. (2025),
## Revista Vivencias 21(43), 499-513, Table 1.
## 6 x 6 grid of 1 m2 BEU per trial, grouped into 15 plot shapes.
## ---------------------------------------------------------------------------
chickpea_X <- c(1, 2, 2, 3, 3, 4, 6, 6, 6, 6, 9, 12, 12, 18, 18)

## number of plots of each size (n = 36 / X); degrees of freedom = n - 1
chickpea_n <- c(36, 18, 18, 12, 12, 9, 6, 6, 6, 6, 4, 3, 3, 2, 2)

chickpea_CV1 <- c(30.40, 19.51, 23.72, 12.89, 21.32, 16.69, 6.71, 10.75,
                  17.58, 14.94, 11.93, 3.18, 8.63, 4.25, 11.41)
chickpea_CV2 <- c(34.79, 29.44, 26.45, 26.35, 31.74, 25.01, 21.44, 22.86,
                  28.93, 16.17, 27.29, 19.89, 16.85, 26.53, 14.70)
chickpea_CV3 <- c(31.12, 20.87, 29.09, 20.02, 28.04, 19.97, 10.74, 19.28,
                  18.18, 26.46, 17.41, 11.51, 17.96, 9.18, 19.17)

chickpea_cv_list <- list(chickpea_CV1, chickpea_CV2, chickpea_CV3)

chickpea_long <- rbind(
  data.frame(x = chickpea_X, cv = chickpea_CV1, trial = "Trial 1"),
  data.frame(x = chickpea_X, cv = chickpea_CV2, trial = "Trial 2"),
  data.frame(x = chickpea_X, cv = chickpea_CV3, trial = "Trial 3"),
  stringsAsFactors = FALSE
)

## Published optima (Table 2 of the article), in m2
chickpea_published <- list(
  LRP = list(Xo = c(7.44, 6.56, 7.57), CVxo = c(7.88, 21.05, 15.05),
             mean_Xo = 7.19),
  QRP = list(Xo = c(10.97, 8.46, 11.31), CVxo = c(7.50, 21.06, 14.76),
             mean_Xo = 10.25),
  MCM = list(Xo = c(5.76, 3.99, 4.68), CVxo = c(12.54, 25.64, 20.12),
             mean_Xo = 4.81)
)

## ---------------------------------------------------------------------------
## Raw grid of chickpea trial 1 (Figure 2 of the article), g m-2.
## Used by the Paranaiba tests, which need the field layout.
## ---------------------------------------------------------------------------
chickpea_grid1 <- matrix(
  c(288.35, 106.10, 213.70,  92.88, 185.88, 310.72,
    262.45, 290.88, 177.45, 162.66, 346.70, 205.90,
    303.74, 196.98, 196.11, 185.08, 137.66, 311.93,
    339.32, 234.04, 273.59, 228.91, 112.36, 215.99,
    311.05, 228.65, 227.15, 259.83, 293.07, 111.10,
    238.23, 204.57, 336.62, 242.58, 173.34, 191.17),
  nrow = 6, byrow = TRUE)

chickpea_grid2 <- matrix(
  c(274.68, 244.44, 296.02, 182.82, 170.22, 270.54,
    272.81, 207.95, 217.12, 307.88, 235.11, 263.63,
    303.56, 184.47, 147.03, 255.88, 210.06, 227.52,
    190.72, 421.05, 528.89, 351.25, 310.37, 156.68,
    294.45, 451.29, 496.12, 446.41, 302.98, 248.26,
    437.76, 434.10, 402.28, 329.27, 278.65, 164.26),
  nrow = 6, byrow = TRUE)

## ---------------------------------------------------------------------------
## Black oat, Cargnelutti Filho et al. (2014), Ciencia Rural 44(10), 1732-1739.
## Table 1 (first evaluation date): rho, variance and mean of the nine trials,
## with the published Xo and CVxo. Used to validate the Paranaiba formulas.
## ---------------------------------------------------------------------------
oat_rho  <- c(0.47, 0.57, 0.49, 0.21, -0.01, 0.14, 0.34, 0.71, 0.65)
oat_s2   <- c(200070.61, 149109.06, 183757.81, 222511.10, 188283.27,
              104665.90, 251591.19, 319318.12, 390063.04)
oat_mean <- c(2165.34, 2162.05, 2297.89, 2272.59, 2167.48,
              2276.14, 2126.17, 1985.56, 2030.86)
oat_Xo_published   <- c(4.05, 3.51, 3.75, 4.35, 4.31, 3.41, 4.62, 4.32, 4.79)
oat_CVxo_published <- c(9.05, 7.84, 8.39, 9.73, 9.64, 7.62, 10.33, 9.66, 10.72)

## CVxo used in Tables 2 and 3 of the same article
oat_CVxo <- 9.25
