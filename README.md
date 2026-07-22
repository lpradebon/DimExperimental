
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DimExp

<img src="man/figures/logo.png" align="right" height="140" />

<!-- badges: start -->

<!-- badges: end -->

DimExp provides tools for experimental design sizing in agricultural
research. It implements several methods for estimating the optimal plot
size from uniformity trial data, and for estimating the optimal number
of replications, with standardized diagnostic statistics and plots:

- `fit_lrp()`: optimal plot size via the Linear Response Plateau model
- `fit_qrp()`: optimal plot size via the Quadratic Response Plateau
  model
- `fit_mcm()`: optimal plot size via the Modified Maximum Curvature
  method (Meier & Lessmann, 1971)
- `calc_paranaiba()`: optimal plot size via the Paranaiba method
  (Paranaiba, Ferreira & Morais, 2009)
- `calc_replicates()`: optimal number of replications (Cargnelutti Filho
  et al., 2014)

## Installation

You can install the development version of DimExp from
[GitHub](https://github.com/lpradebon/DimExp) with:

``` r
# install.packages("pak")
pak::pak("lpradebon/DimExp")
```

## Example

A basic uniformity-trial workflow: fit a Linear Response Plateau model
to coefficient of variation data across plot sizes, and inspect the
estimated optimal plot size. The breakpoint is found by a grid search,
so no starting values are needed.

``` r
library(DimExp)

# Example uniformity trial data: CV (%) decreasing with plot size
plot_size <- c(1, 2, 4, 6, 8, 10, 15, 20, 25, 30)
cv        <- c(22.1, 18.4, 14.2, 11.8, 10.1, 9.3, 8.1, 7.6, 7.4, 7.3)

fit <- fit_lrp(x = plot_size, cv = cv)
fit
#> Linear Response Plateau (LRP) fit
#> Method:                  segment 
#> Breakpoint (Xo):         8.648 
#> CV at breakpoint:        7.940 
#> R2: 0.963  RMSE: 0.935  AIC: 35.0  BIC: 36.3
```

The plot title (and other options) go on the `plot()` method, which
draws the fitted broken line, the breakpoint, and the plateau-model
annotations in the layout commonly used in plot-size articles:

``` r
plot(fit, title = "Uniformity trial example")
```

<img src="man/figures/README-example-plot-1.png" alt="" width="100%" />

### Several trials at once

Pass a data frame and the column names to fit one model per trial. The
result carries a per-trial summary table plus the individual fits:

``` r
trials <- rbind(
  data.frame(plot_size, cv, trial = "Trial A"),
  data.frame(plot_size, cv = cv + 1.5, trial = "Trial B")
)

res <- fit_lrp(trials, x = "plot_size", cv = "cv", trial = "trial")
#> Using x = 'plot_size', cv = 'cv', trial = 'trial' -> 2 trials.
res$summary
#>     trial       a       b breakpoint plateau     R2   RMSE    AIC    BIC
#> 1 Trial A 22.2884 -1.6591      8.648  7.9401 0.9628 0.9354 35.043 36.253
#> 2 Trial B 23.7884 -1.6591      8.648  9.4401 0.9628 0.9354 35.043 36.253
```

If the optimum is known to lie in a given region, restrict the
breakpoint search with `search_range`, for example
`fit_lrp(x = plot_size, cv = cv, search_range = c(5, 20))`.

## Citation

If you use DimExp in your research, please cite the underlying methods
(see each function’s documentation for the corresponding reference) as
well as the package itself.
