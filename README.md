
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
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("lpradebon/DimExp")
```

## Example

A basic uniformity-trial workflow: fit a Linear Response Plateau model
to coefficient of variation data across plot sizes, and inspect the
estimated optimal plot size.

``` r
library(DimExp)

# Example uniformity trial data: CV (%) decreasing with plot size
plot_size <- c(1, 2, 4, 6, 8, 10, 15, 20, 25, 30)
cv        <- c(22.1, 18.4, 14.2, 11.8, 10.1, 9.3, 8.1, 7.6, 7.4, 7.3)

fit <- fit_lrp(x = plot_size, cv = cv, plot_title = "Uniformity trial example")
fit
#> Linear Response Plateau (LRP) fit
#> Breakpoint (Xo):         10.059 
#> CV at breakpoint:        7.600 
#> R2: 0.948  RMSE: 1.108  AIC: 38.4  BIC: 39.6
```

``` r
plot(fit)
#> Warning in ggplot2::geom_segment(ggplot2::aes(x = 0, xend = breakpoint, : All aesthetics have length 1, but the data has 10 rows.
#> ℹ Please consider using `annotate()` or provide this layer with data containing
#>   a single row.
#> Warning in ggplot2::geom_segment(ggplot2::aes(x = breakpoint, xend = breakpoint, : All aesthetics have length 1, but the data has 10 rows.
#> ℹ Please consider using `annotate()` or provide this layer with data containing
#>   a single row.
#> Warning in ggplot2::geom_point(ggplot2::aes(x = breakpoint, y = breakpoint_response), : All aesthetics have length 1, but the data has 10 rows.
#> ℹ Please consider using `annotate()` or provide this layer with data containing
#>   a single row.
#> Warning: Removed 10 rows containing missing values or values outside the scale range
#> (`geom_segment()`).
```

<img src="man/figures/README-example-plot-1.png" alt="" width="100%" />

## Citation

If you use DimExp in your research, please cite the underlying methods
(see each function’s documentation for the corresponding reference) as
well as the package itself.
