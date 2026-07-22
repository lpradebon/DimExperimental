## ============================================================================
## DimExp :: Linear Response Plateau (LRP) model by grid search
## ============================================================================

## ----------------------------------------------------------------------------
## Internal single-series engine (not exported)
## ----------------------------------------------------------------------------
#' Fit one LRP series (internal)
#'
#' Core grid-search fitter for a single \code{x}/\code{cv} series. Users should
#' call [fit_lrp()] instead.
#'
#' @param x,cv numeric vectors of equal length.
#' @param step grid step for the breakpoint search.
#' @param method \code{"segment"} or \code{"ramp"}.
#' @return An object of class \code{"lrp_fit"}.
#' @keywords internal
#' @noRd
.lrp_fit_one <- function(x, cv, step = 0.001, method = "segment",
                         search_range = NULL) {

  if (!is.numeric(x) || !is.numeric(cv))
    stop("`x` and `cv` must be numeric.", call. = FALSE)
  if (length(x) != length(cv))
    stop("`x` and `cv` must have the same length.", call. = FALSE)
  if (anyNA(x) || anyNA(cv))
    stop("`x` and `cv` cannot contain missing values (NA).", call. = FALSE)
  if (length(x) < 4)
    stop("At least 4 observations are required.", call. = FALSE)
  if (length(unique(x)) < 3)
    stop("At least 3 distinct `x` values are required for the breakpoint grid.",
         call. = FALSE)

  x_unique <- sort(unique(x))
  feas_lo  <- x_unique[2]
  feas_hi  <- x_unique[length(x_unique)]

  ## optional restriction of the breakpoint search range
  if (is.null(search_range)) {
    grid_lo <- feas_lo; grid_hi <- feas_hi
  } else {
    if (!is.numeric(search_range) || length(search_range) != 2 ||
        anyNA(search_range) || search_range[1] >= search_range[2])
      stop("`search_range` must be numeric c(lower, upper) with lower < upper.",
           call. = FALSE)
    if (search_range[1] < min(x) || search_range[2] > max(x))
      stop(sprintf("`search_range` must fall within the data range [%.4g, %.4g].",
                   min(x), max(x)), call. = FALSE)
    grid_lo <- max(search_range[1], feas_lo)
    grid_hi <- search_range[2]
    if (grid_lo > grid_hi)
      stop(sprintf(paste0("`search_range` upper bound is below the feasible ",
                          "breakpoint minimum (%.4g); at least 2 points must ",
                          "lie below the breakpoint."), feas_lo), call. = FALSE)
  }

  ## residual sum of squares for a fixed breakpoint
  ss_at <- function(x0) {
    if (method == "segment") {
      m <- x <= x0
      if (sum(m) < 2) return(list(ss = Inf))
      cf <- stats::coef(stats::lm(cv[m] ~ x[m]))
    } else {
      z  <- pmin(x, x0)
      cf <- stats::coef(stats::lm(cv ~ z))
    }
    a <- unname(cf[1]); b <- unname(cf[2]); p <- a + b * x0
    fit <- ifelse(x <= x0, a + b * x, p)
    list(ss = sum((cv - fit)^2), a = a, b = b, p = p)
  }

  ## grid search over the breakpoint
  grid <- seq(grid_lo, grid_hi, by = step)
  best <- list(ss = Inf)
  for (x0 in grid) {
    cand <- ss_at(x0)
    if (cand$ss < best$ss) best <- c(cand, list(x0 = x0))
  }
  if (!is.finite(best$ss))
    stop("Grid search failed to find a valid breakpoint.", call. = FALSE)

  a <- best$a; b <- best$b; x0 <- best$x0; p <- best$p

  ## fitted values and fit statistics
  fitted    <- ifelse(x <= x0, a + b * x, p)
  residuals <- cv - fitted
  n   <- length(cv)
  rss <- sum(residuals^2)
  r_squared <- 1 - rss / sum((cv - mean(cv))^2)
  rmse      <- sqrt(mean(residuals^2))

  ## AIC / BIC via the Gaussian log-likelihood (MLE variance = rss / n).
  ## Parameter count k = a, b, Xo, sigma = 4, matching the nls / AIC.default
  ## convention. Only comparable with a QRP fit that uses the same count.
  loglik <- -0.5 * n * (log(2 * pi) + log(rss / n) + 1)
  k <- 4L
  aic <- -2 * loglik + 2 * k
  bic <- -2 * loglik + log(n) * k

  if (b >= 0)
    warning("Estimated slope is non-negative; the decreasing-plateau ",
            "interpretation may not hold for these data.", call. = FALSE)
  if (isTRUE(all.equal(x0, grid_lo)) || isTRUE(all.equal(x0, grid_hi)))
    warning("Breakpoint is at the edge of the search range; the data (or the ",
            "supplied `search_range`) may not bracket the true breakpoint.",
            call. = FALSE)

  structure(
    list(
      coefficients = c(a = a, b = b),
      parameters   = c(Breakpoint = x0, Breakpoint_Response = p,
                       R2 = r_squared, RMSE = rmse, AIC = aic, BIC = bic),
      fitted = fitted, residuals = residuals,
      data = data.frame(x = x, cv = cv), method = method, step = step,
      search_range = search_range
    ),
    class = "lrp_fit"
  )
}

## ----------------------------------------------------------------------------
## Public fitter
## ----------------------------------------------------------------------------
#' Fit the Linear Response Plateau (LRP) model by grid search
#'
#' Fits the linear-plateau (broken-line) model
#' \deqn{f(x) = a + b\,x \ \text{ if } x \le X_0, \qquad
#'       f(x) = a + b\,X_0 \ \text{ if } x > X_0}
#' by profiling the breakpoint \eqn{X_0} over a fine grid. For each candidate
#' breakpoint the linear coefficients are obtained by least squares and the
#' plateau is set to \eqn{a + b\,X_0}; the breakpoint minimizing the residual
#' sum of squares over all observations is returned. Profiling the breakpoint
#' avoids the starting-value sensitivity and local-minima of a direct
#' \code{\link[stats]{nls}} fit, so no initial values are required.
#'
#' The response is typically the coefficient of variation (CV, percent) and the
#' predictor the plot size. Supply the individual CV values (one per basic-unit
#' form), not the means per plot size; repeated \code{x} values are expected.
#'
#' @section Two ways to call:
#' \describe{
#'   \item{Vectors}{\code{fit_lrp(x, cv)} fits a single series and returns an
#'     \code{"lrp_fit"}.}
#'   \item{Data frame}{\code{fit_lrp(.data, x = "x", cv = "cv", trial = "trial")}
#'     takes a data frame plus column names. With \code{trial}, one model is fit
#'     per trial and an \code{"lrp_multi"} object (summary table plus the
#'     individual fits) is returned. Column names default to \code{"x"},
#'     \code{"cv"} and \code{"trial"}; missing columns raise a clear error.}
#' }
#'
#' @param .data optional data frame. When supplied, \code{x}, \code{cv} and
#'   \code{trial} are interpreted as column names (character strings). When
#'   \code{NULL} (default) the vector interface is used.
#' @param x,cv either numeric vectors (vector interface) or, when \code{.data}
#'   is a data frame, the names of the predictor and response columns.
#' @param trial optional name of the column identifying the trial. When given,
#'   one model is fit per trial.
#' @param step grid step for the breakpoint search. Default \code{0.001}; larger
#'   values run faster with a slightly coarser breakpoint.
#' @param method how the linear coefficients are estimated at each candidate
#'   breakpoint. \code{"segment"} (default) uses only observations with
#'   \code{x <= X0}, reproducing the Paranaiba et al. (2009) procedure.
#'   \code{"ramp"} uses all observations on the basis \code{pmin(x, X0)}, the
#'   standard least-squares LRP.
#' @param search_range optional numeric \code{c(lower, upper)} restricting the
#'   interval (in units of \code{x}) where the breakpoint is searched. Must fall
#'   within the data range. Useful when the optimum is known to lie in a region
#'   and an outlier could otherwise pull the breakpoint outside it. \code{NULL}
#'   (default) searches the full feasible interval.
#'
#' @return
#' For a single series, an object of class \code{"lrp_fit"}: a list with
#' \code{coefficients} (a, b); \code{parameters} (Breakpoint, Breakpoint_Response,
#' R2, RMSE, AIC, BIC); \code{fitted}; \code{residuals}; \code{data};
#' \code{method}; \code{step}. \cr
#' With \code{trial}, an object of class \code{"lrp_multi"}: a list with
#' \code{summary} (one row per trial), \code{fits} (the individual
#' \code{"lrp_fit"} objects) and \code{method}.
#'
#' @references
#' Paranaiba, P. F., Ferreira, D. F. & Morais, A. R. (2009). Tamanho otimo de
#' parcelas experimentais: proposicao de metodos de estimacao.
#' \emph{Revista Brasileira de Biometria}, 27(2), 255-268. \cr
#' Cargnelutti Filho, A. et al. (2025). \emph{Revista Vivencias}, 21(43), 499-513.
#'
#' @examples
#' X   <- c(1, 2, 2, 3, 3, 4, 6, 6, 6, 6, 9, 12, 12, 18, 18)
#' CV1 <- c(30.40, 19.51, 23.72, 12.89, 21.32, 16.69, 6.71, 10.75,
#'          17.58, 14.94, 11.93, 3.18, 8.63, 4.25, 11.41)
#'
#' ## Vector interface (single series)
#' fit <- fit_lrp(X, CV1)
#' fit
#'
#' ## Data-frame interface with one model per trial
#' df <- rbind(
#'   data.frame(x = X, cv = CV1, trial = "T1"),
#'   data.frame(x = X, cv = CV1 * 1.1, trial = "T2")
#' )
#' res <- fit_lrp(df, x = "x", cv = "cv", trial = "trial")
#' res$summary
#'
#' @seealso [plot.lrp_fit()], [predict.lrp_fit()]
#' @export
fit_lrp <- function(.data = NULL, x = NULL, cv = NULL, trial = NULL,
                    step = 0.001, method = c("segment", "ramp"),
                    search_range = NULL) {

  method <- match.arg(method)

  ## Positional vector call: fit_lrp(x_vec, cv_vec)
  if (!is.null(.data) && !is.data.frame(.data)) {
    cv <- x; x <- .data; .data <- NULL
  }

  ## ---- vector interface ----
  if (is.null(.data)) {
    if (is.null(x) || is.null(cv))
      stop("Provide numeric `x` and `cv`, or a data frame as `.data`.",
           call. = FALSE)
    return(.lrp_fit_one(x, cv, step, method, search_range))
  }

  ## ---- data-frame interface ----
  if (!is.data.frame(.data))
    stop("`.data` must be a data frame.", call. = FALSE)

  xcol  <- if (is.null(x))  "x"  else x
  cvcol <- if (is.null(cv)) "cv" else cv
  need  <- c(xcol, cvcol, if (!is.null(trial)) trial)
  missing_cols <- setdiff(need, names(.data))
  if (length(missing_cols))
    stop(sprintf("Column(s) not found in `.data`: %s.\n  Available columns: %s.",
                 paste(missing_cols, collapse = ", "),
                 paste(names(.data), collapse = ", ")), call. = FALSE)

  if (is.null(trial)) {
    message(sprintf("Using x = '%s', cv = '%s' (single series).", xcol, cvcol))
    return(.lrp_fit_one(.data[[xcol]], .data[[cvcol]], step, method,
                        search_range))
  }

  groups <- split(.data, .data[[trial]])
  message(sprintf("Using x = '%s', cv = '%s', trial = '%s' -> %d trials.",
                  xcol, cvcol, trial, length(groups)))

  fits <- lapply(groups, function(g) .lrp_fit_one(g[[xcol]], g[[cvcol]],
                                                  step, method, search_range))
  summ <- do.call(rbind, Map(function(f, nm) data.frame(
    trial      = nm,
    a          = round(unname(f$coefficients["a"]), 4),
    b          = round(unname(f$coefficients["b"]), 4),
    breakpoint = round(unname(f$parameters["Breakpoint"]), 4),
    plateau    = round(unname(f$parameters["Breakpoint_Response"]), 4),
    R2         = round(unname(f$parameters["R2"]), 4),
    RMSE       = round(unname(f$parameters["RMSE"]), 4),
    AIC        = round(unname(f$parameters["AIC"]), 3),
    BIC        = round(unname(f$parameters["BIC"]), 3),
    stringsAsFactors = FALSE), fits, names(fits)))
  rownames(summ) <- NULL

  structure(list(fits = fits, summary = summ, method = method),
            class = "lrp_multi")
}

## ----------------------------------------------------------------------------
## Methods
## ----------------------------------------------------------------------------
#' Predictions from an LRP fit
#'
#' @param object an object of class \code{"lrp_fit"}.
#' @param newx numeric vector of predictor values. Defaults to the fitted data.
#' @param ... ignored.
#' @return A numeric vector of predicted responses.
#' @export
predict.lrp_fit <- function(object, newx = NULL, ...) {
  a  <- object$coefficients["a"]; b <- object$coefficients["b"]
  xo <- object$parameters["Breakpoint"]
  if (is.null(newx)) newx <- object$data$x
  unname(ifelse(newx < xo, a + b * newx, a + b * xo))
}

#' Print an LRP fit
#'
#' @param x an object of class \code{"lrp_fit"}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.lrp_fit <- function(x, ...) {
  cat("Linear Response Plateau (LRP) fit\n")
  cat("Method:                 ", x$method, "\n")
  cat("Breakpoint (Xo):        ", sprintf("%.3f", x$parameters["Breakpoint"]), "\n")
  cat("CV at breakpoint:       ", sprintf("%.3f", x$parameters["Breakpoint_Response"]), "\n")
  cat("R2:", sprintf("%.3f", x$parameters["R2"]),
      " RMSE:", sprintf("%.3f", x$parameters["RMSE"]),
      " AIC:", sprintf("%.1f", x$parameters["AIC"]),
      " BIC:", sprintf("%.1f", x$parameters["BIC"]), "\n")
  invisible(x)
}

#' Summarize an LRP fit
#'
#' @param object an object of class \code{"lrp_fit"}.
#' @param ... ignored.
#' @return \code{object}, invisibly.
#' @export
summary.lrp_fit <- function(object, ...) {
  cat("Model coefficients:\n"); print(object$coefficients)
  cat("\nGoodness of fit:\n");   print(object$parameters)
  invisible(object)
}

#' Plot an LRP fit (article style)
#'
#' Draws the observed points, the fitted broken line, the breakpoint marker and,
#' optionally, the plateau-model annotations (equations, \eqn{R^2}, Xo, CVxo) in
#' the layout commonly used in plot-size articles. Axis limits are taken from the
#' data.
#'
#' @param x an object of class \code{"lrp_fit"}.
#' @param title plot title.
#' @param annotate_model logical; if \code{TRUE} (default) the model equations
#'   and statistics are drawn on the panel.
#' @param base_size base font size passed to the theme. Increase it for large
#'   exports (for example, a high-resolution TIFF).
#' @param label_size size of the annotation text.
#' @param family font family for the theme and annotations.
#' @param ... ignored.
#' @return A \code{ggplot} object.
#' @export
plot.lrp_fit <- function(x, title = "Linear Plateau", annotate_model = TRUE,
                         base_size = 14, label_size = 4, family = "serif", ...) {
  fmt <- function(v) sprintf("%.2f", v)
  d   <- x$data
  a   <- unname(x$coefficients["a"]); b <- unname(x$coefficients["b"])
  xo  <- unname(x$parameters["Breakpoint"])
  p   <- unname(x$parameters["Breakpoint_Response"])
  r2  <- unname(x$parameters["R2"])

  xmax <- max(d$x) * 1.05
  ymax <- max(d$cv) * 1.10
  curve_x <- seq(0, max(d$x), length.out = 400)
  curve   <- data.frame(x = curve_x, cv = predict(x, curve_x))

  g <- ggplot2::ggplot(d, ggplot2::aes(x, cv)) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::geom_line(data = curve, ggplot2::aes(x, cv), linewidth = 0.9) +
    ggplot2::annotate("segment", x = 0, xend = xo, y = p, yend = p, linetype = 5) +
    ggplot2::annotate("segment", x = xo, xend = xo, y = 0, yend = p, linetype = 5) +
    ggplot2::annotate("point", x = xo, y = p, colour = "red", size = 3.2)

  if (annotate_model) {
    sgn <- ifelse(b < 0, "-", "+")
    labs_txt <- c(
      paste0("CV(x)==", fmt(a), sgn, fmt(abs(b)), "*X~~'if'~~X<=", fmt(xo)),
      paste0("CV(x)==", fmt(p), "~~'if'~~X>", fmt(xo)),
      paste0("R^2==", fmt(r2)),
      paste0("Xo==", fmt(xo)),
      paste0("CVxo==", fmt(p))
    )
    xt <- 0.55 * xmax
    yt <- ymax * c(0.97, 0.88, 0.79, 0.70, 0.61)
    for (i in seq_along(labs_txt))
      g <- g + ggplot2::annotate("text", x = xt, y = yt[i], label = labs_txt[i],
                                 parse = TRUE, hjust = 0.5, size = label_size,
                                 family = family)
  }

  g +
    ggplot2::scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::labs(title = title,
                  x = expression("Plot size (" * m^2 * ")"), y = "CV (%)") +
    ggplot2::theme_classic(base_size = base_size, base_family = family) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

#' Print LRP fits for several trials
#'
#' @param x an object of class \code{"lrp_multi"}.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @export
print.lrp_multi <- function(x, ...) {
  cat(sprintf("LRP fits for %d trials (method = %s)\n\n",
              length(x$fits), x$method))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' Summarize LRP fits for several trials
#'
#' @param object an object of class \code{"lrp_multi"}.
#' @param ... ignored.
#' @return \code{object}, invisibly.
#' @export
summary.lrp_multi <- function(object, ...) {
  print(object$summary, row.names = FALSE)
  invisible(object)
}

#' Plot LRP fits for several trials
#'
#' Builds one article-style plot per trial. Combine them with a package such as
#' \pkg{patchwork} or \pkg{ggpubr} for a multi-panel figure.
#'
#' @param x an object of class \code{"lrp_multi"}.
#' @param ... further arguments passed to [plot.lrp_fit()].
#' @return A named list of \code{ggplot} objects, invisibly.
#' @export
plot.lrp_multi <- function(x, ...) {
  plots <- Map(function(f, nm) plot(f, title = as.character(nm), ...),
               x$fits, names(x$fits))
  invisible(plots)
}
