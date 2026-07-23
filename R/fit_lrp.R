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
                         search_range = NULL, start = NULL,
                         local_min_tol = 0.10) {

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
  ss_profile <- vapply(grid, function(x0) ss_at(x0)$ss, numeric(1))

  if (!any(is.finite(ss_profile)))
    stop("Grid search failed to find a valid breakpoint.", call. = FALSE)

  ## global optimum
  i_best <- which.min(ss_profile)
  best <- ss_at(grid[i_best])
  a <- best$a; b <- best$b; x0 <- grid[i_best]; p <- best$p
  ss_best <- ss_profile[i_best]

  ## ---- local minima of the SSE profile -------------------------------------
  ## Competing basins matter: with few distinct plot sizes the profile is
  ## stepped, and gradient-based fitters (nls, nlsLM) settle in whichever basin
  ## their starting value falls into. Reporting them shows the user that more
  ## than one breakpoint is defensible.
  local_minima <- NULL
  fin <- is.finite(ss_profile)
  if (sum(fin) > 3) {
    gg <- grid[fin]; pp <- ss_profile[fin]
    dsign <- sign(diff(pp))
    dsign[dsign == 0] <- 1                    # flat steps count as ascending
    idx <- which(diff(dsign) == 2) + 1        # falling then rising
    if (length(idx)) {
      lm_x0 <- gg[idx]; lm_ss <- pp[idx]
      keep <- abs(lm_x0 - x0) > 10 * step     # drop the global one itself
      lm_x0 <- lm_x0[keep]; lm_ss <- lm_ss[keep]
      if (length(lm_x0)) {
        ord <- order(lm_ss)
        local_minima <- data.frame(
          breakpoint = lm_x0[ord],
          SSE        = lm_ss[ord],
          SSE_excess = lm_ss[ord] / ss_best - 1,
          row.names  = NULL)
      }
    }
  }

  ## ---- optional compatibility fit ------------------------------------------
  ## `start` reproduces what a gradient fitter seeded at that breakpoint would
  ## return: the local minimum of the basin containing `start`.
  compat <- NULL
  if (!is.null(start)) {
    if (!is.numeric(start) || length(start) != 1 || is.na(start))
      stop("`start` must be a single numeric breakpoint value.", call. = FALSE)
    if (start < min(x) || start > max(x))
      stop(sprintf("`start` must fall within the data range [%.4g, %.4g].",
                   min(x), max(x)), call. = FALSE)

    ## walk downhill from the grid point nearest `start`
    i <- which.min(abs(grid - start))
    repeat {
      lo <- if (i > 1) ss_profile[i - 1] else Inf
      hi <- if (i < length(grid)) ss_profile[i + 1] else Inf
      cur <- ss_profile[i]
      if (!is.finite(cur)) { i <- i + 1; next }
      if (is.finite(lo) && lo < cur) { i <- i - 1
      } else if (is.finite(hi) && hi < cur) { i <- i + 1
      } else break
    }
    cs <- ss_at(grid[i])
    compat <- list(start = start, breakpoint = grid[i],
                   coefficients = c(a = cs$a, b = cs$b),
                   plateau = cs$p, SSE = cs$ss,
                   SSE_excess = cs$ss / ss_best - 1)
  }

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

  ## competing local minimum close in fit: the breakpoint is not well identified
  if (!is.null(local_minima) &&
      any(local_minima$SSE_excess <= local_min_tol)) {
    j <- which(local_minima$SSE_excess <= local_min_tol)[1]
    warning(sprintf(paste0("Competing local minimum at breakpoint %.3f ",
                           "(SSE %.1f%% above the optimum at %.3f). The ",
                           "breakpoint is not sharply identified; see ",
                           "$local_minima and consider `search_range`."),
                    local_minima$breakpoint[j],
                    100 * local_minima$SSE_excess[j], x0), call. = FALSE)
  }

  structure(
    list(
      coefficients = c(a = a, b = b),
      parameters   = c(Breakpoint = x0, Breakpoint_Response = p,
                       R2 = r_squared, RMSE = rmse, AIC = aic, BIC = bic),
      fitted = fitted, residuals = residuals,
      data = data.frame(x = x, cv = cv), method = method, step = step,
      search_range = search_range,
      local_minima = local_minima, compat = compat,
      sse_profile = data.frame(breakpoint = grid, SSE = ss_profile)
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
#' @section Competing breakpoints:
#' With few distinct plot sizes the residual sum of squares is a stepped
#' function of the breakpoint, and it often has several local minima. The grid
#' search always returns the global optimum, but a second basin may fit almost
#' as well, in which case the breakpoint is not sharply identified and different
#' implementations legitimately disagree. Any competing minima are reported in
#' \code{$local_minima}, and a warning is issued when one of them fits within
#' \code{local_min_tol} of the optimum. Gradient-based fitters return whichever
#' basin their starting value lands in; \code{start} reproduces that.
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
#' @param start optional single breakpoint value. The fit is unchanged, but the
#'   result also carries a \code{$compat} element holding the local minimum of
#'   the basin containing \code{start}: the solution a gradient-based fitter
#'   (\code{nls}, \code{nlsLM}) seeded there would return. Use it to reproduce
#'   published results obtained with such implementations, and to see how much
#'   worse they fit.
#' @param local_min_tol relative SSE tolerance (default 0.10) for warning about
#'   a competing local minimum. A second basin fitting within 10\% of the
#'   optimum means the breakpoint is not sharply identified.
#'
#' @return
#' For a single series, an object of class \code{"lrp_fit"}: a list with
#' \code{coefficients} (a, b); \code{parameters} (Breakpoint, Breakpoint_Response,
#' R2, RMSE, AIC, BIC); \code{fitted}; \code{residuals}; \code{data};
#' \code{method}; \code{step}; \code{local_minima} (competing basins, with
#' their SSE excess over the optimum, or \code{NULL}); \code{sse_profile} (the
#' SSE at every candidate breakpoint); and \code{compat} when \code{start} was
#' given. \cr
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
                    search_range = NULL, start = NULL,
                    local_min_tol = 0.10) {

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
    return(.lrp_fit_one(x, cv, step, method, search_range, start,
                        local_min_tol))
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
                        search_range, start, local_min_tol))
  }

  groups <- split(.data, .data[[trial]])
  message(sprintf("Using x = '%s', cv = '%s', trial = '%s' -> %d trials.",
                  xcol, cvcol, trial, length(groups)))

  fits <- lapply(groups, function(g) .lrp_fit_one(g[[xcol]], g[[cvcol]],
                                                  step, method, search_range,
                                                  start, local_min_tol))
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
    n_local    = if (is.null(f$local_minima)) 0L else nrow(f$local_minima),
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

  if (!is.null(x$local_minima)) {
    lmin <- x$local_minima
    cat("\nCompeting local minima (", nrow(lmin), "):\n", sep = "")
    show <- utils::head(lmin, 3)
    for (i in seq_len(nrow(show)))
      cat(sprintf("  Xo = %7.3f   SSE %+6.1f%% vs optimum\n",
                  show$breakpoint[i], 100 * show$SSE_excess[i]))
    if (nrow(lmin) > 3) cat("  ... see $local_minima for all\n")
  }

  if (!is.null(x$compat))
    cat(sprintf("\nCompatibility fit (start = %.3f): Xo = %.3f, SSE %+.1f%%\n",
                x$compat$start, x$compat$breakpoint,
                100 * x$compat$SSE_excess))
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

## Save helper (internal) ----------------------------------------------------
.save_lrp <- function(g, file, base, format, dpi, width, height, units,
                      compression) {
  format <- match.arg(format, c("tiff", "png", "jpeg", "pdf", "eps"))
  if (is.null(file)) {
    slug <- gsub("[^A-Za-z0-9]+", "_", base)
    slug <- gsub("^_|_$", "", slug)
    if (!nzchar(slug)) slug <- "lrp_plot"
    file <- paste0(slug, ".", format)
  }
  args <- list(filename = file, plot = g, dpi = dpi,
               width = width, height = height, units = units)
  if (format %in% c("tiff", "tif")) args$compression <- compression
  do.call(ggplot2::ggsave, args)
  message("Saved: ", file)
  invisible(file)
}

## Article-style plot --------------------------------------------------------
#' Plot an LRP fit (publication style)
#'
#' Draws the observed points, the fitted broken line, dotted guides to the
#' breakpoint and the plateau-model annotations, in the layout used in
#' plot-size articles: the model block (equations and \eqn{R^2}) centred at the
#' top and \code{Xo}/\code{CVxo} next to the breakpoint. Axis limits come from
#' the data.
#'
#' @param x an object of class \code{"lrp_fit"}.
#' @param title plot title.
#' @param annotate_model logical; draw the model equations and statistics.
#' @param xlab,ylab axis titles. \code{xlab = NULL} uses "Plot size (m^2)".
#' @param decimal_mark decimal separator for the annotations, "." or ",".
#' @param cond_word conditional word in the equations (e.g. "if" or "se").
#' @param digits_coef,digits_stat decimals for the coefficients (a, b) and for
#'   the statistics (Xo, CVxo, R2).
#' @param point_size,line_size sizes of the points and the fitted line.
#' @param point_colour,line_colour,bp_colour colours of the points, the fitted
#'   line and the breakpoint marker.
#' @param base_size base font size for the theme; axis titles and text scale
#'   from it.
#' @param label_size size of the annotation text (equations, Xo, CVxo).
#' @param title_size size of the plot title. \code{NULL} uses \code{base_size}.
#' @param family font family for the theme and annotations.
#' @param theme optional \pkg{ggplot2} theme object used instead of the default
#'   \code{theme_classic(base_size)}; the axis-line, tick and centred-title
#'   tweaks are applied on top.
#' @param save logical; if \code{TRUE}, write the figure to disk (default FALSE).
#' @param file output file name; if \code{NULL}, derived from \code{title}.
#' @param format one of "tiff", "png", "jpeg", "pdf", "eps". TIFF/PDF/EPS are
#'   preferred for journals; TIFF is written with LZW compression.
#' @param dpi resolution for raster formats.
#' @param width,height,units figure size.
#' @param compression TIFF compression (default "lzw").
#' @param ... ignored.
#' @return A \code{ggplot} object (invisibly when saved).
#' @export
plot.lrp_fit <- function(x, title = "Linear Plateau", annotate_model = TRUE,
                         xlab = NULL, ylab = "CV (%)",
                         decimal_mark = c(".", ","), cond_word = "if",
                         digits_coef = 3, digits_stat = 2,
                         point_size = 2.3, line_size = 0.8,
                         point_colour = "black", line_colour = "black",
                         bp_colour = "red",
                         base_size = 13, label_size = 4.6, title_size = NULL,
                         family = "serif", theme = NULL,
                         save = FALSE, file = NULL,
                         format = c("tiff", "png", "jpeg", "pdf", "eps"),
                         dpi = 300, width = 18, height = 12, units = "cm",
                         compression = "lzw", ...) {

  decimal_mark <- match.arg(decimal_mark)
  format       <- match.arg(format)
  if (is.null(title_size)) title_size <- base_size

  fmtn <- function(v, d) {
    s <- formatC(v, format = "f", digits = d)
    if (decimal_mark == ",") s <- gsub("\\.", ",", s)
    s
  }

  d  <- x$data
  a  <- unname(x$coefficients["a"]); b <- unname(x$coefficients["b"])
  xo <- unname(x$parameters["Breakpoint"])
  p  <- unname(x$parameters["Breakpoint_Response"])
  r2 <- unname(x$parameters["R2"])

  xmax <- max(d$x) * 1.05
  ymax <- max(d$cv) * 1.12
  curve <- data.frame(x = seq(0, max(d$x), length.out = 400))
  curve$cv <- predict(x, curve$x)

  if (is.null(xlab)) xlab <- expression("Plot size (" * m^2 * ")")

  g <- ggplot2::ggplot(d, ggplot2::aes(x, cv)) +
    ggplot2::annotate("segment", x = 0, xend = xo, y = p, yend = p,
                      linetype = 3, linewidth = 0.5) +
    ggplot2::annotate("segment", x = xo, xend = xo, y = 0, yend = p,
                      linetype = 3, linewidth = 0.5) +
    ggplot2::geom_line(data = curve, ggplot2::aes(x, cv),
                       linewidth = line_size, colour = line_colour) +
    ggplot2::geom_point(size = point_size, colour = point_colour) +
    ggplot2::annotate("point", x = xo, y = p, colour = bp_colour,
                      size = point_size * 1.3)

  if (annotate_model) {
    sgn  <- if (b < 0) "-" else "+"
    a_s  <- fmtn(a, digits_coef);  b_s <- fmtn(abs(b), digits_coef)
    xo_s <- fmtn(xo, digits_stat); p_s <- fmtn(p, digits_stat)
    r2_s <- fmtn(r2, digits_stat)

    eq1 <- paste0("CV[(x)]=='", a_s, "'", sgn, "'", b_s, "'*X~~'",
                  cond_word, "'~~X<='", xo_s, "'")
    eq2 <- paste0("CV[(x)]=='", p_s, "'~~'", cond_word, "'~~X>'", xo_s, "'")
    eq3 <- paste0("R^2=='", r2_s, "'")
    lXo <- paste0("X[o]=='", xo_s, "'")
    lCV <- paste0("CV[Xo]=='", p_s, "'")

    ## model block: centred, top
    xm <- 0.55 * xmax
    ym <- ymax * c(0.97, 0.89, 0.81)
    for (i in seq_len(3))
      g <- g + ggplot2::annotate("text", x = xm, y = ym[i],
                                 label = c(eq1, eq2, eq3)[i], parse = TRUE,
                                 hjust = 0.5, size = label_size, family = family)

    ## breakpoint block: right of the vertical guide, below the plateau
    xb <- xo + 0.03 * max(d$x)
    yb <- c(max(p - 0.11 * ymax, 0.10 * ymax),
            max(p - 0.20 * ymax, 0.02 * ymax))
    for (i in 1:2)
      g <- g + ggplot2::annotate("text", x = xb, y = yb[i],
                                 label = c(lXo, lCV)[i], parse = TRUE,
                                 hjust = 0, size = label_size, family = family)
  }

  g <- g +
    ggplot2::scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    (if (is.null(theme))
      ggplot2::theme_classic(base_size = base_size, base_family = family)
     else theme) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = title_size),
      axis.line  = ggplot2::element_line(colour = "black", linewidth = 0.6),
      axis.ticks = ggplot2::element_line(colour = "black"),
      axis.text  = ggplot2::element_text(colour = "black")
    )

  if (save) {
    .save_lrp(g, file, title, format, dpi, width, height, units, compression)
    return(invisible(g))
  }
  g
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

#' Plot LRP fits for several trials (paginated grid)
#'
#' Arranges the per-trial plots in a grid of at most 6 panels (3 rows by 2
#' columns). With more than 6 trials the panels are paginated: each page is a
#' separate 3x2 figure, and with \code{save = TRUE} each page is written to its
#' own file. Requires the \pkg{patchwork} package.
#'
#' @param x an object of class \code{"lrp_multi"}.
#' @param save logical; write the page(s) to disk (default FALSE).
#' @param file base file name; page number and extension are appended when there
#'   is more than one page.
#' @param format,dpi,width,height,units,compression passed to the saver.
#' @param ... styling arguments forwarded to [plot.lrp_fit()]
#'   (for example \code{decimal_mark}, \code{cond_word}, \code{label_size}).
#' @return A list of page objects, invisibly.
#' @export
plot.lrp_multi <- function(x, save = FALSE, file = NULL,
                           format = c("tiff", "png", "jpeg", "pdf", "eps"),
                           dpi = 300, width = 18, height = 24, units = "cm",
                           compression = "lzw", ...) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required to arrange multiple trials. ",
         "Install it with install.packages('patchwork').", call. = FALSE)
  format <- match.arg(format)

  plots <- Map(function(f, nm) plot(f, title = as.character(nm), save = FALSE, ...),
               x$fits, names(x$fits))

  idx    <- seq_along(plots)
  chunks <- split(idx, ceiling(idx / 6))
  pages  <- lapply(chunks, function(ix)
    patchwork::wrap_plots(plots[ix], nrow = 3, ncol = 2))

  if (save) {
    base <- if (is.null(file)) "lrp_trials" else sub("\\.[^.]*$", "", file)
    n <- length(pages)
    for (i in seq_len(n)) {
      fn <- if (n == 1) paste0(base, ".", format)
      else sprintf("%s_p%d.%s", base, i, format)
      .save_lrp(pages[[i]], fn, base, format, dpi, width, height, units,
                compression)
    }
  } else {
    for (pg in pages) print(pg)
  }
  invisible(pages)
}
