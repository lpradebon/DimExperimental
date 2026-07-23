## ============================================================================
## DimExp :: Quadratic Response Plateau (QRP) model by grid search
## ----------------------------------------------------------------------------
## Mirrors fit_lrp(): same interfaces, S3 classes, publication plot and saving.
## The quadratic-plateau joins smoothly (vertex at the breakpoint), so fixing
## the breakpoint X0 makes the model linear in (A, c):
##   x <= X0:  cv = A + c (x - X0)^2      x > X0:  cv = A   (plateau)
## The breakpoint is profiled over a grid, avoiding nls starting values.
## Relation to the a + b x + c x^2 form: a = A + c X0^2, b = -2 c X0.
## Uses the shared internal saver .save_lrp() defined in fit_lrp.R.
## ============================================================================

## Internal single-series engine ---------------------------------------------
.qrp_fit_one <- function(x, cv, step = 0.001, search_range = NULL,
                         start = NULL, local_min_tol = 0.10) {

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
                          "breakpoint minimum (%.4g)."), feas_lo), call. = FALSE)
  }

  ## residual sum of squares for a fixed breakpoint (linear in A and c)
  ss_at <- function(x0) {
    z  <- ifelse(x <= x0, (x - x0)^2, 0)
    cf <- stats::coef(stats::lm(cv ~ z))
    A  <- unname(cf[1]); cc <- unname(cf[2])
    fit <- A + ifelse(x <= x0, cc * (x - x0)^2, 0)
    list(ss = sum((cv - fit)^2), A = A, cc = cc)
  }

  grid <- seq(grid_lo, grid_hi, by = step)
  ss_profile <- vapply(grid, function(x0) ss_at(x0)$ss, numeric(1))

  if (!any(is.finite(ss_profile)))
    stop("Grid search failed to find a valid breakpoint.", call. = FALSE)

  i_best <- which.min(ss_profile)
  best <- ss_at(grid[i_best])
  x0 <- grid[i_best]; A <- best$A; cc <- best$cc
  ss_best <- ss_profile[i_best]

  ## ---- local minima of the SSE profile -------------------------------------
  ## The quadratic-plateau joins smoothly, so this profile is usually a single
  ## smooth basin (unlike the LRP, whose kink makes it stepped). Competing
  ## minima are still reported when they occur.
  local_minima <- NULL
  fin <- is.finite(ss_profile)
  if (sum(fin) > 3) {
    gg <- grid[fin]; pp <- ss_profile[fin]
    dsign <- sign(diff(pp)); dsign[dsign == 0] <- 1
    idx <- which(diff(dsign) == 2) + 1
    if (length(idx)) {
      lm_x0 <- gg[idx]; lm_ss <- pp[idx]
      keep <- abs(lm_x0 - x0) > 10 * step
      lm_x0 <- lm_x0[keep]; lm_ss <- lm_ss[keep]
      if (length(lm_x0)) {
        ord <- order(lm_ss)
        local_minima <- data.frame(
          breakpoint = lm_x0[ord], SSE = lm_ss[ord],
          SSE_excess = lm_ss[ord] / ss_best - 1, row.names = NULL)
      }
    }
  }

  ## ---- optional compatibility fit ------------------------------------------
  compat <- NULL
  if (!is.null(start)) {
    if (!is.numeric(start) || length(start) != 1 || is.na(start))
      stop("`start` must be a single numeric breakpoint value.", call. = FALSE)
    if (start < min(x) || start > max(x))
      stop(sprintf("`start` must fall within the data range [%.4g, %.4g].",
                   min(x), max(x)), call. = FALSE)
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
                   plateau = cs$A, SSE = cs$ss,
                   SSE_excess = cs$ss / ss_best - 1)
  }
  a <- A + cc * x0^2; b <- -2 * cc * x0; c_par <- cc
  breakpoint <- x0; plateau <- A

  fitted    <- A + ifelse(x <= x0, cc * (x - x0)^2, 0)
  residuals <- cv - fitted
  n   <- length(cv)
  p   <- 3L
  sse <- sum(residuals^2)
  mse <- mean(residuals^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(residuals))
  r2      <- 1 - sse / sum((cv - mean(cv))^2)
  r2_adj  <- 1 - (1 - r2) * (n - 1) / (n - p)
  loglik  <- -0.5 * n * (log(2 * pi) + log(sse / n) + 1)
  k   <- 4L
  aic <- -2 * loglik + 2 * k
  bic <- -2 * loglik + log(n) * k

  if (cc <= 0)
    warning("Estimated quadratic coefficient is non-positive; the ",
            "decreasing-then-plateau shape may not hold for these data.",
            call. = FALSE)
  if (isTRUE(all.equal(x0, grid_lo)) || isTRUE(all.equal(x0, grid_hi)))
    warning("Breakpoint is at the edge of the search range; the data (or the ",
            "supplied `search_range`) may not bracket the true breakpoint.",
            call. = FALSE)

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
      coefficients = c(a = a, b = b, c = c_par),
      parameters   = c(Breakpoint = breakpoint, Breakpoint_Response = plateau,
                       R2 = r2, R2_adj = r2_adj, RMSE = rmse, MAE = mae,
                       AIC = aic, BIC = bic, SSE = sse, MSE = mse),
      fitted = fitted, residuals = residuals,
      data = data.frame(x = x, cv = cv), step = step,
      search_range = search_range,
      local_minima = local_minima, compat = compat,
      sse_profile = data.frame(breakpoint = grid, SSE = ss_profile)
    ),
    class = "qrp_fit"
  )
}

## Public fitter --------------------------------------------------------------
#' Fit the Quadratic Response Plateau (QRP) model by grid search
#'
#' Fits the quadratic-plateau (smooth broken-line) model
#' \deqn{f(x) = a + b x + c x^2 \ \text{ if } x \le X_0, \qquad
#'       f(x) = a - b^2/(4c) \ \text{ if } x > X_0,}
#' with \eqn{X_0 = -b/(2c)}. The breakpoint is profiled over a grid: for each
#' candidate \eqn{X_0} the model is linear in the plateau level and the
#' curvature, so it is fit by least squares with no starting values, and the
#' \eqn{X_0} minimizing the residual sum of squares is returned. This mirrors
#' [fit_lrp()] and avoids the convergence problems of a direct nls fit.
#'
#' @section Two ways to call:
#' \describe{
#'   \item{Vectors}{\code{fit_qrp(x, cv)} returns a \code{"qrp_fit"}.}
#'   \item{Data frame}{\code{fit_qrp(.data, x = "x", cv = "cv", trial = "trial")};
#'     with \code{trial}, one model per trial is fit and a \code{"qrp_multi"}
#'     object is returned. Column names default to "x", "cv", "trial".}
#' }
#'
#' @param .data optional data frame; when supplied, \code{x}/\code{cv}/\code{trial}
#'   are column names.
#' @param x,cv numeric vectors, or column names when \code{.data} is a data frame.
#' @param trial optional column name identifying the trial.
#' @param step grid step for the breakpoint search (default 0.001).
#' @param search_range optional \code{c(lower, upper)} restricting the breakpoint
#'   search; must fall within the data range.
#' @param start optional single breakpoint value. The fit is unchanged, but the
#'   result also carries \code{$compat}: the local minimum of the basin
#'   containing \code{start}, i.e. what a gradient-based fitter seeded there
#'   would return. Rarely needed for the QRP, whose SSE profile is usually a
#'   single smooth basin.
#' @param local_min_tol relative SSE tolerance (default 0.10) for warning about
#'   a competing local minimum.
#' @return A \code{"qrp_fit"} (single series) or \code{"qrp_multi"} (per trial).
#'   The fit also carries \code{local_minima} (competing basins with their SSE
#'   excess over the optimum, or \code{NULL}), \code{sse_profile}, and
#'   \code{compat} when \code{start} was given. Because the quadratic-plateau
#'   joins smoothly, this profile is typically a single basin, unlike the
#'   stepped profile of [fit_lrp()].
#' @seealso [fit_lrp()], [plot.qrp_fit()]
#' @examples
#' X   <- c(1, 2, 2, 3, 3, 4, 6, 6, 6, 6, 9, 12, 12, 18, 18)
#' CV1 <- c(30.40, 19.51, 23.72, 12.89, 21.32, 16.69, 6.71, 10.75,
#'          17.58, 14.94, 11.93, 3.18, 8.63, 4.25, 11.41)
#' fit <- fit_qrp(X, CV1)
#' fit
#' @export
fit_qrp <- function(.data = NULL, x = NULL, cv = NULL, trial = NULL,
                    step = 0.001, search_range = NULL, start = NULL,
                    local_min_tol = 0.10) {

  if (!is.null(.data) && !is.data.frame(.data)) {
    cv <- x; x <- .data; .data <- NULL
  }

  if (is.null(.data)) {
    if (is.null(x) || is.null(cv))
      stop("Provide numeric `x` and `cv`, or a data frame as `.data`.",
           call. = FALSE)
    return(.qrp_fit_one(x, cv, step, search_range, start, local_min_tol))
  }

  if (!is.data.frame(.data)) stop("`.data` must be a data frame.", call. = FALSE)

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
    return(.qrp_fit_one(.data[[xcol]], .data[[cvcol]], step, search_range,
                        start, local_min_tol))
  }

  groups <- split(.data, .data[[trial]])
  message(sprintf("Using x = '%s', cv = '%s', trial = '%s' -> %d trials.",
                  xcol, cvcol, trial, length(groups)))

  fits <- lapply(groups, function(g) .qrp_fit_one(g[[xcol]], g[[cvcol]],
                                                  step, search_range, start,
                                                  local_min_tol))
  summ <- do.call(rbind, Map(function(f, nm) data.frame(
    trial      = nm,
    a          = round(unname(f$coefficients["a"]), 4),
    b          = round(unname(f$coefficients["b"]), 4),
    c          = round(unname(f$coefficients["c"]), 4),
    breakpoint = round(unname(f$parameters["Breakpoint"]), 4),
    plateau    = round(unname(f$parameters["Breakpoint_Response"]), 4),
    R2         = round(unname(f$parameters["R2"]), 4),
    RMSE       = round(unname(f$parameters["RMSE"]), 4),
    AIC        = round(unname(f$parameters["AIC"]), 3),
    BIC        = round(unname(f$parameters["BIC"]), 3),
    n_local    = if (is.null(f$local_minima)) 0L else nrow(f$local_minima),
    stringsAsFactors = FALSE), fits, names(fits)))
  rownames(summ) <- NULL

  structure(list(fits = fits, summary = summ), class = "qrp_multi")
}

## Methods --------------------------------------------------------------------
#' Predictions from a QRP fit
#' @param object a \code{"qrp_fit"} object.
#' @param newx numeric predictor values; defaults to the fitted data.
#' @param ... ignored.
#' @return numeric vector of predicted responses.
#' @export
predict.qrp_fit <- function(object, newx = NULL, ...) {
  a  <- object$coefficients["a"]; b <- object$coefficients["b"]
  cc <- object$coefficients["c"]
  bp <- object$parameters["Breakpoint"]; pl <- object$parameters["Breakpoint_Response"]
  if (is.null(newx)) newx <- object$data$x
  unname(ifelse(newx <= bp, a + b * newx + cc * newx^2, pl))
}

#' @export
print.qrp_fit <- function(x, ...) {
  cat("Quadratic Response Plateau (QRP) fit\n")
  cat("Breakpoint (Xo):        ", sprintf("%.3f", x$parameters["Breakpoint"]), "\n")
  cat("CV at breakpoint:       ", sprintf("%.3f", x$parameters["Breakpoint_Response"]), "\n")
  cat("R2:", sprintf("%.3f", x$parameters["R2"]),
      " R2 adj:", sprintf("%.3f", x$parameters["R2_adj"]),
      " RMSE:", sprintf("%.3f", x$parameters["RMSE"]),
      " MAE:", sprintf("%.3f", x$parameters["MAE"]), "\n")

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

#' @export
summary.qrp_fit <- function(object, ...) {
  cat("Model coefficients:\n"); print(object$coefficients)
  cat("\nGoodness of fit:\n");   print(object$parameters)
  invisible(object)
}

#' Plot a QRP fit (publication style)
#'
#' Same layout and options as [plot.lrp_fit()], with the quadratic descending
#' arm and the quadratic equation in the annotation.
#'
#' @param x a \code{"qrp_fit"} object.
#' @param title plot title.
#' @param annotate_model draw the model equations and statistics.
#' @param xlab,ylab axis titles. \code{xlab = NULL} uses "Plot size (m^2)".
#' @param decimal_mark decimal separator, "." or ",".
#' @param cond_word conditional word in the equations (e.g. "if" or "se").
#' @param digits_coef,digits_c,digits_stat decimals for a/b, for c, and for the
#'   statistics.
#' @param point_size,line_size sizes of points and fitted line.
#' @param point_colour,line_colour,bp_colour point, line and breakpoint colours.
#' @param base_size,label_size,title_size,family theme and annotation sizes and
#'   font family.
#' @param theme optional ggplot2 theme used instead of the default.
#' @param save,file,format,dpi,width,height,units,compression saving controls;
#'   see [plot.lrp_fit()].
#' @param ... ignored.
#' @return A \code{ggplot} object (invisibly when saved).
#' @export
plot.qrp_fit <- function(x, title = "Quadratic Plateau", annotate_model = TRUE,
                         xlab = NULL, ylab = "CV (%)",
                         decimal_mark = c(".", ","), cond_word = "if",
                         digits_coef = 3, digits_c = 4, digits_stat = 2,
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
  cc <- unname(x$coefficients["c"])
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
    sgn_b <- if (b  < 0) "-" else "+"
    sgn_c <- if (cc < 0) "-" else "+"
    a_s  <- fmtn(a, digits_coef);      b_s <- fmtn(abs(b), digits_coef)
    c_s  <- fmtn(abs(cc), digits_c)
    xo_s <- fmtn(xo, digits_stat); p_s <- fmtn(p, digits_stat)
    r2_s <- fmtn(r2, digits_stat)

    eq1 <- paste0("CV[(x)]=='", a_s, "'", sgn_b, "'", b_s, "'*X", sgn_c,
                  "'", c_s, "'*X^2~~'", cond_word, "'~~X<='", xo_s, "'")
    eq2 <- paste0("CV[(x)]=='", p_s, "'~~'", cond_word, "'~~X>'", xo_s, "'")
    eq3 <- paste0("R^2=='", r2_s, "'")
    lXo <- paste0("X[o]=='", xo_s, "'")
    lCV <- paste0("CV[Xo]=='", p_s, "'")

    xm <- 0.55 * xmax
    ym <- ymax * c(0.97, 0.89, 0.81)
    for (i in seq_len(3))
      g <- g + ggplot2::annotate("text", x = xm, y = ym[i],
                                 label = c(eq1, eq2, eq3)[i], parse = TRUE,
                                 hjust = 0.5, size = label_size, family = family)

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

#' @export
print.qrp_multi <- function(x, ...) {
  cat(sprintf("QRP fits for %d trials\n\n", length(x$fits)))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @export
summary.qrp_multi <- function(object, ...) {
  print(object$summary, row.names = FALSE)
  invisible(object)
}

#' Plot QRP fits for several trials (paginated grid)
#'
#' Arranges the per-trial plots in a grid of at most 6 panels (3 x 2), paginating
#' when there are more. Requires \pkg{patchwork}. See [plot.lrp_multi()].
#'
#' @param x a \code{"qrp_multi"} object.
#' @param save,file,format,dpi,width,height,units,compression saving controls.
#' @param ... styling arguments forwarded to [plot.qrp_fit()].
#' @return A list of page objects, invisibly.
#' @export
plot.qrp_multi <- function(x, save = FALSE, file = NULL,
                           format = c("tiff", "png", "jpeg", "pdf", "eps"),
                           dpi = 300, width = 18, height = 24, units = "cm",
                           compression = "lzw", ...) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required to arrange multiple trials.",
         call. = FALSE)
  format <- match.arg(format)

  plots <- Map(function(f, nm) plot(f, title = as.character(nm), save = FALSE, ...),
               x$fits, names(x$fits))
  idx    <- seq_along(plots)
  chunks <- split(idx, ceiling(idx / 6))
  pages  <- lapply(chunks, function(ix)
    patchwork::wrap_plots(plots[ix], nrow = 3, ncol = 2))

  if (save) {
    base <- if (is.null(file)) "qrp_trials" else sub("\\.[^.]*$", "", file)
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
