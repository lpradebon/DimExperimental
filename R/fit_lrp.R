#' Fit a Linear Response Plateau (LRP) model by grid search
#'
#' Estimates the linear-plateau (broken-line) model
#' \deqn{f(x) = a + b x \text{ for } x \le X_0, \qquad f(x) = a + b X_0 \text{ for } x > X_0}
#' by profiling the breakpoint \eqn{X_0} over a grid. For each candidate
#' breakpoint the linear coefficients are obtained by least squares and the
#' plateau level is set to \eqn{a + b X_0}. The breakpoint that minimizes the
#' residual sum of squares over all observations is returned. This avoids the
#' starting-value and local-minimum problems of a direct \code{nls} fit, so no
#' initial values are needed.
#'
#' The response is typically the coefficient of variation (CV, percent) and the
#' predictor the plot size. Individual CV values (one per basic-unit form), not
#' means per plot size, must be supplied. Repeated \code{x} values are expected
#' and intended.
#'
#' @param x numeric predictor (for example, plot size). Repeats allowed.
#' @param cv numeric response (for example, CV in percent), same length as \code{x}.
#' @param step grid step for the breakpoint search. Default 0.001, matching the
#'   validation routine. Larger values are faster with slightly coarser Xo.
#' @param method how the linear coefficients are estimated at each candidate
#'   breakpoint. \code{"segment"} (default) uses only observations with
#'   \code{x <= X0}, reproducing the Paranaiba et al. (2009) validation.
#'   \code{"ramp"} uses all observations on the basis \code{pmin(x, X0)}, the
#'   standard least-squares LRP.
#'
#' @return An object of class \code{"lrp_fit"}: a list with \code{coefficients}
#'   (a, b); \code{parameters} (Breakpoint, Breakpoint_Response, R2, RMSE, AIC,
#'   BIC); \code{fitted}; \code{residuals}; \code{data}; \code{method}; \code{step}.
#' @export
fit_lrp <- function(x, cv, step = 0.001, method = c("segment", "ramp")) {

  method <- match.arg(method)

  ## ---- input validation ----
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

  ## ---- residual sum of squares for a fixed breakpoint ----
  ss_at <- function(x0) {
    if (method == "segment") {
      m <- x <= x0
      if (sum(m) < 2) return(list(ss = Inf))
      cf <- stats::coef(stats::lm(cv[m] ~ x[m]))
    } else {
      z  <- pmin(x, x0)
      cf <- stats::coef(stats::lm(cv ~ z))
    }
    a <- unname(cf[1]); b <- unname(cf[2])
    p <- a + b * x0
    fit <- ifelse(x <= x0, a + b * x, p)
    list(ss = sum((cv - fit)^2), a = a, b = b, p = p)
  }

  ## ---- grid search over the breakpoint ----
  grid <- seq(x_unique[2], x_unique[length(x_unique)], by = step)
  best <- list(ss = Inf)
  for (x0 in grid) {
    cand <- ss_at(x0)
    if (cand$ss < best$ss) best <- c(cand, list(x0 = x0))
  }
  if (!is.finite(best$ss))
    stop("Grid search failed to find a valid breakpoint.", call. = FALSE)

  a <- best$a; b <- best$b; x0 <- best$x0; p <- best$p

  ## ---- fitted values and fit statistics ----
  fitted    <- ifelse(x <= x0, a + b * x, p)
  residuals <- cv - fitted
  n         <- length(cv)
  rss       <- sum(residuals^2)
  r_squared <- 1 - rss / sum((cv - mean(cv))^2)
  rmse      <- sqrt(mean(residuals^2))

  ## AIC / BIC via the Gaussian log-likelihood (MLE variance = rss / n).
  ## Parameter count k = a, b, Xo, sigma = 4, to match the nls / AIC.default
  ## convention. Only comparable with the QRP fit if it uses the same count.
  loglik <- -0.5 * n * (log(2 * pi) + log(rss / n) + 1)
  k      <- 4L
  aic    <- -2 * loglik + 2 * k
  bic    <- -2 * loglik + log(n) * k

  ## ---- plausibility warnings ----
  if (b >= 0)
    warning("Estimated slope is non-negative; the decreasing-plateau ",
            "interpretation may not hold for these data.", call. = FALSE)
  if (isTRUE(all.equal(x0, x_unique[2])) ||
      isTRUE(all.equal(x0, x_unique[length(x_unique)])))
    warning("Breakpoint is at the edge of the search range; the data may not ",
            "capture the full descent or the plateau.", call. = FALSE)

  structure(
    list(
      coefficients = c(a = a, b = b),
      parameters   = c(Breakpoint = x0, Breakpoint_Response = p,
                       R2 = r_squared, RMSE = rmse, AIC = aic, BIC = bic),
      fitted       = fitted,
      residuals    = residuals,
      data         = data.frame(x = x, cv = cv),
      method       = method,
      step         = step
    ),
    class = "lrp_fit"
  )
}

#' Predict from an LRP fit
#'
#' @param object an object of class \code{"lrp_fit"}.
#' @param newx numeric vector of predictor values. Defaults to the fitted data.
#' @param ... ignored.
#' @return numeric vector of predicted responses.
#' @export
predict.lrp_fit <- function(object, newx = NULL, ...) {
  a  <- object$coefficients["a"]
  b  <- object$coefficients["b"]
  xo <- object$parameters["Breakpoint"]
  if (is.null(newx)) newx <- object$data$x
  unname(ifelse(newx < xo, a + b * newx, a + b * xo))
}

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

#' @export
summary.lrp_fit <- function(object, ...) {
  cat("Model coefficients:\n")
  print(object$coefficients)
  cat("\nGoodness of fit:\n")
  print(object$parameters)
  invisible(object)
}

#' @export
plot.lrp_fit <- function(x, ...) {
  d  <- x$data
  xo <- unname(x$parameters["Breakpoint"])
  p  <- unname(x$parameters["Breakpoint_Response"])
  r2 <- unname(x$parameters["R2"])

  curve_x <- seq(min(d$x), max(d$x), length.out = 400)
  curve   <- data.frame(x = curve_x, cv = predict(x, curve_x))

  g <- ggplot2::ggplot(d, ggplot2::aes(x, cv)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(data = curve, ggplot2::aes(x, cv), linewidth = 0.9) +
    ggplot2::geom_segment(x = 0, xend = xo, y = p, yend = p,
                          linetype = 2) +
    ggplot2::geom_segment(x = xo, xend = xo, y = -Inf, yend = p,
                          linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(x = xo, y = p), colour = "red", size = 3) +
    ggplot2::labs(
      x = expression("Plot size (" * m^2 * ")"), y = "CV (%)",
      subtitle = sprintf("Xo = %.2f    CV(Xo) = %.2f    R2 = %.2f", xo, p, r2)
    ) +
    ggplot2::theme_classic()

  print(g)
  invisible(g)
}
