#' Fit a Quadratic Response Plateau (QRP) model for plot size optimization
#'
#' @description
#' Fits a quadratic-plateau model to coefficient of variation (CV) values
#' across a range of plot sizes, returning the breakpoint (optimal plot
#' size), fit statistics, and a standardized ggplot2 visualization. Mirrors
#' the structure of \code{fit_lrp()} for consistency within the package.
#'
#' @param x Numeric vector of plot sizes.
#' @param cv Numeric vector of CV values (percent), same length as \code{x}.
#' @param plot_title Character string used as the plot title. Default \code{NULL}.
#' @param start_b,start_c Starting values for the \code{b} and \code{c}
#'   coefficients in \code{nlsLM}. Defaults \code{-1} and \code{0.05}, matching
#'   the original implementation. NOTE: unlike \code{a} (initialized from
#'   \code{max(cv)}), these are fixed constants, not derived from the data.
#'   Flag for review once tested on datasets with a different CV/plot-size
#'   scale, see Details.
#' @param base_size Base font size for the plot theme. Default \code{14}.
#' @param font_family Font family for all plot text. Default \code{"sans"}.
#' @param axis_expand_x,axis_expand_y Fractional margin added beyond the data
#'   range when building axis limits. Default \code{0.05} (5%) each.
#'
#' @details
#' The fitting method is unchanged from the original script:
#' \code{minpack.lm::nlsLM} with \code{start = list(a = max(cv), b = start_b,
#' c = start_c)}. This was kept deliberately "as is" for this round of
#' standardization, same approach used for \code{fit_lrp()}.
#'
#' @return An object of class \code{qrp_fit}, a list with:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{model}{The fitted \code{nlsLM} object.}
#'   \item{coefficients}{Coefficient table from \code{summary(model)}.}
#'   \item{parameters}{Named vector: Breakpoint, Breakpoint_Response, R2,
#'     R2_adj, RMSE, MAE, AIC, BIC, SSE, MSE.}
#'   \item{data}{Original \code{x}/\code{cv} data as a data frame.}
#'   \item{curve}{Fitted curve, 500 points, for plotting.}
#'   \item{plot}{A ggplot2 object. Not printed automatically.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- fit_qrp(x = trial1$x, cv = trial1$cv, plot_title = "Trial 1")
#' fit
#' plot(fit)
#' }
#'
#' @export
fit_qrp <- function(x, cv,
                    plot_title = NULL,
                    start_b = -1,
                    start_c = 0.05,
                    base_size = 14,
                    font_family = "sans",
                    axis_expand_x = 0.05,
                    axis_expand_y = 0.05) {

  call <- match.call()

  # ---- Input validation ----
  if (length(x) != length(cv)) {
    stop("`x` and `cv` must have the same length.", call. = FALSE)
  }
  if (anyNA(x) || anyNA(cv)) {
    stop("`x` and `cv` cannot contain missing values (NA).", call. = FALSE)
  }
  if (length(x) < 4) {
    stop("At least 4 distinct points are required to fit a 3-parameter model.",
         call. = FALSE)
  }
  if (!requireNamespace("minpack.lm", quietly = TRUE)) {
    stop("Package `minpack.lm` is required for `fit_qrp()`. Install it with install.packages('minpack.lm').",
         call. = FALSE)
  }

  # ---- QRP model definition ----
  qrp_model <- function(x, a, b, c) {
    breakpoint <- -b / (2 * c)
    plateau <- a - (b^2 / (4 * c))
    ifelse(x <= breakpoint, a + b * x + c * x^2, plateau)
  }

  obs_data <- data.frame(x = x, cv = cv)

  # ---- Model fitting (unchanged: single nlsLM call) ----
  fit <- tryCatch(
    minpack.lm::nlsLM(
      cv ~ qrp_model(x, a, b, c),
      data = obs_data,
      start = list(a = max(cv), b = start_b, c = start_c)
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    stop(sprintf(
      "Model did not converge. Initial values used: a = %.3f, b = %.3f, c = %.3f. Original error: %s",
      max(cv), start_b, start_c, conditionMessage(fit)
    ), call. = FALSE)
  }

  # ---- Fit statistics ----
  coefs <- coef(fit)
  a <- coefs["a"]; b <- coefs["b"]; c_par <- coefs["c"]
  breakpoint <- -b / (2 * c_par)
  breakpoint_response <- a - (b^2 / (4 * c_par))

  pred <- predict(fit)
  resid <- cv - pred
  n <- length(cv)
  p <- length(coefs)

  sse <- sum(resid^2)
  mse <- mean(resid^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(resid))
  r2 <- 1 - sse / sum((cv - mean(cv))^2)
  r2_adj <- 1 - (1 - r2) * (n - 1) / (n - p)
  aic_val <- AIC(fit)
  bic_val <- BIC(fit)

  # ---- Fitted curve (for plotting) ----
  curve_data <- data.frame(x = seq(min(x), max(x), length.out = 500))
  curve_data$cv <- qrp_model(curve_data$x, a, b, c_par)

  # ---- Dynamic axis limits (replaces hardcoded 0-45 / 0-25) ----
  x_range <- range(c(x, breakpoint))
  y_range <- range(c(cv, breakpoint_response, curve_data$cv))
  x_margin <- diff(x_range) * axis_expand_x
  y_margin <- diff(y_range) * axis_expand_y

  xlim <- c(max(0, x_range[1] - x_margin), x_range[2] + x_margin)
  ylim <- c(max(0, y_range[1] - y_margin), y_range[2] + y_margin)

  # ---- Formatted equations for annotation ----
  fmt <- function(v) sprintf("%.2f", v)
  eq1 <- paste0("CV(x) == '", fmt(a), "'",
                ifelse(b >= 0, " + '", " - '"), fmt(abs(b)), "'*X",
                ifelse(c_par >= 0, " + '", " - '"), fmt(abs(c_par)),
                "'*X^2 ~~~ 'if' ~~~ X <= '", fmt(breakpoint), "'")
  eq2 <- paste0("CV(x) == '", fmt(breakpoint_response), "' ~~~ 'if' ~~~ X > '", fmt(breakpoint), "'")
  eq3 <- paste0("R^2 == '", fmt(r2), "'")
  label_breakpoint <- paste0("Xo == '", fmt(breakpoint), "'")
  label_response <- paste0("CVxo == '", fmt(breakpoint_response), "'")

  x_eq <- mean(xlim)
  y_eqs <- ylim[2] - diff(ylim) * c(0.05, 0.17, 0.28, 0.38, 0.48)

  # ---- Plot ----
  g <- ggplot2::ggplot(obs_data, ggplot2::aes(x, cv)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(data = curve_data, ggplot2::aes(x, cv), linewidth = 1) +
    ggplot2::geom_segment(ggplot2::aes(x = breakpoint, xend = max(xlim),
                                       y = breakpoint_response, yend = breakpoint_response),
                          linewidth = 1) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = breakpoint,
                                       y = breakpoint_response, yend = breakpoint_response),
                          linetype = 5) +
    ggplot2::geom_segment(ggplot2::aes(x = breakpoint, xend = breakpoint,
                                       y = 0, yend = breakpoint_response),
                          linetype = 5) +
    ggplot2::geom_point(ggplot2::aes(x = breakpoint, y = breakpoint_response),
                        colour = "red", size = 3) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[1], label = eq1,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[2], label = eq2,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[3], label = eq3,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[4], label = label_breakpoint,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[5], label = label_response,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    ggplot2::labs(title = plot_title,
                  x = expression("Plot size (" * m^2 * ")"),
                  y = "CV (%)") +
    ggplot2::theme_classic(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.line = ggplot2::element_line(colour = "black"),
      axis.ticks = ggplot2::element_line(colour = "black")
    )

  # ---- Structured output (harmonized with fit_lrp) ----
  result <- list(
    call = call,
    model = fit,
    coefficients = summary(fit)$coefficients,
    parameters = c(
      Breakpoint = unname(breakpoint),
      Breakpoint_Response = unname(breakpoint_response),
      R2 = r2,
      R2_adj = r2_adj,
      RMSE = rmse,
      MAE = mae,
      AIC = aic_val,
      BIC = bic_val,
      SSE = sse,
      MSE = mse
    ),
    data = obs_data,
    curve = curve_data,
    plot = g
  )

  class(result) <- "qrp_fit"
  result
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
  invisible(x)
}

#' @export
summary.qrp_fit <- function(object, ...) {
  cat("Model coefficients:\n")
  print(object$coefficients)
  cat("\nGoodness of fit:\n")
  print(object$parameters)
  invisible(object)
}

#' @export
plot.qrp_fit <- function(x, ...) {
  print(x$plot)
  invisible(x$plot)
}
