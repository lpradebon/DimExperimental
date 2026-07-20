#' Fit the Modified Maximum Curvature (MCM) plot-size model
#'
#' @description
#' Estimates the optimal plot size using the Modified Maximum Curvature
#' method proposed by Meier & Lessmann (1971), from a set of grouped basic
#' experimental units (\code{x}) and their corresponding coefficients of
#' variation (\code{cv}), obtained from uniformity trials.
#'
#' The method identifies the point beyond which increasing the plot size
#' yields negligible gains in experimental precision, and is widely used in
#' experimental sizing studies.
#'
#' @param x Numeric vector. Grouped basic experimental unit sizes (Xo),
#'   obtained from uniformity trials.
#' @param cv Numeric vector. Coefficient of variation of the grouped plots
#'   (%), same length as \code{x}.
#' @param plot_title Character string used as the plot title. Default
#'   \code{"Meier & Lessmann (1971)"}.
#' @param df_correction Logical. If \code{TRUE}, applies a degrees-of-freedom
#'   correction to the breakpoint estimate (some authors apply this
#'   correction, others do not). Default \code{FALSE}, which reproduces the
#'   original uncorrected Meier & Lessmann formula. See Details: the
#'   correction itself is not yet implemented, only the insertion point is
#'   marked in the source.
#' @param base_size Base font size for the plot theme. Default \code{14}.
#' @param font_family Font family for all plot text. Default \code{"sans"}.
#' @param axis_expand_x,axis_expand_y Fractional margin added beyond the data
#'   range when building axis limits. Default \code{0.05} (5%) each.
#'
#' @details
#' The fitting method is unchanged from the original implementation: initial
#' values from a log-log linear regression, followed by a single \code{nls}
#' call on \code{cv ~ a * x^(-b)}.
#'
#' \strong{Degrees-of-freedom correction}: when \code{df_correction = TRUE},
#' the function currently still returns the uncorrected breakpoint. The
#' insertion point for the correction is marked with a \code{TODO} comment
#' in the source, inside \code{fit_mcm()}, right after the uncorrected
#' \code{breakpoint} is computed.
#'
#' @references
#' Meier, V. D.; Lessman, K. J. (1971). Estimation of optimum field plot
#' shape and size for testing yield in Crambe abyssinica hochst. Crop
#' Science, v. 11, p. 648-650.
#'
#' @return An object of class \code{mcm_fit}, a list with:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{model}{The fitted \code{nls} object.}
#'   \item{coefficients}{Coefficient table from \code{summary(model)}.}
#'   \item{parameters}{Named vector: Breakpoint, Breakpoint_Response, R2.}
#'   \item{data}{Original \code{x}/\code{cv} data as a data frame.}
#'   \item{curve}{Fitted curve, 500 points, for plotting.}
#'   \item{plot}{A ggplot2 object. Not printed automatically.}
#' }
#'
#' @examples
#' \dontrun{
#' grouped_data <- data.frame(
#'   x  = c(1, 2, 4, 8, 2, 4, 8, 16, 3, 6, 12, 24, 6, 12, 24),
#'   cv = c(19.55309038, 15.23512639, 10.0533745, 8.203086387,
#'          14.021032, 12.02564594, 7.707005924, 7.207523603,
#'          11.51373415, 9.876254603, 11.59399511, 6.690432438,
#'          7.504701872, 5.791576016, 3.084453732)
#' )
#'
#' # NOTE: the model needs the full set of grouped points, not their means.
#' fit <- fit_mcm(x = grouped_data$x, cv = grouped_data$cv)
#' fit
#' plot(fit)
#' }
#'
#' @export
fit_mcm <- function(x, cv,
                    plot_title = "Meier & Lessmann (1971)",
                    df_correction = FALSE,
                    base_size = 14,
                    font_family = "sans",
                    axis_expand_x = 0.05,
                    axis_expand_y = 0.05) {

  call <- match.call()

  # ---- Input validation ----
  if (length(x) != length(cv)) {
    stop("`x` and `cv` must have the same length.", call. = FALSE)
  }
  if (any(x <= 0) || any(cv <= 0)) {
    stop("`x` and `cv` must be strictly positive.", call. = FALSE)
  }
  if (length(x) < 3) {
    stop("At least 3 points are required to fit a 2-parameter model with residual degrees of freedom.",
         call. = FALSE)
  }

  # ---- Initial values (log-log linear regression, unchanged) ----
  ini <- lm(log(cv) ~ log(x))
  a_ini <- exp(coef(ini)[1])
  b_ini <- -coef(ini)[2]

  # ---- Model fitting (unchanged: single nls call) ----
  fit <- tryCatch(
    nls(cv ~ a * x^(-b), start = list(a = a_ini, b = b_ini)),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    stop(sprintf(
      "Model did not converge. Initial values used: a = %.3f, b = %.3f. Original error: %s",
      a_ini, b_ini, conditionMessage(fit)
    ), call. = FALSE)
  }

  coefs <- coef(fit)
  a <- coefs["a"]; b <- coefs["b"]

  # ---- Modified Maximum Curvature (uncorrected, original formula) ----
  breakpoint <- ((a^2 * b^2 * (2 * b + 1)) / (b + 2))^(1 / (2 * b + 2))

  # ---- Optional degrees-of-freedom correction ----
  # Some authors apply a correction to the breakpoint above based on the
  # error degrees of freedom of the uniformity trial; others use the
  # uncorrected Meier & Lessmann formula as-is (the default here).
  #
  # TODO (Willyan): implement the df correction below when
  # `df_correction = TRUE`. Replace `breakpoint` (and, if the correction
  # also affects it, `breakpoint_response` further down) with the
  # corrected formula.
  if (isTRUE(df_correction)) {
    # Placeholder: no correction implemented yet, `breakpoint` unchanged.
  }

  breakpoint_response <- a * breakpoint^(-b)

  # ---- Goodness of fit ----
  pred <- predict(fit)
  r2 <- 1 - sum((cv - pred)^2) / sum((cv - mean(cv))^2)

  # ---- Fitted curve (for plotting) ----
  curve_data <- data.frame(x = seq(min(x), max(x), length.out = 500))
  curve_data$cv <- a * curve_data$x^(-b)

  obs_data <- data.frame(x = x, cv = cv)

  # ---- Dynamic axis limits ----
  x_range <- range(c(x, breakpoint))
  y_range <- range(c(cv, breakpoint_response, curve_data$cv))
  x_margin <- diff(x_range) * axis_expand_x
  y_margin <- diff(y_range) * axis_expand_y

  xlim <- c(max(0, x_range[1] - x_margin), x_range[2] + x_margin)
  ylim <- c(max(0, y_range[1] - y_margin), y_range[2] + y_margin)

  # ---- Formatted equation for annotation ----
  fmt <- function(v) sprintf("%.2f", v)
  eq1 <- paste0("CV(x) == '", fmt(a), "'*X^{-'", fmt(b), "'}")
  eq2 <- paste0("R^2 == '", fmt(r2), "'")
  label_breakpoint <- paste0("Xo == '", fmt(breakpoint), "'")
  label_response <- paste0("CVxo == '", fmt(breakpoint_response), "'")

  x_eq <- mean(xlim)
  y_eqs <- ylim[2] - diff(ylim) * c(0.05, 0.16, 0.27, 0.38)

  # ---- Plot ----
  g <- ggplot2::ggplot(obs_data, ggplot2::aes(x, cv)) +
    ggplot2::geom_point(size = 3, colour = "darkgreen") +
    ggplot2::geom_line(data = curve_data, ggplot2::aes(x, cv), linewidth = 1, colour = "blue") +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = breakpoint,
                                       y = breakpoint_response, yend = breakpoint_response),
                          linetype = 2, colour = "red") +
    ggplot2::geom_segment(ggplot2::aes(x = breakpoint, xend = breakpoint,
                                       y = 0, yend = breakpoint_response),
                          linetype = 2, colour = "red") +
    ggplot2::geom_point(ggplot2::aes(x = breakpoint, y = breakpoint_response),
                        colour = "red", size = 3) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[1], label = eq1,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[2], label = eq2,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[3], label = label_breakpoint,
                      parse = TRUE, hjust = 0.5, family = font_family, size = base_size / 3.5) +
    ggplot2::annotate("text", x = x_eq, y = y_eqs[4], label = label_response,
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
      axis.text = ggplot2::element_text(colour = "black")
    )

  # ---- Structured output (harmonized with fit_lrp / fit_qrp) ----
  result <- list(
    call = call,
    model = fit,
    coefficients = summary(fit)$coefficients,
    parameters = c(
      Breakpoint = unname(breakpoint),
      Breakpoint_Response = unname(breakpoint_response),
      R2 = r2
    ),
    data = obs_data,
    curve = curve_data,
    plot = g
  )

  class(result) <- "mcm_fit"
  result
}

#' @export
print.mcm_fit <- function(x, ...) {
  cat("Modified Maximum Curvature (MCM) fit\n")
  cat("Breakpoint (Xo):        ", sprintf("%.3f", x$parameters["Breakpoint"]), "\n")
  cat("CV at breakpoint:       ", sprintf("%.3f", x$parameters["Breakpoint_Response"]), "\n")
  cat("R2:", sprintf("%.3f", x$parameters["R2"]), "\n")
  invisible(x)
}

#' @export
summary.mcm_fit <- function(object, ...) {
  cat("Model coefficients:\n")
  print(object$coefficients)
  cat("\nGoodness of fit:\n")
  print(object$parameters)
  invisible(object)
}

#' @export
plot.mcm_fit <- function(x, ...) {
  print(x$plot)
  invisible(x$plot)
}
