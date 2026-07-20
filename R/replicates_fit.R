#' Estimate the optimal number of replications
#'
#' @description
#' Estimates the optimal number of replications for agricultural
#' experiments, following the methodology described by Cargnelutti Filho
#' et al. (2014), considering different numbers of treatments, experimental
#' coefficient of variation, least significant difference (LSD), and
#' significance level, based on the critical value of the Tukey test. The
#' calculation is performed iteratively for completely randomized designs
#' (CRD) and randomized complete block designs (RCBD).
#'
#' @param treatments Numeric vector. Number of treatments evaluated in the
#'   experiment.
#' @param cv_percent Numeric. Experimental coefficient of variation (%),
#'   typically obtained from the optimal plot size (Xo) estimated from
#'   uniformity trials (see \code{fit_lrp}/\code{fit_qrp}).
#' @param lsd_percent Numeric vector. Least significant difference (LSD),
#'   expressed as a percentage of the experimental mean.
#' @param alpha Numeric, optional. Significance level used in the Tukey
#'   test. Default \code{0.05}.
#' @param design Character, optional. Experimental design, either
#'   \code{"CRD"} (completely randomized design) or \code{"RCBD"}
#'   (randomized complete block design).
#' @param r_init Numeric, optional. Initial value for the number of
#'   replications used in the iterative process. Default \code{2}.
#' @param tol Numeric, optional. Numerical tolerance used as the
#'   convergence criterion between successive iterations. Default \code{1e-10}.
#' @param max_iter Integer, optional. Maximum number of iterations allowed
#'   in the iterative procedure. Default \code{1000}.
#' @param plot Logical, optional. If \code{TRUE}, builds diagnostic plots
#'   relating the optimal number of replications to the number of
#'   treatments and to the Tukey critical value. Default \code{TRUE}. The
#'   plot is stored in the returned object, not printed automatically.
#'
#' @details
#' The optimal number of replications is obtained from:
#' \deqn{r = \left(\frac{q_{\alpha} \cdot CV}{LSD}\right)^2}
#' where \eqn{q_{\alpha}} is the Tukey critical value, \eqn{CV} is the
#' experimental coefficient of variation (%), and \eqn{LSD} is the least
#' significant difference (%). Since \eqn{q_{\alpha}} depends on the error
#' degrees of freedom, which in turn depend on the number of replications,
#' the procedure is solved iteratively until convergence, following
#' Cargnelutti Filho et al. (2014).
#'
#' If a given combination of treatments/LSD does not converge within
#' \code{max_iter} iterations, a warning is issued and the corresponding
#' row is flagged in the \code{converged} column of the result.
#'
#' @references
#' Cargnelutti Filho, A., Alves, B. M., Toebe, M., Burin, C., Santos, G. O.,
#' Facco, G., Neu, I. M. M., & Stefanello, R. B. (2014). Tamanho de parcela e
#' numero de repeticoes em aveia preta. Ciencia Rural, 44(10), 1732-1739.
#'
#' @return An object of class \code{replicates_fit}, a list with:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{data}{A data frame with, for each combination of number of
#'     treatments and LSD: \code{Treatments}, \code{CV_percent},
#'     \code{LSD_percent}, \code{Alpha}, \code{Design}, \code{r_continuous},
#'     \code{r_optimal}, \code{df_error}, \code{q_tukey}, \code{converged}.}
#'   \item{plot}{A combined ggplot2/patchwork object, or \code{NULL} if
#'     \code{plot = FALSE}. Not printed automatically.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- calc_replicates(
#'   treatments  = 3:50,
#'   cv_percent  = 8.97,
#'   lsd_percent = c(10, 20, 30, 40, 50),
#'   design      = "CRD"
#' )
#' fit
#' plot(fit)
#' }
#'
#' @export
calc_replicates <- function(treatments,
                            cv_percent,
                            lsd_percent,
                            alpha = 0.05,
                            design = c("CRD", "RCBD"),
                            r_init = 2,
                            tol = 1e-10,
                            max_iter = 1000,
                            plot = TRUE) {

  call <- match.call()
  design <- match.arg(design)

  # ---- Input validation ----
  if (any(treatments < 2)) {
    stop("`treatments` must be 2 or greater.", call. = FALSE)
  }
  if (any(cv_percent <= 0)) {
    stop("`cv_percent` must be strictly positive.", call. = FALSE)
  }
  if (any(lsd_percent <= 0)) {
    stop("`lsd_percent` must be strictly positive.", call. = FALSE)
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be between 0 and 1.", call. = FALSE)
  }

  # ---- Iterative calculation for one (n_treat, lsd) combination ----
  calc_r_iterative <- function(n_treat, lsd) {

    r_old <- r_init
    r_new <- NA_real_
    converged <- FALSE

    for (iter in seq_len(max_iter)) {

      if (design == "CRD") {
        df_error <- n_treat * (r_old - 1)
      } else {
        df_error <- (n_treat - 1) * (r_old - 1)
      }

      # df_error too low breaks qtukey; force r_old up
      if (df_error <= 1) {
        r_old <- r_old + 1
        next
      }

      q_alpha <- stats::qtukey(1 - alpha, nmeans = n_treat, df = df_error)

      if (is.nan(q_alpha)) {
        r_old <- r_old + 1
        next
      }

      r_new <- ((q_alpha * cv_percent) / lsd)^2

      if (is.nan(r_new)) {
        r_old <- r_old + 1
        next
      }

      if (abs(r_new - r_old) < tol) {
        converged <- TRUE
        break
      }

      r_old <- (r_old + r_new) / 2
    }

    if (is.nan(r_new)) r_new <- r_old

    r_optimal <- ceiling(r_new)

    if (design == "CRD") {
      df_error_final <- n_treat * (r_optimal - 1)
    } else {
      df_error_final <- (n_treat - 1) * (r_optimal - 1)
    }

    q_final <- stats::qtukey(1 - alpha, nmeans = n_treat, df = df_error_final)

    data.frame(
      Treatments  = n_treat,
      CV_percent  = cv_percent,
      LSD_percent = lsd,
      Alpha       = alpha,
      Design      = design,
      r_continuous = r_new,
      r_optimal   = r_optimal,
      df_error    = df_error_final,
      q_tukey     = q_final,
      converged   = converged
    )
  }

  # ---- Full calculation ----
  result_data <- do.call(
    rbind,
    lapply(lsd_percent, function(lsd) {
      do.call(rbind, lapply(treatments, calc_r_iterative, lsd = lsd))
    })
  )

  if (any(!result_data$converged)) {
    n_failed <- sum(!result_data$converged)
    warning(sprintf(
      "%d combination(s) of treatments/LSD did not converge within max_iter = %d iterations. See the `converged` column.",
      n_failed, max_iter
    ), call. = FALSE)
  }

  # ---- Plot ----
  g_combined <- NULL

  if (plot) {

    result_data$LSD_percent <- factor(result_data$LSD_percent)

    p1 <- ggplot2::ggplot(result_data,
                          ggplot2::aes(Treatments, r_optimal, colour = LSD_percent)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
        x = "Number of treatments",
        y = "Optimal number of replications",
        colour = "LSD (%)",
        title = "Optimal replications by number of treatments"
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    p2 <- ggplot2::ggplot(result_data,
                          ggplot2::aes(q_tukey, r_optimal, colour = LSD_percent)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
        x = "Tukey critical value (q)",
        y = "Optimal number of replications",
        colour = "LSD (%)",
        title = "Optimal replications by Tukey critical value"
      ) +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    if (requireNamespace("patchwork", quietly = TRUE)) {
      g_combined <- patchwork::wrap_plots(p1, p2, guides = "collect")
    } else {
      warning("Package `patchwork` not installed; returning the two plots as a list instead of a combined figure.",
              call. = FALSE)
      g_combined <- list(treatments_plot = p1, tukey_plot = p2)
    }
  }

  # ---- Structured output (harmonized with fit_lrp / fit_qrp) ----
  result <- list(
    call = call,
    data = result_data,
    plot = g_combined
  )

  class(result) <- "replicates_fit"
  result
}

#' @export
print.replicates_fit <- function(x, ...) {
  cat("Optimal number of replications\n")
  cat("Rows:", nrow(x$data), " | Non-converged:", sum(!x$data$converged), "\n\n")
  print(utils::head(x$data))
  invisible(x)
}

#' @export
summary.replicates_fit <- function(object, ...) {
  cat("Summary by LSD (%):\n")
  print(stats::aggregate(r_optimal ~ LSD_percent, data = object$data,
                         FUN = function(v) c(min = min(v), max = max(v), mean = mean(v))))
  invisible(object)
}

#' @export
plot.replicates_fit <- function(x, ...) {
  if (is.null(x$plot)) {
    stop("This object was created with `plot = FALSE`; no plot available.", call. = FALSE)
  }
  print(x$plot)
  invisible(x$plot)
}
