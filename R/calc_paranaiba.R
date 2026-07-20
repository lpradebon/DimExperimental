#' Estimate optimal plot size using the Paranaiba method
#'
#' @description
#' Estimates the optimal plot size for agricultural experiments using the
#' maximum curvature method as a function of the coefficient of variation,
#' following Paranaiba, Ferreira & Morais (2009). The method uses the
#' spatial layout of a uniformity trial (rows and columns of basic
#' experimental units) to estimate the spatial correlation between
#' neighboring units, and from that, the optimal plot size.
#'
#' @param data A data frame or tibble with the raw uniformity trial data.
#'   Each row represents one basic experimental unit (BEU). The initial
#'   columns should hold spatial identification (e.g. row, column,
#'   replication), followed by the quantitative variable of interest.
#' @param n_row Integer. Number of rows in the trial (longitudinal
#'   direction), i.e. the number of basic units along the crop rows.
#' @param n_col Integer. Number of columns in the trial (transversal
#'   direction), i.e. the number of basic units along the experiment width.
#' @param n_rep Integer. Number of simulated replications/splits of the
#'   trial used to compute the between-plot variability statistics.
#' @param start_col Integer, optional. Column index in \code{data} where
#'   the quantitative variable of interest starts; preceding columns are
#'   treated as metadata. Default \code{4}.
#' @param digits Integer, optional. Number of decimal places used when
#'   rounding the reported results. Default \code{2}.
#' @param plot_title Character string used as the plot title. Default
#'   \code{NULL}.
#' @param base_size Base font size for the plot theme. Default \code{14}.
#' @param font_family Font family for all plot text. Default \code{"sans"}.
#'
#' @details
#' The Paranaiba method is based on the relationship between the
#' coefficient of variation and plot size, allowing the estimation of the
#' optimal plot size and its associated experimental variability. The
#' function uses the spatial structure defined by \code{n_row} and
#' \code{n_col} to simulate aggregations of basic experimental units,
#' walking a serpentine path across rows and columns to estimate the
#' spatial autocorrelation (\code{rho}).
#'
#' If a given replication has zero variance (a degenerate block), its
#' \code{rho}, optimal size, and optimal CV are returned as \code{NA} and a
#' warning is issued, rather than silently propagating \code{NaN}.
#'
#' @references
#' Paranaiba, P. F.; Ferreira, D. F.; Morais, A. R. (2009). Tamanho otimo de
#' parcelas experimentais: proposicao de metodos de estimacao. Revista
#' Brasileira de Biometria, 27:255-268.
#'
#' @return An object of class \code{paranaiba_fit}, a list with:
#' \describe{
#'   \item{call}{The matched call.}
#'   \item{data}{A data frame with one row per replication: \code{Repetition},
#'     \code{Mean}, \code{Rho_row}, \code{Rho_col}, \code{Rho_mean},
#'     \code{Optimal_size}, \code{CV_percent}, \code{CV_optimal}, \code{valid}.}
#'   \item{plot}{A ggplot2 object showing optimal plot size by replication.
#'     Not printed automatically.}
#' }
#'
#' @examples
#' \dontrun{
#' fit <- calc_paranaiba(
#'   data      = trial_data,
#'   n_row     = 6,
#'   n_col     = 8,
#'   n_rep     = 3,
#'   start_col = 4,
#'   digits    = 3
#' )
#' fit
#' plot(fit)
#' }
#'
#' @export
calc_paranaiba <- function(data,
                           n_row,
                           n_col,
                           n_rep,
                           start_col = 4,
                           digits = 2,
                           plot_title = NULL,
                           base_size = 14,
                           font_family = "sans") {

  call <- match.call()

  # ---- Input validation ----
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }
  if (start_col + n_col - 1 > ncol(data)) {
    stop(sprintf(
      "`start_col` + `n_col` - 1 = %d exceeds the number of columns in `data` (%d).",
      start_col + n_col - 1, ncol(data)
    ), call. = FALSE)
  }
  if (n_rep * n_row > nrow(data)) {
    stop(sprintf(
      "`n_rep` * `n_row` = %d exceeds the number of rows in `data` (%d).",
      n_rep * n_row, nrow(data)
    ), call. = FALSE)
  }

  # ---- Internal calculation for one block (one replication) ----
  calc_block <- function(block, rep_id) {

    block_mean <- mean(block)
    block_sd   <- sd(block)
    cv         <- (block_sd / block_mean) * 100

    err   <- block - block_mean
    ss    <- sum(err^2)
    var1  <- block_sd^2

    if (ss == 0) {
      warning(sprintf("Replication %d has zero variance; rho and optimal size set to NA.", rep_id),
              call. = FALSE)
      return(data.frame(
        Repetition = rep_id, Mean = block_mean, Rho_row = NA_real_,
        Rho_col = NA_real_, Rho_mean = NA_real_, Optimal_size = NA_real_,
        CV_percent = cv, CV_optimal = NA_real_, valid = FALSE
      ))
    }

    # ---- Row-wise serpentine walk ----
    err_row <- c()
    for (i in seq_len(nrow(block))) {
      if (i %% 2 == 1) {
        err_row <- c(err_row, err[i, ])
      } else {
        err_row <- c(err_row, err[i, ncol(block):1])
      }
    }
    rho_row <- sum(err_row[-1] * err_row[-length(err_row)]) / ss

    # ---- Column-wise serpentine walk ----
    err_col <- c()
    for (j in seq_len(ncol(block))) {
      if (j %% 2 == 1) {
        err_col <- c(err_col, err[, j])
      } else {
        err_col <- c(err_col, err[nrow(block):1, j])
      }
    }
    rho_col <- sum(err_col[-1] * err_col[-length(err_col)]) / ss

    rho_mean <- mean(c(rho_row, rho_col))

    optimal_size <- 10 * (2 * (1 - rho_mean^2) * var1 * block_mean)^(1 / 3) / block_mean

    cv_optimal <- 100 *
      sqrt((1 - rho_mean^2) * var1 / block_mean^2) /
      sqrt(optimal_size)

    data.frame(
      Repetition = rep_id, Mean = block_mean, Rho_row = rho_row,
      Rho_col = rho_col, Rho_mean = rho_mean, Optimal_size = optimal_size,
      CV_percent = cv, CV_optimal = cv_optimal, valid = TRUE
    )
  }

  # ---- Loop over replications ----
  results <- vector("list", n_rep)

  for (rep in seq_len(n_rep)) {

    row_start <- (rep - 1) * n_row + 1
    row_end   <- rep * n_row

    block <- matrix(
      unlist(data[row_start:row_end, start_col:(start_col + n_col - 1)]),
      nrow = n_row,
      ncol = n_col,
      byrow = FALSE
    )

    results[[rep]] <- calc_block(block, rep_id = rep)
  }

  result_data <- do.call(rbind, results)
  rownames(result_data) <- NULL

  numeric_cols <- c("Mean", "Rho_row", "Rho_col", "Rho_mean", "Optimal_size",
                    "CV_percent", "CV_optimal")
  result_data[numeric_cols] <- lapply(result_data[numeric_cols], round, digits = digits)

  # ---- Plot: optimal plot size by replication ----
  g <- ggplot2::ggplot(result_data, ggplot2::aes(factor(Repetition), Optimal_size)) +
    ggplot2::geom_col(fill = "steelblue", na.rm = TRUE) +
    ggplot2::labs(
      title = plot_title,
      x = "Replication",
      y = expression("Optimal plot size (" * m^2 * ")")
    ) +
    ggplot2::theme_classic(base_size = base_size, base_family = font_family) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # ---- Structured output ----
  result <- list(
    call = call,
    data = result_data,
    plot = g
  )

  class(result) <- "paranaiba_fit"
  result
}

#' @export
print.paranaiba_fit <- function(x, ...) {
  cat("Paranaiba optimal plot size estimate\n")
  cat("Replications:", nrow(x$data), " | Invalid (zero variance):", sum(!x$data$valid), "\n\n")
  print(x$data)
  invisible(x)
}

#' @export
summary.paranaiba_fit <- function(object, ...) {
  cat("Optimal plot size across replications:\n")
  print(summary(object$data$Optimal_size))
  invisible(object)
}

#' @export
plot.paranaiba_fit <- function(x, ...) {
  print(x$plot)
  invisible(x$plot)
}
