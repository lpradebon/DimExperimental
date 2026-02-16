#' Tamanho ótimo de parcela pelo modelo linear de resposta com platô (LRP)
#'
#' Estima o tamanho ótimo de parcela utilizando o modelo linear de resposta
#' com platô, proposto por Paranaíba, Ferreira e Morais (2009), a partir do
#' agrupamento de parcelas de diferentes tamanhos e formas e de seus
#' respectivos coeficientes de variação (CV).
#'
#' O modelo é composto por dois segmentos: o primeiro descreve uma relação
#' linear crescente ou decrescente entre o coeficiente de variação e o
#' tamanho da parcela até um determinado ponto de quebra (Xo). A partir
#' desse ponto, observa-se um platô, no qual o coeficiente de variação
#' assume um valor constante.
#'
#' @param data \code{data.frame} contendo as variáveis de tamanho da parcela
#'   e coeficiente de variação.
#' @param x_col Nome da coluna do \code{data.frame} que representa o tamanho
#'   da parcela (padrão: \code{"x"}).
#' @param cv_col Nome da coluna do \code{data.frame} que representa o coeficiente
#'   de variação (padrão: \code{"cv"}).
#' @param plot se \code{TRUE}, gera o gráfico do modelo ajustado.
#'
#' @return Uma lista contendo:
#' \describe{
#'   \item{modelo}{Objeto do ajuste não linear (\code{nls}).}
#'   \item{parametros}{\code{data.frame} com os parâmetros estimados do modelo,
#'   incluindo o ponto ótimo de parcela (Xo) e o coeficiente de determinação (R²).}
#'   \item{predicoes}{\code{data.frame} com valores observados e preditos de CV.}
#' }
#'
#' @references
#' Paranaíba, P.F.; Ferreira, D. F.; Morais, A. R. (2009) Tamanho ótimo de parcelas
#' experimentais: proposição de métodos de estimaçãoo. Revista Brasileira de Biometria 27:255–268.
#'
#' @examples
#' dados <- data.frame(
#'   x = c(1,2,4,8,2,4,8,16,3,6,12,24,6,12,24),
#'   cv = c(18.4442670112, 13.7417134862, 10.3923088782, 9.4085177252,
#'   13.4948645984, 10.4858272400, 6.6816059018, 6.6972274507,
#'   10.0527295678, 8.0042609703, 11.9884481410, 5.3277132704,
#'   6.1318167805, 3.8125811769, 0.4788116057))
#'
#'   a <- LRP(dados)
#'
#'   #visualizar os resultados
#'   a
#'   a$modelo
#'
#' @export


LRP <- function(data, x_col = "x", cv_col = "cv", plot = TRUE) {

  # ---------------------------
  # Extrair variáveis
  # ---------------------------
  X  <- data[[x_col]]
  cv <- data[[cv_col]]

  # ---------------------------
  # Modelo linear-platô
  # ---------------------------
  linear_plateau <- function(x, b0, b1, Xo) {
    ifelse(x < Xo, b0 + b1 * x, b0 + b1 * Xo)
  }

  # ---------------------------
  # Ajuste do modelo
  # ---------------------------
  fit <- nls(cv ~ linear_plateau(X, b0, b1, Xo),
             start = list(b0 = max(cv), b1 = -1, Xo = median(X)),
             algorithm = "port",
             lower = c(-Inf, -Inf, min(X)),
             upper = c(Inf, 0, max(X)))

  # ---------------------------
  # Parâmetros estimados
  # ---------------------------
  coef_fit <- coef(fit)
  b0 <- coef_fit["b0"]
  b1 <- coef_fit["b1"]
  Xo <- coef_fit["Xo"]

  CV_plato <- linear_plateau(Xo, b0, b1, Xo)

  # ---------------------------
  # Predição e R²
  # ---------------------------
  cv_pred <- predict(fit)
  R2 <- 1 - sum((cv - cv_pred)^2) / sum((cv - mean(cv))^2)

  # ---------------------------
  # Tabela organizada de parâmetros
  # ---------------------------
  tabela_parametros <- data.frame(
    Parametro = c("Intercepto (b0)",
                  "Inclinação (b1)",
                  "Ponto de quebra (Xo)",
                  "CV no Platô",
                  "R² do modelo"),
    Valor = c(b0, b1, Xo, CV_plato, R2)
  )

  # ---------------------------
  # Tabela de predições
  # ---------------------------
  tabela_predicoes <- data.frame(
    Tamanho_Parcela = X,
    CV_Observado = cv,
    CV_Predito = cv_pred
  )

  # ---------------------------
  # Gráfico (opcional)
  # ---------------------------
  if (plot) {

    plot(X, cv, pch = 19, col = "blue",
         xlab = "Tamanho da parcela (m²)",
         ylab = "CV (%)",
         xlim = c(0, max(X) + 5),
         ylim = c(0, max(cv) + 6),
         xaxt = "n", yaxt = "n",
         main = "MODELO LINEAR COM RESPOSTA PLATÔ\n(Paranaíba et al., 2009)")

    axis(2, at = seq(0, max(cv) + 6, by = 2))
    axis(1, at = seq(0, max(X) + 5, by = 2))

    # Reta linear
    x_reg <- seq(0, Xo, length.out = 200)
    lines(x_reg, b0 + b1 * x_reg, lwd = 2)

    # Platô
    segments(Xo, CV_plato, max(X) + 5, CV_plato, col = "black", lwd = 2)

    # Ponto ótimo
    points(Xo, CV_plato, pch = 19, col = "red")

    # Linhas auxiliares
    segments(0, CV_plato, Xo, CV_plato, col = "red", lty = 2)
    segments(Xo, 0, Xo, CV_plato, col = "red", lty = 2)

    # Equação no gráfico
    eq_text <- paste0(
      "CV(X) = ", round(b0, 2), " - ", round(abs(b1), 2), "X   se X ≤ ", round(Xo, 2), "\n",
      "CV(X) = ", round(CV_plato, 2), "   se X > ", round(Xo, 2), "\n",
      "R² = ", round(R2,2)
    )

    lim <- par("usr")
    text(x = mean(lim[1:2]), y = lim[4] * 0.88, labels = eq_text, cex = 1.05)
  }

  # ---------------------------
  # Retorno organizado
  # ---------------------------
  return(list(
    modelo = fit,
    parametros = tabela_parametros,
    predicoes = tabela_predicoes
  ))
}
