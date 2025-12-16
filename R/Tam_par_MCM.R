#' Tamanho ótimo de parcela pelo método da Máxima Curvatura Modificada (MCM)
#'
#' Estima o tamanho ótimo de parcela utilizando o método da Máxima Curvatura
#' Modificada, proposto por Meier & Lessmann (1971), a partir do tamanho ótimo
#' de parcela de referência (X0) e do coeficiente de variação do ensaio.
#'
#' O método baseia-se na identificação do ponto em que o aumento do tamanho
#' da parcela passa a resultar em ganhos irrisórios na precisão experimental,
#' sendo amplamente utilizado em estudos de dimensionamento experimental.
#'
#' @param Xo valor numérico correspondente ao tamanho ótimo de parcela de
#'   referência (X₀), obtidos pelo agrupamento de unidades experimentais basicas, 
#'   obtidas em ensaios de uniformidade.
#' @param CV coeficiente de variação das parcelas agrupadas, expresso em porcentagem (%).
#' @param main título do gráfico gerado pela função.
#' @param plot lógico; se \code{TRUE}, gera o gráfico da relação entre o
#'   tamanho da parcela e o coeficiente de variação estimado.
#'
#' @return Um valor numérico correspondente ao tamanho ótimo de parcela
#'   estimado pelo método da Máxima Curvatura Modificada.
#'
#' @references
#'MEIER, V. D.; LESSMAN, K. J. (1971).Estimation of optimum field plot shape and size 
#'for testing yield in Crambe abyssinica hochst. Crop Science, v. 11, p. 648-650.
#' @examples
#' # ---------------------------------------
#' # Conjunto de dados de exemplo
#' # ---------------------------------------
#' dados_curvatura <- data.frame(
#'   Xo = c(1, 2, 4, 8,
#'          2, 4, 8, 16,
#'          3, 6, 12, 24,
#'          6, 12, 24),
#'   CV = c(19.55309038, 15.23512639, 10.0533745,  8.203086387,
#'          14.021032,   12.02564594, 7.707005924, 7.207523603,
#'          11.51373415, 9.876254603, 11.59399511, 6.690432438,
#'          7.504701872, 5.791576016, 3.084453732)
#' )
#'
#' # Valores médios utilizados no método
#' Xo_med <- mean(dados_curvatura$Xo)
#' CV_med <- mean(dados_curvatura$CV)
#'
#' # Estimativa do tamanho ótimo de parcela
#' cal_max_curv_mod(
#'   Xo = Xo_med,
#'   CV = CV_med,
#'   plot = TRUE
#' )


cal_max_curv_mod <- function(Xo, CV,
                             main = "Meier & Lessmann (1971)",
                             plot = TRUE) {
  
  # =============================
  # Checagens
  # =============================
  if (length(Xo) != length(CV))
    stop("Dimensões incompatíveis!")
  
  if (any(Xo <= 0) || any(CV <= 0))
    stop("Xo e CV devem ser maiores que zero.")
  
  # =============================
  # Estimativas iniciais (log-log)
  # =============================
  ini <- lm(log(CV) ~ log(Xo))
  a_ini <- exp(coef(ini)[1])
  b_ini <- -coef(ini)[2]
  
  # =============================
  # Ajuste não linear
  # =============================
  fit <- nls(
    CV ~ a * Xo^(-b),
    start = list(a = a_ini, b = b_ini)
  )
  
  a <- coef(fit)["a"]
  b <- coef(fit)["b"]
  
  # =============================
  # Máxima Curvatura Modificada
  # =============================
  Xo_mc <- ((a^2 * b^2 * (2*b + 1)) / (b + 2))^(1 / (2*b + 2))
  CV_Xo <- a * Xo_mc^(-b)
  
  # =============================
  # Qualidade do ajuste
  # =============================
  pred <- predict(fit)
  R2 <- 1 - sum((CV - pred)^2) / sum((CV - mean(CV))^2)
  
  # =============================
  # Gráfico
  # =============================
  if (plot) {
    
    X_lim_max <- max(Xo) * 1.05
    Y_lim_max <- max(CV) * 1.10
    
    plot(Xo, CV,
         pch = 19,
         col = "darkgreen",
         xlab = "Tamanho da parcela (m²)",
         ylab = "Coeficiente de variação (%)",
         main = main,
         xaxt = "n",
         yaxt = "n",
         xlim = c(0, X_lim_max),
         ylim = c(0, Y_lim_max))
    
    axis(1, at = pretty(Xo))
    axis(2, at = pretty(CV))
    
    curve(a * x^(-b),
          from = min(Xo) * 1.02,
          to   = max(Xo) * 0.98,
          add = TRUE,
          lwd = 2,
          col = "blue")
    
    points(Xo_mc, CV_Xo, pch = 19, col = "red", cex = 1.4)
    
    segments(Xo_mc, 0, Xo_mc, CV_Xo, col = "red", lty = 2)
    segments(0, CV_Xo, Xo_mc, CV_Xo, col = "red", lty = 2)
    
    legend_text <- bquote(
      atop(
        Y == .(round(a,3)) / x^{.(round(b,3))},
        R^2 == .(round(R2,3)) ~ ~ 
          X[o] == .(round(Xo_mc,2)) ~ ~ 
          CV(X[o]) == .(round(CV_Xo,2))
      )
    )
    
    text(x = X_lim_max * 0.6,
         y = Y_lim_max * 0.9,
         labels = legend_text,
         cex = 0.9)
  }
  
  # =============================
  # Tabela de resultados
  # =============================
  resultados <- data.frame(
    Parametro = c("a", "b", "Xo (MC)", "CV(Xo)", "R²"),
    Valor = c(a, b, Xo_mc, CV_Xo, R2)
  )
  
  print(resultados)
  
  # =============================
  # Retorno
  # =============================
  invisible(list(
    model = fit,
    resultados = resultados
  ))
}