#'Número de Repetições
#'@description
#' Função para estimar o número ótimo de repetições em experimentos agrícolas,
#' conforme a metodologia descrita por Cargnelutti Filho et al. (2014),
#' considerando diferentes números de tratamentos, coeficiente de variação,
#' diferença mínima significativa (DMS) e nível de significância, com base
#' no valor crítico do teste de Tukey. O cálculo é realizado de forma iterativa
#' para delineamentos inteiramente casualizado (DIC) e em blocos casualizados (DBC).
#' Os parâmetros podem ser definidos conforme o resquisitos do pesquisador.
#'
#' @param tratamentos vetor numérico. Números de tratamentos avaliados
#' no experimento.
#'
#' @param CV_percent numérico. Coeficiente de variação experimental (%),
#' obtido pelo tamanho ótimo de parcela (Xo) a partir de ensaios de uniformidade.
#'
#' @param d_percent vetor numérico. Diferença mínima significativa (DMS),
#' expressa em porcentagem da média do experimento.
#'
#' @param alpha numérico, opcional. Nível de significância adotado no teste
#' de Tukey. O valor padrão é \code{0.05}.
#'
#' @param design caractere, opcional. Delineamento experimental considerado,
#' podendo ser \code{"DIC"} (delineamento inteiramente casualizado) ou
#' \code{"DBC"} (delineamento em blocos casualizados).
#'
#' @param r_init numérico, opcional. Valor inicial do número de repetições
#' utilizado no processo iterativo. O valor padrão é \code{2}.
#'
#' @param tol numérico, opcional. Tolerância numérica utilizada como critério
#' de convergência entre iterações sucessivas do número de repetições.
#' O valor padrão é \code{1e-10}.
#'
#' @param max_iter inteiro, opcional. Número máximo de iterações permitidas
#' no procedimento iterativo. O valor padrão é \code{1000}.
#'
#' @param plot lógico, opcional. Se \code{TRUE}, gera gráficos ilustrativos
#' relacionando o número ótimo de repetições com o número de tratamentos e
#' com o valor crítico do teste de Tukey. O valor padrão é \code{TRUE}.
#'
#' @return
#' Um \code{data.frame} contendo, para cada combinação de número de tratamentos
#' e DMS:
#' \itemize{
#'   \item \code{Tratamentos}: número de tratamentos, ;
#'   \item \code{CV_percent}: coeficiente de variação experimental (%);
#'   \item \code{DMS_percent}: diferença mínima significativa (%);
#'   \item \code{Alpha}: nível de significância;
#'   \item \code{Delineamento}: tipo de delineamento experimental;
#'   \item \code{r_continuo}: número de repetições estimado em escala contínua;
#'   \item \code{r_otimo}: número ótimo de repetições (valor inteiro);
#'   \item \code{GLE}: graus de liberdade do erro;
#'   \item \code{q_Tukey}: valor crítico do teste de Tukey.
#' }
#'
#' @details
#' O número ótimo de repetições é obtido a partir da relação:
#' \deqn{r = \left(\frac{q_{\alpha} \cdot CV}{DMS}\right)^2}
#' em que \eqn{q_{\alpha}} é o valor crítico do teste de Tukey,
#' \eqn{CV} é o coeficiente de variação experimental (%) e \eqn{DMS}
#' é a diferença mínima significativa (%). Como o valor de
#' \eqn{q_{\alpha}} depende dos graus de liberdade do erro,
#' que por sua vez dependem do número de repetições,
#' o procedimento é resolvido por iteração até convergência,
#' conforme descrito por Cargnelutti Filho et al. (2014).
#'
#' @references
#'Cargnelutti Filho, A., Alves, B. M., Toebe, M., Burin, C., Santos, G. O.,
#'Facco, G., Neu, I. M. M., & Stefanello, R. B. (2014). Tamanho de parcela e
#'número de repetições em aveia preta. Ciência Rural, 44(10), 1732–1739.

#' @examples
#' # ===============================
#' # Dados de entrada
#' # ===============================
#' cv <- 8.96596
#' alpha <- 0.05
#' tratamentos <- 3:50
#' DMS <- c(10, 20, 30, 40, 50)
#'
#' # ===============================
#' # Cálculo para DIC
#' # ===============================
#' resultado_DIC <- calc_repeticoes(
#'   tratamentos = tratamentos,
#'   CV_percent  = cv,
#'   d_percent   = DMS,
#'   alpha       = alpha,
#'   design      = "DIC",
#'   plot        = TRUE
#' )
#'
#' # ===============================
#' # Cálculo para DBC
#' # ===============================
#' resultado_DBC <- calc_repeticoes(
#'   tratamentos = tratamentos,
#'   CV_percent  = cv,
#'   d_percent   = DMS,
#'   alpha       = alpha,
#'   design      = "DBC",
#'   plot        = TRUE
#' )
#'
#' # Visualizar parte dos resultados
#' head(resultado_DIC)
#' head(resultado_DBC)
#'
#' @export

calc_repeticoes <- function(tratamentos,
                            CV_percent,
                            d_percent,
                            alpha = 0.05,
                            design = c("DIC","DBC"),
                            r_init = 2,
                            tol = 1e-10,
                            max_iter = 1000,
                            plot = TRUE){

  design <- match.arg(design)

  calc_r_iterativo <- function(i, dms){

    r_old <- r_init
    r_new <- NA_real_

    for(iter in 1:max_iter){

      if(design == "DIC"){
        GLE <- i * (r_old - 1)
      } else {
        GLE <- (i - 1) * (r_old - 1)
      }

      # GLE muito baixo quebra a qtukey; força aumentar r_old
      if(GLE <= 1){
        r_old <- r_old + 1
        next
      }

      q_alpha <- stats::qtukey(1 - alpha, nmeans = i, df = GLE)

      # Se qtukey não conseguir calcular, tenta aumentar r_old
      if (is.nan(q_alpha)) {
        r_old <- r_old + 1
        next
      }

      r_new <- ((q_alpha * CV_percent) / dms)^2

      # Proteção contra NaN em r_new
      if (is.nan(r_new)) {
        r_old <- r_old + 1
        next
      }

      if(abs(r_new - r_old) < tol) break

      r_old <- (r_old + r_new) / 2
    }

    # Caso extremo: se nunca conseguiu r_new válido, usa r_old
    if (is.nan(r_new)) r_new <- r_old

    r_otimo <- ceiling(r_new)

    if(design == "DIC"){
      GLE_final <- i * (r_otimo - 1)
    } else {
      GLE_final <- (i - 1) * (r_otimo - 1)
    }

    q_final <- stats::qtukey(1 - alpha, nmeans = i, df = GLE_final)

    data.frame(
      Tratamentos  = i,
      CV_percent   = CV_percent,
      DMS_percent  = dms,
      Alpha        = alpha,
      Delineamento = design,
      r_continuo   = r_new,
      r_otimo      = r_otimo,
      GLE          = GLE_final,
      q_Tukey      = q_final
    )
  }

  # ======================
  # Cálculo completo
  # ======================
  resultado <- do.call(
    rbind,
    lapply(d_percent, function(dms){
      do.call(rbind,
              lapply(tratamentos, calc_r_iterativo, dms = dms))
    })
  )

  # ======================
  # GRÁFICOS
  # ======================
  if(plot){

    dms_levels <- unique(resultado$DMS_percent)
    cores <- seq_along(dms_levels)
    names(cores) <- dms_levels

    # Layout: 2 gráficos + legenda
    layout(matrix(c(1, 2, 3), ncol = 1), heights = c(4, 4, 1))

    # ---------- Gráfico 1 ----------
    par(mar = c(4, 5, 3, 2))
    plot(NA,
         xlim = range(tratamentos),
         ylim = range(resultado$r_otimo),
         xlab = "Número de tratamentos",
         ylab = "Número ótimo de repetições",
         main = "Número ótimo de repetições para diferentes DMS")

    grid(col = "gray80", lty = "dotted")

    for(dms in dms_levels){
      dados_dms <- subset(resultado, DMS_percent == dms)
      lines(dados_dms$Tratamentos,
            dados_dms$r_otimo,
            type = "o",
            pch = 16,
            lwd = 2,
            col = cores[as.character(dms)])
    }

    # ---------- Gráfico 2 ----------
    par(mar = c(5, 5, 3, 2))

    # limites seguros
    x_q <- range(resultado$q_Tukey, na.rm = TRUE)
    y_r <- range(resultado$r_otimo, na.rm = TRUE)

    # evita eixo com amplitude zero (caso comum no DBC)
    if(diff(x_q) == 0){
      x_q <- x_q + c(-0.1, 0.1)
    }

    plot(NA,
         xlim = x_q,
         ylim = y_r,
         xlab = "Valor crítico de Tukey (q)",
         ylab = "Número ótimo de repetições",
         main = "Valor crítico de Tukey em função do número ótimo de repetições")

    grid(col = "gray80", lty = "dotted")

    for(dms in dms_levels){
      dados_dms <- subset(resultado, DMS_percent == dms)
      points(dados_dms$q_Tukey,
             dados_dms$r_otimo,
             pch = 16,
             col = cores[as.character(dms)])
    }
    # ---------- Legenda ----------
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center",
           legend = paste("DMS =", dms_levels, "%"),
           col = cores,
           pch = 16,
           lwd = 2,
           horiz = TRUE,
           bty = "n")

    layout(1)
  }

  return(resultado)
}
