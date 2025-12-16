#'Tamanho de parcela
#'@description
#'Dimensiona o tamanho ótiomo de parcela para experimentos agrícolas pelo
#'método da curvatura máxima em função do coeficiente de variação
#' (PARANAÍBA; FERREIRA; MORAIS, 2009)
#' @param dados data.frame ou tibble contendo os dados brutos do ensaio de uniformidade.
#' Cada linha deve representar uma unidade básica experimental (UEB).
#' As colunas iniciais devem conter informações de identificação espacial
#' (por exemplo, linha, coluna e repetição), seguidas pelas variáveis
#' quantitativas de interesse.
#'
#' @param nlin inteiro. Número de linhas do ensaio (direção longitudinal),
#' correspondente ao número de unidades básicas no sentido das linhas de cultivo.
#'
#' @param ncol inteiro. Número de colunas do ensaio (direção transversal),
#' correspondente ao número de unidades básicas no sentido da largura do experimento.
#'
#' @param nrep inteiro. Número de repetições simuladas ou divisões do ensaio,
#' utilizado para o cálculo das estatísticas de variabilidade entre parcelas.
#'
#' @param col_inicio inteiro, opcional. Índica a coluna em que se iniciam
#' as variáveis quantitativas de interesse em \code{dados}.
#' As colunas anteriores são tratadas como metadados.
#' O valor padrão é \code{4}.
#'
#' @param digits inteiro, opcional. Número de casas decimais utilizadas
#' para o arredondamento dos resultados apresentados.
#' O valor padrão é \code{2}.
#'
#' @return Um data.frame contendo as estatísticas do ensaio de uniformidade,
#' incluindo médias, coeficientes de variação e parâmetros associados
#' ao método de Paranaíba para diferentes tamanhos de parcela.
#'
#' @details
#' O método de Paranaíba baseia-se na relação entre o coeficiente de variação
#' e o tamanho da parcela, permitindo a estimativa do tamanho ótimo de parcela
#' e da variabilidade experimental associada. A função utiliza a estrutura
#' espacial definida por \code{nlin} e \code{ncol} para simular agregações
#' de unidades básicas experimentais.
#'
#' @references
#' Paranaíba, P F; Ferreira, D F; Morais, A R. (2009). Tamanho ótimo de parcelas
#' experimentais: proposição de métodos de estimação. \emph{Revista Brasileira
#' de Biometria}, 27:255–268.
#' @examples
#' # Carregar os dados de exemplo do pacote
#' data(dados_ensaio_C1)
#'
#' # Aplicar o método de Paranaíba
#' resultado_C1 <- calc_paranaiba(
#'   dados      = dados_ensaio_C1,
#'   nlin       = 6,
#'   ncol       = 8,
#'   nrep       = 3,
#'   col_inicio = 4,
#'   digits     = 3
#' )
#'
#' # Visualizar os resultados
#' resultado_C1
#'
#' @export

calc_paranaiba <- function(dados,
                           nlin,
                           ncol,
                           nrep,
                           col_inicio = 4,
                           digits = 2){

  # =====================================================
  # Função interna: cálculo para um bloco (1 repetição)
  # =====================================================
  calcula_bloco <- function(x){

    media  <- mean(x)
    desvio <- sd(x)
    cv     <- (desvio / media) * 100

    err   <- x - mean(x)
    soma2 <- sum(err^2)
    var1  <- desvio^2

    # -----------------------------
    # Caminhamento em LINHA
    # -----------------------------
    erro_linha <- c()
    for(i in 1:nrow(x)){
      if(i %% 2 == 1){
        erro_linha <- c(erro_linha, err[i, ])
      } else {
        erro_linha <- c(erro_linha, err[i, ncol(x):1])
      }
    }

    rho_l <- sum(erro_linha[-1] * erro_linha[-length(erro_linha)]) / soma2

    Xo_l <- 10 * (2 * (1 - rho_l^2) * var1 * media)^(1/3) / media

    # -----------------------------
    # Caminhamento em COLUNA
    # -----------------------------
    erro_col <- c()
    for(j in 1:ncol(x)){
      if(j %% 2 == 1){
        erro_col <- c(erro_col, err[, j])
      } else {
        erro_col <- c(erro_col, err[nrow(x):1, j])
      }
    }

    rho_c <- sum(erro_col[-1] * erro_col[-length(erro_col)]) / soma2

    # -----------------------------
    # Médias e valores ótimos
    # -----------------------------
    rho_med <- mean(c(rho_l, rho_c))

    tamanho_otimo <- 10 * (2 * (1 - rho_med^2) * var1 * media)^(1/3) / media

    cv_otimo <- 100 *
      sqrt((1 - rho_med^2) * var1 / media^2) /
      sqrt(tamanho_otimo)

    # -----------------------------
    # Tabela do bloco
    # -----------------------------
    data.frame(
      Estatistica = c(
        "média",
        "rho_l",
        "rho_c",
        "rho_médio",
        "Tamanho ótimo da parcela",
        "Coeficiente de variação (%)",
        "Coeficiente de variação da parcela ótima (%)"
      ),
      Valor = c(
        media,
        rho_l,
        rho_c,
        rho_med,
        tamanho_otimo,
        cv,
        cv_otimo
      ),
      stringsAsFactors = FALSE
    )
  }

  # =====================================================
  # Loop nas repetições
  # =====================================================
  resultados <- vector("list", nrep)

  for(rep in 1:nrep){

    lin_ini <- (rep - 1) * nlin + 1
    lin_fim <- rep * nlin

    bloco <- matrix(
      unlist(dados[lin_ini:lin_fim,
                   col_inicio:(col_inicio + ncol - 1)]),
      nrow = nlin,
      ncol = ncol,
      byrow = FALSE
    )

    res <- calcula_bloco(bloco)
    res$Repeticao <- rep

    resultados[[rep]] <- res
  }

  # =====================================================
  # Tabela final
  # =====================================================
  tabela_final <- do.call(rbind, resultados)

  tabela_final$Valor <- round(tabela_final$Valor, digits)
  rownames(tabela_final) <- NULL

  return(tabela_final)
}
