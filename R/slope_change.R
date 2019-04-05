#' This function fits to a numerical vector sorted in the non decreasing order two simple linear regressions and returns the index corresponding to the estimated change between the two regression models.
#'
#' @param  Y numerical vector sorted in the non decreasing order.
#' @return K the index corresponding to the estimated change between the two linear regression models.
#' @importFrom Matrix Matrix
#' @importFrom dplyr arrange filter mutate cummean
#' @importFrom tibble tibble rowid_to_column
#' @examples
#' n <- 30
#' q <- 100
#' Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)
#' Matrix::image(Sigma)
#' E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))
#' corE <- cor(as.matrix(E))
#' vec_up_emp <- corE[upper.tri(corE)]
#' G <- matrix(0, ncol = (q - 1), nrow = (q - 1))
#' G[upper.tri(G, diag = TRUE)] <- vec_up_emp
#' G[lower.tri(G)] <- t(as.matrix(G))[lower.tri(t(as.matrix(G)))]
#' res_svd <- svd(G)
#' vp <- res_svd$d
#' slope_change(vp)
#' @export
slope_change <- function(Y) {
  tb <- tibble(y = Y) %>%
    arrange(desc(y)) %>%
    rowid_to_column(var = "x") %>%
    mutate(
      mx = (x + 1) / 2,
      my = cummean(y),
      pxy = x * y,
      spxy = cumsum(pxy),
      s2xy = spxy - mx * x * my,
      s2x = (x * (x + 1) * (x - 1)) / 12,
      b = s2xy / s2x,
      a = my - b * mx,
      y2 = cumsum(y^2),
      e = y2 + (b^2 * (x * (x + 1) * (2 * x + 1)) / 6) +
        2 * a * b * x * mx + x * a^2 -
        2 * b * spxy - 2 * a * my * x
    ) %>%
    filter(x != n())

  n <- length(Y)
  tb2 <- tibble(y = Y) %>%
    arrange(desc(y)) %>%
    rowid_to_column(var = "x") %>%
    arrange(y) %>%
    rowid_to_column(var = "k") %>%
    mutate(
      mx = (x + n) / 2,
      my = cummean(y),
      pxy = x * y,
      spxy = cumsum(pxy),
      s2xy = spxy - mx * k * my,
      s2xi = (n * (n + 1) * (2 * n + 1) - x * (x - 1) * (2 * x - 1)) / 6,
      s2x = s2xi - k * (n + x)^2 / 4,
      b = s2xy / s2x,
      a = my - b * mx,
      y2 = cumsum(y^2),
      e = y2 + (b^2 * s2xi) +
        2 * a * b * k * mx + k * a^2 -
        2 * b * spxy - 2 * a * my * k
    ) %>%
    filter(x != n())
  errors <- tb$e + rev(tb2$e)
  which.min(errors) - 1
}
