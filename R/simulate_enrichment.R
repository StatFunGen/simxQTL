#' Simulate Causal Configuration
#'
#' This function simulates a causal configuration based on annotations and odds ratios.
#'
#' @param A Annotation matrix where each column is a vector of real numbers.
#' @param odds_ratio Vector of odds ratios for each column in the annotation matrix.
#' @param baseline_odds Baseline odds for calculating probabilities.
#' @param n_samples Number of samples to draw.
#' @return A list with two elements: a vector of causal probabilities and a matrix of binary samples.
#' @examples
#' A <- matrix(runif(100), ncol = 5)
#' odds_ratio <- runif(5, 1, 2)
#' baseline_odds <- 0.5
#' n_samples <- 100
#' result <- simulate_causal_config(A, odds_ratio, baseline_odds, n_samples)
#' @export
simulate_causal_config <- function(A, odds_ratio, baseline_odds, n_samples) {
  if (ncol(A) != length(odds_ratio)) stop("Length of odds_ratio must match the number of columns in A")
  
  # Calculate the product of odds ratios for each variant
  odds_product <- apply(A, 1, function(row) prod(row ^ odds_ratio))
  
  # Convert to causal probabilities
  causal_prob <- baseline_odds * odds_product / (1 + baseline_odds * odds_product)
  
  # Generate binary samples
  samples <- matrix(rbinom(n_samples * nrow(A), 1, causal_prob), nrow = n_samples, ncol = nrow(A))

  return(list(causal_probabilities = causal_prob, samples = samples))
}