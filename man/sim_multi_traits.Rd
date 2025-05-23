% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_linreg.R
\name{sim_multi_traits}
\alias{sim_multi_traits}
\title{Simulate Multiple Traits from Genotype and Effect Sizes}
\usage{
sim_multi_traits(
  G,
  B,
  h2g,
  is_h2g_total = TRUE,
  max_h2g = 1,
  residual_corr = NULL
)
}
\arguments{
\item{G}{Genotype matrix.}

\item{B}{Matrix of effect sizes for multiple traits.}

\item{h2g}{Heritability (proportion of variance explained by genetics).}

\item{is_h2g_total}{Logical indicating if h2g is total (TRUE) or per-SNP (FALSE).}

\item{residual_corr}{Matrix of residual correlations (NULL for independent samples).}
}
\value{
A list containing the simulated phenotypes matrix (t * n, t = trait number, n = sample size) (`P`) and residual variance (`residual_var`).
}
\description{
This function simulates multiple traits (phenotypes) based on genotype data, 
a matrix of effect sizes, and heritability. It allows specifying if the heritability
is total or per-SNP, optionally scales the phenotypes, and can handle residual correlations.
per eQTL heritability is 0.05 to 0.07 according to https://www.nature.com/articles/ng.3506
}
\examples{
G = matrix(rbinom(1000, 2, 0.5), nrow = 1000, ncol = 50) 
# Simulating effect sizes for two traits
B = sim_beta(G, ncausal = 5, ntrait = 3, is_h2g_total = F, shared_pattern = "all")
P = sim_multi_traits(G, B, h2g = 0.1, is_h2g_total = T, max_h2g = 1)
}
