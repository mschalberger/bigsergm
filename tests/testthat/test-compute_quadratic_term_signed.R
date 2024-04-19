test_that("quadratic term calculation without features works", {
  # Number of nodes
  N <- 6
  # Number of clusters
  K <- 3
  # Create an adjacency matrix
  # Create a signed adjacency matrix
  adj <- matrix(sample(c(-1, 0, 1), N^2, replace = TRUE), nrow = N)
  adj[lower.tri(adj)] = t(adj)[lower.tri(adj)]
  diag(adj) <- 0

  positive_adj <- adj
  positive_adj[positive_adj == -1] <- 0  # Set negative entries to 0

  negative_adj <- adj
  negative_adj[negative_adj == 1] <- 0   # Set positive entries to 0
  negative_adj[negative_adj == -1] <- 1  # Convert negative entries to 1

  # Create a N x K matrix whose (i, k) element represents the probability that node i belongs to block k.
   tau <-
     matrix(c(
       0.2, 0.5, 0.3,
       0.4, 0.4, 0.2,
       0.1, 0.4, 0.5,
       0.4, 0.4, 0.2,
       0.1, 0.1, 0.8,
       0.05, 0.05, 0.9
     ),
     nrow = K, ncol = N
     )
  tau <- t(tau)

  # Create a K x K matrix whose (k, l) element represents Pr(D_ij = 1 | Z_i = k, Z_j = l).
  sumTaus <- compute_sumTaus(N, K, tau)
  # Compute pi for positive entries
  pi_positive <- (t(tau) %*% positive_adj %*% tau) / sumTaus

  # Compute pi for negative entries
  pi_negative <- (t(tau) %*% negative_adj %*% tau) / sumTaus

  # Compute gamma (parameter of multinomial distribution)
  alpha <- colSums(tau)

  # Compute the true quadratic term in a naive way
  A <- matrix(0, nrow = N, ncol = K)

  for (i in 1:N) {
    for (k in 1:K) {
      for (j in 1:N) {
        if (i != j) {
          for (l in 1:K) {
            if (adj[i, j] == 0) {
              pi_ij <- 1 - pi_positive - pi_negative
            } else if (adj[i, j] == 1) {
              pi_ij <- pi_positive
            } else {
              pi_ij <- pi_negative
            }
            a_ij <- tau[j, l] * log(pi_ij[k, l])
            A[i, k] <- A[i, k] + a_ij
          }
        }
      }
    }
  }

  A <- 1 - A / 2

  # Divide A by alpha_{ik}
  A <- A / tau

  adj <- as(adj, "dgCMatrix")
  A_cpp <- compute_quadratic_term_signed(N, K, alpha, tau, adj, LB = 0)

  # Check if computation works as expected
  expect_equal(A, A_cpp, check.attributes = FALSE, tolerance = 1e-10)
})
