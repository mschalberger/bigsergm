test_that("signed probabilities between blocks", {
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
  tau <- matrix(0, nrow = N, ncol = K)
  for (i in 1:N) {
    tau[i, sample(1:K, 1)] <- 1
  }

  true_membership <- apply(tau, 1, which.max)

  # Create a K x K matrix whose (k, l) element represents Pr(D_ij = 1 | Z_i = k, Z_j = l)
  sumTaus <- compute_sumTaus(N, K, tau)

  # Compute pi for positive entries
  pi_positive <- (t(tau) %*% positive_adj %*% tau) / sumTaus

  # Compute pi for negative entries
  pi_negative <- (t(tau) %*% negative_adj %*% tau) / sumTaus

  # Compute pi for zero entries
  pi_zero <- 1 - pi_positive - pi_negative

  # Initialize matrices to store relative frequencies
  frequency_1 <- matrix(0, nrow = k, ncol = k)
  frequency_0 <- matrix(0, nrow = k, ncol = k)
  frequency_minus_1 <- matrix(0, nrow = k, ncol = k)

  # Loop through each block pair
  for (i in 1:K) {
    for (j in 1:K) {
      # Get nodes in block i and block j
      tmp <- adj[which(true_membership == i), which(true_membership == j)]

      # Get upper triangle of the adjacency matrix if i = j
      if (i == j) {
        tmp <- tmp[upper.tri(tmp)]
      }

      # Calculate the number of edges between nodes in block i and block j
      edges_1 <- sum(tmp == 1)
      edges_0 <- sum(tmp == 0)
      edges_minus_1 <- sum(tmp == -1)

      # Calculate the total number of possible edges (excluding diagonal and lower triangle)
      total_possible_edges <- length(tmp)

      # Calculate relative frequencies
      frequency_1[i, j] <- edges_1 / total_possible_edges
      frequency_0[i, j] <- edges_0 / total_possible_edges
      frequency_minus_1[i, j] <- edges_minus_1 / total_possible_edges
    }
  }

  expect_equal(pi_positive,frequency_1, check.attributes = FALSE, tolerance = 1e-10)
  expect_equal(pi_negative,frequency_minus_1, check.attributes = FALSE, tolerance = 1e-10)
  expect_equal(pi_zero,frequency_0, check.attributes = FALSE, tolerance = 1e-10)
})
