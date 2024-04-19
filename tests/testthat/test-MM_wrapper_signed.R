test_that("Test if estimated block membership is correct", {
  # Simulate a network with local dependence
  k <- 4
  n <- 100
  nodes_data <- data.frame(node_id = 1:n, block = sample(1:k, size=n, replace = T))
  sim_dep <- simulate_hergm(formula_for_simulation = net ~ Pos(~edges) + Neg(~edges),
                            data_for_simulation = nodes_data,
                            colname_vertex_id = "node_id",
                            colname_block = "block",
                            coef_within_block =  c(3,0),
                            coef_between_block = c(-3,-3),
                            use_fast_between_simulation = F)
  sim_dep_sgl <- UnLayer(sim_dep)

  # Estimate block membership
  mm <- MM_wrapper_signed(network = sim_dep, formula = sim_dep ~ Pos(~edges) + Neg(~edges),
                          n_clusters = k, n_MM_step_max = 100, cache = NULL, tol_MM_step = 0.0001, compute_pi = T,
                          check_alpha_update = F)

  # Get indices of unique values
  indices_mm <- lapply(unique(mm$z_memb_final), function(x) which(mm$z_memb_final == x))
  indices_block <- lapply(unique(sim_dep_sgl%v%"block"), function(x) which(sim_dep_sgl%v%"block" == x))

  # Check if blockmemberships are identical
  expect_identical(indices_mm, indices_block)
  # Check if lower bound is increasing
  expect_true(all(diff(mm$lower_bound) > 0))
})
