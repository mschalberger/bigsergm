## ----include = FALSE----------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bigergm)
library(ergm)
library(dplyr)


## ----message=FALSE------------------------------------------------------------
# Load an embedded network object.
data(toyNet)
# Draw the network.
plot(toyNet, vertex.col = rep(c("tomato", "steelblue", "darkgreen", "black"),
                        each = toyNet$gal$n/4))

## ----message=FALSE------------------------------------------------------------
model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle

hergm_res <-
  bigergm::hergm(
    # The model you would like to estiamte
    object = model_formula,
    # The number of blocks
    n_clusters = 4,
    # The maximum number of MM algorithm steps
    n_MM_step_max = 100,
    # Perform parameter estimation after the block recovery step
    estimate_parameters = TRUE,
    # Indicate that clustering must take into account nodematch on characteristics
    clustering_with_features = TRUE,
    # Keep track of block memberships at each EM iteration
    check_block_membership = TRUE
)


## -----------------------------------------------------------------------------
plot(1:length(hergm_res$MM_lower_bound),
     hergm_res$MM_lower_bound, type = "l", xlab = "Iterations", ylab = "Lower Bound")

## -----------------------------------------------------------------------------
# Number of nodes per recovered block:
table(hergm_res$partition)

## -----------------------------------------------------------------------------
# For the between networks
summary(hergm_res$est_between)

## -----------------------------------------------------------------------------
# For the within networks
summary(hergm_res$est_within)

## ----message=FALSE------------------------------------------------------------
# Prepare a data frame that contains nodal id and covariates.
nodes_data <-
  tibble::tibble(
    node_id = network::network.vertex.names(toyNet),
    block = hergm_res$partition,
    x = network::get.vertex.attribute(toyNet, "x"),
    y = network::get.vertex.attribute(toyNet, "y")
    )

# The feature adjacency matrices
list_feature_matrices <- bigergm::get_list_sparse_feature_adjmat(toyNet, model_formula)

# The MCMC settings
sim_ergm_control <- ergm::control.simulate.formula(
  MCMC.burnin = 1000000,
  MCMC.interval = 100000
)

# The feature adjacency matrices
list_feature_matrices <- bigergm::get_list_sparse_feature_adjmat(toyNet, model_formula)

gof_res <- bigergm::gof_bigergm(
  toyNet,
  # The feature adjacency matrices
  list_feature_matrices = list_feature_matrices,
  # A dataframe containing the nodes data.
  data_for_simulation = nodes_data,
  # The name of the nodes_data column containing the node IDs
  # which are used within the network g
  colname_vertex_id = 'node_id',
  # The name of the nodes_data column containing the block ID.
  colname_block_membership = 'block',
  # The object returned by bigergm::hergm()
  bigergm_results = hergm_res,
  # The MCMC settings
  ergm_control = sim_ergm_control,
  # The number of simulations to use
  n_sim = 100
)

## ----message=FALSE, warning=FALSE---------------------------------------------
degree_gof <- 
  gof_res$simulated$degree_dist %>%
  dplyr::group_by(degree) %>%
  dplyr::summarise(log_mean_share = mean(log(share)),
                   log_sd_share = sd(log(share))) %>%
  dplyr::ungroup()
plot(degree_gof$degree, degree_gof$log_mean_share,
     xlab = "Degree", ylab = "Log Prop. of Nodes",
     ylim = c(-5.5,-1.8), xlim = c(6,20), type = "l", lty = 2)
lines(degree_gof$degree, degree_gof$log_mean_share+ 1.96 * degree_gof$log_sd_share, type = "l", lty = 2)
lines(degree_gof$degree, degree_gof$log_mean_share- 1.96 * degree_gof$log_sd_share, type = "l", lty = 2)
tmp_info <- gof_res$original$degree_dist %>% 
  dplyr::filter(share > 0 & degree < 22)
lines(tmp_info$degree, log(tmp_info$share), lty = 1)

## ----message=FALSE, echo=TRUE-------------------------------------------------
# Estimated coefficients for the between-community connections
coef_between_block <- coef(hergm_res$est_between)

# Estimated coefficients for the within-community connections
coef_within_block <- coef(hergm_res$est_within)

sim_net <- bigergm::simulate_hergm(
  # Formula for between-blocks
  formula_for_simulation = model_formula,
  # Same as for gof, a dataframe containing nodes attributes
  data_for_simulation = nodes_data,
  # Name of the column containing node IDs
  colname_vertex_id = "node_id",
  # Name of the column containing block IDs
  colname_block_membership = "block",
  # The coefficients for the between connections
  coef_between_block = coef_between_block,
   # The coefficients for the within connections
  coef_within_block = coef_within_block,
  # The MCMC settings
  ergm_control = sim_ergm_control,
  # Number of simulations to return
  n_sim = 1,
  # If `stats` a list with network statistics 
  # for the between and within connections is returned
  output = "network",
  # Simulates between connections by drawing from a logistic distribution. 
  # If FALSE, draws between connections by MCMC.
  use_fast_between_simulation = TRUE,
  # The feature adjacency matrices
  list_feature_matrices = list_feature_matrices
)

## -----------------------------------------------------------------------------
plot(sim_net)

## ----message=FALSE------------------------------------------------------------
hergm_res_second <-
  bigergm::hergm(object = hergm_res)

