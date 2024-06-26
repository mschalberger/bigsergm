% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_hergm.R
\name{simulate_hergm}
\alias{simulate_hergm}
\title{Simulate a network}
\usage{
simulate_hergm(
  formula_for_simulation,
  data_for_simulation,
  colname_vertex_id,
  colname_block_membership,
  seed_edgelist = NULL,
  coef_within_block,
  coef_between_block,
  ergm_control = ergm::control.simulate.formula(),
  seed = NULL,
  directed = FALSE,
  n_sim = 1,
  output = "network",
  prevent_duplicate = TRUE,
  use_fast_between_simulation = FALSE,
  list_feature_matrices = NULL,
  verbose = 0,
  ...
)
}
\arguments{
\item{formula_for_simulation}{formula for simulating a network}

\item{data_for_simulation}{a data frame that contains vertex id, block membership, and vertex features.}

\item{colname_vertex_id}{a column name in the data frame for the vertex id}

\item{colname_block_membership}{a column name in the data frame for the block membership}

\item{seed_edgelist}{an edgelist used for creating a seed network. It should have the "edgelist" class}

\item{coef_within_block}{a vector of within-block parameters. The order of the parameters should match that of the formula.}

\item{coef_between_block}{a vector of between-block parameters. The order of the parameters should match that of the formula without externality terms.}

\item{ergm_control}{auxiliary function as user interface for fine-tuning ERGM simulation}

\item{seed}{seed value (integer) for network simulation.}

\item{directed}{whether the simulated network is directed}

\item{n_sim}{number of networks generated}

\item{output}{Normally character, one of "network" (default), "stats", "edgelist", to determine the output format.}

\item{prevent_duplicate}{If `TRUE`, the coefficient on nodematch("block") is set to be a very large negative number in drawing between-block links, so that there will be (almost) no within-block links.}

\item{use_fast_between_simulation}{If `TRUE`, this function uses an effcient way to simulate a between-block network. If the network is very large, you should consider using this option.
Note that when you use this, the first element of `coef_between_block` must be the edges parameter.}

\item{list_feature_matrices}{a list of feature adjacency matrices. If `use_fast_between_simulation`, this must be given.}

\item{verbose}{If this is TRUE/1, the program will print out additional information about the progress of simulation.}

\item{...}{Additional arguments, to be passed to lower-level functions}
}
\value{
Simulated networks, the output form depends on the parameter \code{output}
(default is a list of networks).
}
\description{
Simulate a network
}
\examples{
data(toyNet)

# Specify the model that you would like to estimate.
model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle

# Prepare a data frame that contains nodal id and covariates.
nodes_data <-
  data.frame(
    node_id = network::network.vertex.names(toyNet),
    block = network::get.vertex.attribute(toyNet, "block"),
    x = network::get.vertex.attribute(toyNet, "x"),
    y = network::get.vertex.attribute(toyNet, "y")
  )
# The feature adjacency matrices
list_feature_matrices <- bigergm::get_list_sparse_feature_adjmat(toyNet, model_formula)

# Simulate network stats
sim_stats <- bigergm::simulate_hergm(
  formula_for_simulation = model_formula,
  data_for_simulation = nodes_data,
  # Nodal data
  colname_vertex_id = "node_id",
  # Name of the column containing node IDs
  colname_block_membership = "block",
  # Name of the column containing block IDs
  coef_between_block = c(-4.5,0.8, 0.4),
  # The coefficients for the between connections
  coef_within_block = c(-1.7,0.5,0.6,0.15),
  n_sim = 10,
  # Number of simulations to return
  output = "stats",
  # Type of output
  list_feature_matrices = list_feature_matrices
  # Information on the covariates
)

}
