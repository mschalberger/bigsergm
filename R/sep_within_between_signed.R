sep_within_between <- function(network, block_membership){

  # Check if the network is signed
  if(!any(class(network) %in% c("static.sign","dynamic.sign"))){
    stop("The formula does not contain a signed network.")
  }
  # Get vertex ids
  net_ids <- network.vertex.names(network)
  block <- block_membership[match(net_ids, block_membership[,1]),2]
  network%v%"block" <- block

  # Get edgelist
  edgelist <- as.edgelist(network, attrname = "sign")
  colnames(edgelist) <- c("sender", "receiver", "sign")

  # Merge block membership into the edge list
  edgelist <- cbind(edgelist,
                    sender_block = block[edgelist[, 1]],
                    receiver_block = block[edgelist[, 2]])

  # Divide edgelist into within and between block edges
  within_edges <- edgelist[edgelist[, 4] == edgelist[, 5], 1:3]
  between_edges <- edgelist[edgelist[, 4] != edgelist[, 5], 1:3]

  # Turn the edgelist into a adjacecncy matrix
  n <- length(net_ids)
  mat <- matrix(NA, nrow = n, ncol = n)
  colnames(mat) <- rownames(mat) <- net_ids

  mat[within_edges[, 1:2]] <- "within"
  mat[between_edges[, 1:2]]  <- "between"

  # Add block membership to the network object
  network%e%"block" <- mat

  # Create within and between block networks
  within_network <- network::delete.edges(network, which(network%e%"block" == "between"))
  between_network <- network::delete.edges(network, which(network%e%"block" == "within"))

  return(list(within_network = within_network, between_network = between_network))
}

estimate_between <- function(formula) {
  # Create a formula that contains only dyad-independent terms. i.e. exclude externality terms like triangle.
  terms <- ergm::ergm_model(formula)$terms

  varnames <-
    list_rhs.formula(formula) %>%
    as.character()
  dep_terms <-
    terms %>% purrr::map(function(t) {
      dep <- t$dependence
      is_dep <- is.null(dep) || dep
    }) %>% unlist()
  between_rhs <- varnames[!dep_terms]
  between_rhs <- between_rhs[!is.na(between_rhs)]
  between_formula <- as.formula(glue::glue("g_logit ~ {paste(between_rhs, collapse = '+')}"))

  # Estimate logit
  between_logit <- mple_sign(
    formula = between_formula
  )

  # Remove unnecessary network objects
  between_logit$newnetwork <- NULL

  # Return the output
  return(between_logit)
}
