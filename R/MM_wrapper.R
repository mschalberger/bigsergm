# Function that checks if the algorithm converged or not
check_has_converged <- function(old_value, new_value, tolerance) {
  realized_tol <- abs(new_value - old_value) / abs(old_value)
  if (realized_tol < tolerance) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Wrapper function that conducts clustering.
MM_wrapper_signed <-
  function(network,
           formula,
           n_clusters,
           n_MM_step_max,
           tol_MM_step,
           min_size = 2,
           initialization_method = 1,
           use_infomap_python = FALSE,
           virtualenv_python = "r-bigergm",
           seed_infomap = NULL,
           initialized_cluster_data = NULL,
           clustering_with_features = FALSE,
           list_multiplied_feature_matrices = NULL,
           verbose = 0,
           weight_for_initialization = 1000,
           compute_pi = FALSE,
           check_alpha_update = FALSE,
           check_block_membership = FALSE,
           MM_restart_object = NULL,
           minAlpha = 1e-6,
           cache,
           ...) {
    network <- UnLayer(network)
    n_nodes <- as.integer(network$gal$n)
    is_undirected <- !network$gal$directed

    # Cache eigenvectors computation in disk to free memory for other computations
    if (is.null(cache)) {
      eigenvectors_sparse_fn <- eigenvectors_sparse
    } else {
      eigenvectors_sparse_fn <- memoise::memoise(eigenvectors_sparse, cache = cache)
    }

    # If restart_object is NULL, initialize memberships.
    if (is.null(MM_restart_object)) {
      # Step 0: Initialization
      n_nodes <- as.integer(network$gal$n)
      is_undirected <- !network$gal$directed

      if (verbose > 0) message("Converting network to edgelist...")
      edgelist <- network::as.edgelist(network, attrname = "sign")

      if (verbose > 0) message("Converting edgelist to sparse matrix...")
      g <- Matrix::sparseMatrix(i = edgelist[, 1], j = edgelist[, 2], x = edgelist[,3], dims = c(n_nodes, n_nodes), symmetric = is_undirected)
      g_pos <- g
      g_pos[g_pos < 0] <- 0
      g_neg <- g
      g_neg[g_neg > 0] <- 0

      ## Calculate network statistics needed for variational approximation of p(Z|X)
      if (verbose > 0) message("\nStep 1: Initialize z")

      # Step 1: Get initial estimate of Z memberships
      # If you already have an initialized clustering result, start the MM iterations from the given clustering result.
      if (!is.null(initialized_cluster_data)) {
        z_memb_init <- factor(initialized_cluster_data, levels = 1:n_clusters)
        z_memb <- check_clusters_sparse(z_memb_init, g, n_clusters, eigenvectors_sparse_fn, min_size, verbose = verbose)
      }

      # Start cluster initialization
      else if (initialization_method == 1) # Infomap
      {
        if (verbose > 0) message("Using Infomap to initialize clustering step...")

        if(use_infomap_python){

          # first check if an env is available
          if(!reticulate::py_available())
          {
            message("No Python Environemt available. Use check_and_install() ",
                    "to install recommended environment.")
            invisible(return(NULL))
          }
          reticulate::use_virtualenv(virtualenv_python, required = TRUE)
          have_infomap <- reticulate::py_module_available("infomap")
          have_numpy <- reticulate::py_module_available("numpy")

          if (have_numpy + have_infomap != 2) {
            use_infomap_python <- F
            message("Since either Numpy or Infomap is not available in Pythin the igraph implementation is used.")
          }
        }

        if (use_infomap_python) {
          python_path <- system.file("python", package = "bigergm")
          opt <- reticulate::import_from_path("infomap_R", path = python_path)

          z_memb_init <- factor(opt$infomap_python(matrix = edgelist, n_clusters = n_clusters), levels = 1:n_clusters)
          z_memb <- check_clusters_sparse(z_memb_init, g, n_clusters, eigenvectors_sparse_fn, min_size, verbose = verbose)
        } else {
          # Use igraph's infomap
          set.seed(seed_infomap)
          b <- igraph::cluster_infomap(intergraph::asIgraph(network))
          z_memb_init <- factor(b$membership, levels = 1:n_clusters)
          z_memb <- check_clusters_sparse(z_memb_init, g, n_clusters, eigenvectors_sparse_fn, min_size, verbose = verbose)
        }
      } else if (initialization_method == 2) {
        if (verbose > 0) message("Initializing with uniformly distributed clusters...")
        z_memb <- z_memb_init <- factor(sample.int(n_clusters, size = n_nodes, replace = TRUE))
      } else # Spectral clustering
      {
        z_memb <- spec_clust_sparse(g, n_clusters, eigenvectors_sparse_fn)
        z_memb_init <- z_memb
      }

      # Step 2: Find A(Z=z) ~ P(Z=z|X=x)
      if (verbose > 0) {
        message(paste("\n\nStep 2: Find variational approximation A(Z=z) ~ P(Z=z|X=x)", sep = ""))
      }
      block_sizes <- table(z_memb)
      memb_inds <- as.numeric(names(block_sizes))
      n_blocks <- length(block_sizes)

      # Initialize alpha
      if (verbose > 0) message("Initializing alpha...")
      alpha <- matrix(1, n_nodes, n_clusters)
      if (verbose > 0) message("Replacing zeros in alpha...")
      for (i in which(!is.na(z_memb))) {
        alpha[i, z_memb[i]] <- weight_for_initialization
      }
      if (verbose > 0) message("Normalizing alpha...")
      alpha <- alpha / Matrix::rowSums(alpha)

      # Keep track of alpha, z, changes in alpha, and the lower bound
      list_alpha <- list()
      list_z <- list()
      change_in_alpha <- c()
      lower_bound <- c()
      if (check_alpha_update) {
        list_alpha[[1]] <- alpha
      }
      if (check_block_membership) {
        list_z[[1]] <- z_memb
      }
      # Set the counter of MM iteration
      counter_e_step <- 0
    }
    # If restart_object is given, restore the values.
    else {
      if (verbose > 0) message("Making alpha matrix dense...")
      alpha <- as.matrix(MM_restart_object$alpha)
      alpha[alpha == 0] <- minAlpha
      list_alpha <- MM_restart_object$list_alpha
      z_memb_init <- MM_restart_object$z_memb_init
      list_z <- MM_restart_object$list_z
      change_in_alpha <- MM_restart_object$change_in_alpha
      lower_bound <- MM_restart_object$lower_bound
      counter_e_step <- MM_restart_object$counter_e_step
      g <- MM_restart_object$adjacency_matrix
      n_blocks <- n_clusters
    }

    # Get a sparse feature adjacency matrix if clustering with features
    if (clustering_with_features) {
      # If a list of feature matrices is given from outside, skip this step.
      if (!is.null(list_multiplied_feature_matrices)) {
        if (verbose > 0) {
          message("Skipped forming features matrix")
        }
      } else {
        if (verbose > 0) {
          message("Forming features matrix")
        }
        list_feature_matrices <- get_list_sparse_feature_adjmat(network, formula)
        if (verbose > 0) {
          message("Got features matrix")
        }
      }
    }

    # Get multiplied feature matrices if clustering with features
    if (clustering_with_features) {
      if (is.null(list_multiplied_feature_matrices)) {
        if (verbose > 0) {
          message("Getting multiplied matrices used for E step")
        }
        list_multiplied_feature_matrices <- get_elementwise_multiplied_matrices(g, list_feature_matrices)
        if (verbose > 0) {
          message("Got multiplied matrices")
        }
        rm(list_feature_matrices)
      } else {
        message("Skipped getting multiplied matrices")
      }
    }

    # Compute gamma (parameter for multinomial distribution)
    gamma <- colMeans(alpha)


    # If clustering with features:
    if (verbose > 0) {
      message(paste("MM algorithm started at ", counter_e_step))
      MM_start <- Sys.time()
      message(paste("Started at ", MM_start))
    }

    # This is for internal use
    inner_counter <- 0

    if (clustering_with_features) {
      # MM iterations start
      repeat {
        inner_counter <- inner_counter + 1
        counter_e_step <- counter_e_step + 1

        if (verbose > 0) {
          message(paste("MM algorithm iteration: ", counter_e_step))
          iter_start <- Sys.time()
          message(glue::glue("MM algorithm iteration {counter_e_step} started at {iter_start}"))
        }

        # If you would like to see the change in \alpha during MM iterations, keep \alpha_{t-1} before updating.
        if (check_alpha_update) {
          alpha_prev <- rlang::duplicate(alpha)
        }

        # Update alpha
        alpha_LB <-
          run_MM_with_features(
            as.integer(n_nodes),
            as.integer(n_clusters),
            as.double(gamma),
            list_multiplied_feature_matrices,
            alpha,
            verbose
          )

        alpha <- alpha_LB[[1]]
        gamma <- colMeans(alpha)

        # Keep track of the lower bound:
        lower_bound <- c(lower_bound, alpha_LB[[2]])
        # Check if the model is converged
        if (inner_counter > 1) {
          if (check_has_converged(old_value = lower_bound[inner_counter - 1],
                                  new_value = lower_bound[inner_counter],tolerance = tol_MM_step)) {
            break
          }
        }

        # If you would like to see the change in \alpha during MM iterations, print max_{i, k} |\alpha_{t}(i, k) - \alpha_{t-1}(i, k)|.
        if (check_alpha_update) {
          list_alpha[[counter_e_step + 1]] <- alpha
          change_in_alpha <- c(change_in_alpha, max(abs(alpha - alpha_prev)))
        }

        # If you would like to keep track on block memberships over MM iterations:
        if (check_block_membership) {
          z <- factor(apply(alpha, 1, which.max), levels = 1:n_clusters)
          list_z[[counter_e_step + 1]] <- z
        }

        if (verbose > 0) {
          iter_end <- Sys.time()
          message(glue::glue("MM algorithm iteration {counter_e_step} finished at {iter_end}"))
          message(glue::glue("MM iterations {counter_e_step} took {signif(difftime(iter_end, iter_start, units = 'mins'), digits = 3)} minutes."))
        }

        if (inner_counter %% 5) {
          gc()
        }

        if (inner_counter >= n_MM_step_max) {
          break
        }
      }
      # End of MM iteration
      if (verbose > 0) {
        MM_end <- Sys.time()
        message("Finished the MM cycle at ", MM_end)
        diff <- MM_end - MM_start
        message(glue::glue("{inner_counter} MM iterations took {signif(difftime(MM_end, MM_start, units = 'mins'), digits = 3)} minutes.\n"))
        message(glue::glue("Total MM iterations: {counter_e_step}"))
      }

      # Compute pi if necessary
      Pi <- NULL
      if (compute_pi) {
        Pi <- compute_pi_with_features(n_nodes, n_clusters, list_multiplied_feature_matrices, alpha)
      }
    } else { # If clustering without features:
      list_multiplied_feature_matrices <- NULL
      # MM iterations start
      repeat {
        counter_e_step <- counter_e_step + 1
        inner_counter <- inner_counter + 1

        if (verbose > 0) {
          message(paste("MM algorithm iteration: ", counter_e_step))
          iter_start <- Sys.time()
          message(glue::glue("MM algorithm iteration {counter_e_step} started at {iter_start}"))
        }

        # If you would like to see the change in \alpha during MM iterations, keep \alpha_{t-1} before updating.
        if (check_alpha_update) {
          alpha_prev <- rlang::duplicate(alpha)
        }

        alpha_LB <-
          run_MM_without_features_signed(
            n_nodes,
            n_blocks,
            gamma,
            alpha,
            g,
            verbose
          )

        # \alpha and \gamma
        alpha <- alpha_LB[[1]]
        gamma <- colMeans(alpha)

        # Keep track of the lower bound:
        lower_bound <- c(lower_bound, alpha_LB[[2]])
        # Check if the model is converged
        if (inner_counter > 1) {
          if (check_has_converged(old_value = lower_bound[inner_counter - 1],
                                  new_value = lower_bound[inner_counter],tolerance = tol_MM_step)) {
            break
          }
        }

        # If you would like to see the change in \alpha during MM iterations, print max_{i, k} |\alpha_{t}(i, k) - \alpha_{t-1}(i, k)|.
        if (check_alpha_update) {
          list_alpha[[counter_e_step + 1]] <- alpha
          change_in_alpha <- c(change_in_alpha, max(abs(alpha - alpha_prev)))
        }

        # If you would like to keep track on block memberships over MM iterations:
        if (check_block_membership) {
          z <- factor(apply(alpha, 1, which.max), levels = 1:n_clusters)
          list_z[[counter_e_step + 1]] <- z
        }

        if (verbose > 0) {
          iter_end <- Sys.time()
          message(glue::glue("MM algorithm iteration {counter_e_step} finished at {iter_end}"))
          message(glue::glue("MM iterations {counter_e_step} took {signif(difftime(iter_end, iter_start, units = 'mins'), digits = 3)} minutes."))
        }

        if (inner_counter %% 5) {
          gc()
        }

        if (inner_counter >= n_MM_step_max) {
          break
        }
      }

      # End of MM iteration
      if (verbose > 0) {
        MM_end <- Sys.time()
        message("Finished the MM cycle at ", MM_end)
        diff <- MM_end - MM_start
        message(glue::glue("{inner_counter} MM iterations took {signif(difftime(MM_end, MM_start, units = 'mins'), digits = 3)} minutes.\n"))
        message(glue::glue("Total MM iterations: {counter_e_step}"))
      }

      # Compute pi if necessary
      Pi_pos <- NULL
      Pi_neg <- NULL
      if (compute_pi) {
        Pi_pos <- (t(alpha) %*% g_pos %*% alpha) / compute_sumTaus(n_nodes, n_blocks, alpha)
        Pi_neg <- (t(alpha) %*% g_neg %*% alpha) / compute_sumTaus(n_nodes, n_blocks, alpha)
      }
    }

    # Estimate the block membership
    z_memb <-
      factor(apply(alpha, 1, which.max), levels = 1:n_clusters)
    # Keep the final clustering result before check_cluster_sparse()
    z_memb_final_before_kmeans <- z_memb

    if (verbose > 0) message("Making alpha matrix sparse...")
    alpha[alpha <= minAlpha] <- 0
    alpha <- as(alpha, "dgCMatrix")

    # Check bad clusters
    z_memb_final <-
      factor(check_clusters_sparse(z_memb, g, n_clusters, eigenvectors_sparse_fn, min_size, verbose = verbose),
             levels = 1:n_clusters
      )

    # Return the output
    return(
      list(
        Pi_pos = Pi_pos,
        Pi_neg = Pi_neg,
        z_memb_init = z_memb_init,
        list_z = list_z,
        z_memb_final_before_kmeans = z_memb_final_before_kmeans,
        z_memb_final = z_memb_final,
        list_alpha = list_alpha,
        change_in_alpha = change_in_alpha,
        lower_bound = lower_bound,
        counter_e_step = counter_e_step,
        adjacency_matrix = g,
        alpha = alpha,
        list_multiplied_feature_matrices = list_multiplied_feature_matrices
      )
    )
  }

# ------------------------------------------------------------------------
# -------------- Auxiliary functions -------------------------------------
# ------------------------------------------------------------------------

permute_tau <- function(tau, labels) {
  new_tau <- tau
  for (i in 1:length(labels)) {
    # temp[z_memb == i] <- labels[as.numeric(names(labels))==i]
    new_tau[, i] <- tau[, labels[i]]
  }
  new_tau
}


# Function for spectral clustering
# @param network a sparse adjacency matrix
# @param n_clusters number of specified clusters
# @param eigenvectors_fn a function that performs eigenvector decomposition
spec_clust_sparse <- function(network, n_clusters, eigenvectors_fn) {
  n <- nrow(network)
  n_vec <- ceiling(sqrt(n))
  message("Calculating eigenvectors...")
  b <- eigenvectors_fn(network, n_vec)
  message("K-means clustering...")
  c <- kmeans(
    b,
    centers = n_clusters,
    nstart = 100,
    iter.max = 20,
    algorithm = "Hartigan-Wong"
  )
  c.ind <- c$cluster
  z_memb <- factor(c.ind, levels = 1:n_clusters)
  z_memb
}


check_clusters_sparse <-
  function(z_memb, network, n_clusters, eigenvectors_fn, min_size = 2, verbose = verbose) {
    n <- nrow(network)
    n_vec <- ceiling(sqrt(n))

    if (verbose > 0) {
      message("Eigenvalue decomposition")
    }

    b <- eigenvectors_fn(network, n_vec)
    z_memb <- factor(z_memb, levels = 1:n_clusters)

    if (verbose > 0) {
      message("Checking for bad clusters...")
    }

    counter_bad_cluster <- 0
    length_prev_bad_clusters <- NULL
    repeat {
      z_memb_temp <- as.vector(z_memb)
      z_memb <- factor(z_memb_temp, levels = 1:n_clusters)
      bad_clusters <- which(table(z_memb) < min_size)

      if (verbose > 0) {
        message(paste(c("Remaining bad clusters:", length(bad_clusters))))
      }

      # If the number of bad clusters doesn't change over one iteration:
      if (!is.null(length_prev_bad_clusters)) {
        if (length(bad_clusters) == length_prev_bad_clusters) {
          counter_bad_cluster <- counter_bad_cluster + 1
        } else {
          counter_bad_cluster <- 0
        }
      }
      length_prev_bad_clusters <- length(bad_clusters)

      # If removing bad clusters seems to fail, stop the process.
      if (counter_bad_cluster >= 100) {
        stop("\nRemoving bad clusters failed. The number of pre-specified blocks might be too high.")
      }

      if (length(bad_clusters) == 0) {
        break
      }

      bad_cluster <- which(table(z_memb) < 2)[1]
      donor <- which.max(table(z_memb))
      indeces <- c(bad_clusters, donor)
      temp <-
        kmeans(
          b[z_memb %in% indeces, ],
          centers = 2,
          nstart = 100,
          iter.max = 20,
          algorithm = "Hartigan-Wong"
        )
      for (i in 1:length(indeces)) {
        z_memb_temp[z_memb %in% indeces][temp$cluster == i] <- indeces[i]
      }
      z_memb <- factor(z_memb_temp, levels = 1:n_clusters)
    }

    if (verbose > 0) {
      message("Done checking clusters")
    }

    return(z_memb)
  }
