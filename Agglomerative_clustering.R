# By Dr. Rajen Shah (http://www.statslab.cam.ac.uk/~rds37/)
# Sent to Sahir upon request via email August 12, 2015
# Returns a sequence of groupings and corresponding maximal canonical correlations, and also
# a possibly "sign-corrected" version of the the design matrix x.
# Note the group numbers in the groupings are arbitrary and should be thought of as factors.

clustering <- function(x, switch_signs = TRUE, normalise = TRUE, verbose = FALSE, max_iter = ncol(x)) {
  n <- nrow(x); p <- ncol(x); min_np <- min(n, p)
  
  if (normalise) x <- scale(x)
  # Sigma will contain all canonical correlations
  Sigma <- abs(cor(x))
  diag(Sigma) <- -1  
  
  iter <- 1
  cor_max <- max(Sigma)
  which_cor_max <- which(Sigma == cor_max, arr.ind = TRUE)
  # Records pair with highest correlation
  group_pair <- as.integer(which_cor_max[1, ])
  max_group_size <- 1
  cur_groups <- 1:p
  distinct_groups <- 1:p
  # Records group membership at each iteration
  groups <- matrix(cur_groups, max_iter, p, byrow = TRUE)
  # Records rho hat max as a function of iterations
  # in order to select optimal clustering later
  cors <- c(cor_max, rep(0, max_iter-1))
  
  repeat {
    if (verbose) {
      print(iter)
    }
    iter <- iter + 1
    # Identifies which group will be grouped
    group1 <- which(cur_groups == group_pair[1])
    group2 <- which(cur_groups == group_pair[2])
    new_group_indices <- c(group1, group2)
    # Redefine indices
    new_group_num <- min(group_pair[1], group_pair[2])
    removed_group_num <- max(group_pair[1], group_pair[2])
    
    # Update signs - to understand what it does, see Dr. Shah's website
    if (switch_signs) {
      # Note: redefining rowMeans is actually faster than using the
      # original function on a non-coerced matrix
      group1_mean <- rowMeans2(x[, group1])
      group2_mean <- rowMeans2(x[, group2])
      x[, group2] <- x[, group2] * sign(sum(group1_mean * group2_mean))
    }
    # Update cur_groups, max_group_size, distinct_groups
    cur_groups[new_group_indices] <- new_group_num
    if (verbose) old_max_group_size <- max_group_size
    max_group_size <- max(length(new_group_indices), max_group_size)
    if (verbose && (max_group_size > old_max_group_size)) {
      print("New max_group_size")
      print(max_group_size)
    }
    distinct_groups <- distinct_groups[distinct_groups != removed_group_num]
    
    # Update Sigma
    new_group_x <- x[, new_group_indices]
    other_groups <- distinct_groups[distinct_groups != new_group_num]
    for (other_group in other_groups) {
      other_group_indices <- cur_groups == other_group
      cur_cor <- cancor(x[, other_group_indices], new_group_x)$cor[1]
      Sigma[other_group_indices, new_group_indices] <- cur_cor
      Sigma[new_group_indices, other_group_indices] <- cur_cor
    }
    Sigma[new_group_indices, new_group_indices] <- -1
    
    # Update group_pair
    cor_max <- max(Sigma)
    which_cor_max <- which(Sigma == cor_max, arr.ind = TRUE)
    group_pair <- c(cur_groups[which_cor_max[1, 1]], cur_groups[which_cor_max[1, 2]])
    
    # Update groups, cors, iter
    if ((max_group_size >= min_np) || iter >= max_iter) break
    groups[iter, ] <- cur_groups
    cors[iter] <- cor_max
  }
  return(list("groupings" = groups[seq_len(iter-1), ],
              "correlations" = cors[seq_len(iter-1)],
              "x" = x))
}

#########################################

rowMeans2 <- function(x, ...) {
  if (is.vector(x)) {
    return(x)
  } else {
    return(rowMeans(x, ...))
  }
}
