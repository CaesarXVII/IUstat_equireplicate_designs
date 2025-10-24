###             Functions used in                     ### 
### Incomplete U-Statistics of Equireplicate Designs: ###
### Berry-Esseen Bound and Efficient Construction     ###

## Remark 6, Section 4 ##

equirep_even <- function(n,r) {
  
  if (n %% 2 != 0) {
    stop("n must be even.")
  }
  
  if (r > (n-1)) {
    stop("r cannot be strictly greater than (n-1)")
  }
  
  # "g" and "l" are the "indices"
  # "g" indices the group, and "l" indices the elements of the group
  gs = rep(1:r,each=(n/2))
  ls = rep(0:(n/2-1),(r))
  # modular arithmetic for 1st column
  m1 = (gs+ls)%%(n-1)
  # replace zeros with n-1
  m1[which(m1==0)] = n-1
  # modular arithmetic for 2nd column
  m2 = (gs-ls)%%(n-1)
  # replace zeros with n-1
  m2[which(m2==0)]=n-1
  # the fixed points (i,i) need to be replaced with (i,n)
  m2[which(m1==m2)]=n
  
  # final design matrix, which is r-equireplicate
  design_mat = cbind(m1,m2)
  
  design_list = split(design_mat, row(design_mat))
  
  output <- list(design_list,design_mat)
  
  return(output)
  
}

equirep_odd <- function(n,r) {
  
  if (n %% 2 == 0) {
    stop("n must be odd.")
  }
  
  if (r > (n-1)) {
    stop("r cannot be strictly greater than (n-1).")
  }
  
  if (r %% 2 != 0) {
    stop("r must be even.")
  }
  
  # first column of the design matrix
  m1 = rep(1:n,r/2)
  # initial modular arithmetic, but may have zeros
  m2 = (m1+rep(1:(r/2),each=n))%%n
  # replace zeros with n
  m2[which(m2==0)]=n
  
  # final design matrix, which is r-equireplicate
  design_mat = cbind(m1,m2)
  
  design_list = split(design_mat, row(design_mat))
  
  output <- list(design_list,design_mat)
  
  return(output)
  
}

## Simulation Study 1 ##

exact_medheur <- function(d, tau = 1) {
  # d: Number of dimensions
  # tau: sd of the normal random variable
  
  # Compute sigma
  sigma <- tau* sqrt(qchisq(0.5, df = d))
  
  return(sigma)
}

median_heur_MMD <- function(X, Y) {
  comb_data <- rbind(X, Y)
  L2_dist <- dist(comb_data, method = "euclidean")
  sigma <- median(L2_dist)/sqrt(2)
  return(sigma)
}

# Linear kernel (works for pairs of points)
linear_kernelpp <- function(a, b, params) {
  return(a%*%b)  # Dot product between two vectors
}

# Gaussian RBF kernel (works for pairs of points)
rbf_kernelpp <- function(A, B, params) {
  sigma <- params$sigma
  dist_squared <- sum((A - B)^2)  # Squared Euclidean distance between two vectors
  return(exp(-dist_squared / (2 * sigma^2)))
}

# Polynomial kernel (works for pairs of points)
polynomial_kernelpp <- function(a, b, params) {
  degree <- params$degree
  return((a%*%b)^degree)  # Dot product raised to the power of 'degree'
}

mmd_iustat <- function(X, Y, D, kernel_function, kernel_params = NULL) {
  
  # Check if X and Y are vectors, if so convert them to matrices (column vectors)
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  # Get the number of samples
  m <- nrow(X)
  n <- nrow(Y)
  
  if (n != m) {
    warning("Same PMMat will be used for both samples, even if this option is suboptimal.")
  }
  
  iustat_vec <- apply(D,1, FUN = function(x) 
    
    kernel_function(X[x[1],], X[x[2],], kernel_params) +
      
      kernel_function(Y[x[1],], Y[x[2],], kernel_params) -
      
      kernel_function(X[x[1],], Y[x[2],], kernel_params) -
      
      kernel_function(X[x[2],], Y[x[1],], kernel_params)
    
  ) 
  
  return(iustat_vec)
}

## Simulation Study 2 ##

randes_full_k2 <- function(n, B) {
  idx <- sample(choose(n, 2), B)           # choose B columns
  t(combn(n, 2, simplify = TRUE)[, idx])   # each column is an unordered pair
}

## Simulation Study 3 ##

gcd <- function(a, b) {
  a <- abs(a);  b <- abs(b)
  while (b != 0L) {
    tmp <- b
    b   <- a %% b
    a   <- tmp
  }
  a
}

coprimes <- function(n) {
  if (n <= 1L) stop("n must be ≥ 2")
  Filter(function(g) gcd(g, n) == 1L, 1L:(n - 1L))
}

cyclic_design <- function(n = 37,k = 3, r=3, eta= 2^x) {
  if (k < 2L) {stop("k must be ≥ 2")}
  
  if (r > k*length(coprimes(n))) {stop("r must be smaller than k*phi(n)")}
  
  if (r %% k != 0) {
    stop("r must be divisible by k.")
  }
  
  G  <- coprimes(n)[1:(r/k)]                     # all admissible multipliers (injective case only)
  
  Js <- 0L:(k - 1L)                     # exponents {0,…,k-1}
  
  blocks <- list()
  blk_id <- 1L
  
  for (g in G) {
    for (i in 0L:(n - 1L)) {
      block <- (i + g * (eta(Js) - eta(0))) %% n
      blocks[[blk_id]] <- block
      blk_id <- blk_id + 1L
    }
  }
  attr(blocks, "n") <- n
  attr(blocks, "k") <- k
  attr(blocks, "multipliers") <- G
  
  k <- length(blocks[[1]])
  n_blocks <- length(blocks)
  
  # replace zeros with n
  blocks <- lapply(blocks, function(x) replace(x, x == 0, n))
  
  mat <- matrix(NA_integer_, nrow = n_blocks, ncol = k)
  
  for (i in seq_len(n_blocks)) {
    
    mat[i, ] <- blocks[[i]]
  }

  return(list(blocks, mat))
}

quad_mmd <- function(a1, a2, a3, a4, kernel_function, kernel_params = NULL) {
  
  kernel_function(a1, a3, kernel_params) +
    
    kernel_function(a2, a4, kernel_params) -
    
    kernel_function(a1, a4, kernel_params) -
    
    kernel_function(a2, a3, kernel_params)
  
}

hsic_iustat <- function(X, Y, D, kernel_function, kernel_paramsX = NULL,kernel_paramsY = NULL) {
  
  # Check if X and Y are vectors, if so convert them to matrices (column vectors)
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  # Get the number of samples
  m <- nrow(X)
  n <- nrow(Y)
  
  if (n != m) {
    warning("Same PMMat will be used for both samples, even if this option is suboptimal.")
  }
  
  iustat_vec <- apply(D,1, FUN = function(x) {
    
    d1 <- quad_mmd(X[x[1],],X[x[2],],X[x[3],],X[x[4],], kernel_function, kernel_params = kernel_paramsX)
    d1l <- quad_mmd(Y[x[1],],Y[x[2],],Y[x[3],],Y[x[4],], kernel_function, kernel_params = kernel_paramsY)
    d2 <- quad_mmd(X[x[1],],X[x[3],],X[x[2],],X[x[4],], kernel_function, kernel_params = kernel_paramsX) 
    d2l <- quad_mmd(Y[x[1],],Y[x[3],],Y[x[2],],Y[x[4],], kernel_function, kernel_params = kernel_paramsY)
    d3 <- quad_mmd(X[x[1],],X[x[4],],X[x[2],],X[x[3],], kernel_function, kernel_params = kernel_paramsX) 
    d3l <- quad_mmd(Y[x[1],],Y[x[4],],Y[x[2],],Y[x[3],], kernel_function, kernel_params = kernel_paramsY)
    
    return((d1*d1l + d2*d2l + d3*d3l)/12)
    
  }
  ) 
  
  return(iustat_vec)
}

## Simulation Study 4 ##

# Random design generation: no duplicates for sure but computationally intensive

randes_full_k4 <- function(n, B) {
  idx <- sample(choose(n, 4), B)           # choose B columns
  t(combn(n, 4, simplify = TRUE)[, idx])   # each column is an unordered quadruple
}

median_heur_HSIC <- function(Z) {
  Z  <- as.matrix(Z)
  d2 <- as.matrix(dist(Z))^2              # all pairwise squared distances
  m2 <- median(d2[d2 > 0])                # ignore exact zeros (ties)
  if (!is.finite(m2) || m2 == 0) m2 <- 1  # safe fallback
  sqrt(m2 / 2)
}

## Real Data Example ##

# Build two datasets X and Y from CIFAR-10 per-class matrices P$P_0..P_9
# - kx: integer vector of labels used to build X (and the (1-c) part of Y)
# - ky: integer vector of labels used to build the c part of Y
# - n : total rows in X and total rows in Y (each)
# - c : fraction of Y coming from ky (balanced across ky); (1-c) from kx (balanced across kx)
# - seed: optional RNG seed for reproducibility
# - shuffle: if TRUE, shuffles row order of X and Y at the end (recommended)
#
# Notes:
# * Balanced = as equal as possible per label. If n isn't divisible by |kx| or |ky|,
#   remainders are distributed (randomly) so counts differ by at most 1.
# * If a label appears in BOTH kx and ky, the function still works:
#   it draws disjoint indices for X, Y-from-kx, and Y-from-ky within that label.

# helper: split a total T as evenly as possible across 'labs' i.e., labels
split_counts <- function(T, labs) {
  g <- length(labs)
  if (g == 0L || T == 0L) return(setNames(integer(length(labs)), as.character(labs)))
  base <- T %/% g
  r    <- T %%  g
  cnt  <- rep(base, g)
  if (r > 0) {
    # randomly pick r labels to get +1; relies on set.seed() above if provided
    bump_idx <- sample.int(g, r, replace = FALSE)
    cnt[bump_idx] <- cnt[bump_idx] + 1L
  }
  setNames(cnt, as.character(labs))
}

build_XY_multi <- function(P, kx, ky, n, c, seed = NULL, shuffle = TRUE) {
  stopifnot(is.numeric(c), c >= 0, c <= 1)
  if (!is.null(seed)) set.seed(seed)
  
  kx <- as.integer(kx); ky <- as.integer(ky)
  if (length(kx) == 0L) stop("kx must contain at least one label.")
  if (length(ky) == 0L && c > 0) stop("ky must contain labels when c > 0.")
  
  # total rows for Y from kx / ky
  n1 <- floor((1 - c) * n)  # from kx into Y
  n2 <- n - n1              # from ky into Y
  
  # per-label counts
  cnt_X   <- split_counts(n,  kx)  # X only from kx
  cnt_Ykx <- split_counts(n1, kx)  # (1-c) part of Y from kx
  cnt_Yky <- split_counts(n2, ky)  # c     part of Y from ky
  
  # union of labels we'll touch
  all_labs <- sort(unique(c(kx, ky)))
  
  # storage
  X_blocks  <- list()
  Y_blocks  <- list()
  
  for (lab in all_labs) {
    key <- paste0("P_", lab)
    if (!key %in% names(P)) stop(sprintf("Missing matrix for label %d (expected P$%s).", lab, key))
    M <- P[[key]]
    n_avail <- nrow(M)
    
    x_need   <- if (lab %in% kx) unname(cnt_X[as.character(lab)])   else 0L
    ykx_need <- if (lab %in% kx) unname(cnt_Ykx[as.character(lab)]) else 0L
    yky_need <- if (lab %in% ky) unname(cnt_Yky[as.character(lab)]) else 0L
    
    total_need <- x_need + ykx_need + yky_need
    if (total_need > n_avail) {
      stop(sprintf("Not enough rows for label %d: need %d, have %d.", lab, total_need, n_avail))
    }
    
    if (total_need == 0L) next
    
    # draw disjoint indices for this label, then split into X, Y-from-kx, Y-from-ky
    idx <- sample.int(n_avail, size = total_need, replace = FALSE)
    take <- 0L
    if (x_need > 0L) {
      X_blocks[[length(X_blocks) + 1L]] <- M[idx[(take + 1L):(take + x_need)], , drop = FALSE]
      take <- take + x_need
    }
    if (ykx_need > 0L) {
      Y_blocks[[length(Y_blocks) + 1L]] <- M[idx[(take + 1L):(take + ykx_need)], , drop = FALSE]
      take <- take + ykx_need
    }
    if (yky_need > 0L) {
      Y_blocks[[length(Y_blocks) + 1L]] <- M[idx[(take + 1L):(take + yky_need)], , drop = FALSE]
      take <- take + yky_need
    }
  }
  
  # bind results
  X <- if (length(X_blocks)) do.call(rbind, X_blocks) else matrix(numeric(0), nrow = 0)
  Y <- if (length(Y_blocks)) do.call(rbind, Y_blocks) else matrix(numeric(0), nrow = 0)
  
  # sanity: exact sizes
  if (nrow(X) != n) stop(sprintf("Built X has %d rows; expected n = %d.", nrow(X), n))
  if (nrow(Y) != n) stop(sprintf("Built Y has %d rows; expected n = %d.", nrow(Y), n))
  
  # optional shuffle (order doesn't matter for MMD; this just avoids block structure)
  if (shuffle && nrow(X) > 0) X <- X[sample.int(nrow(X)), , drop = FALSE]
  if (shuffle && nrow(Y) > 0) Y <- Y[sample.int(nrow(Y)), , drop = FALSE]
  
  list(
    X = X,
    Y = Y,
    counts = list(
      X_per_label   = cnt_X,
      Y_from_kx     = cnt_Ykx,
      Y_from_ky     = cnt_Yky
    ),
    meta = list(kx = kx, ky = ky, n = n, c = c)
  )
}

#tst: two sample test

MMD_equirep_tst <- function(X,Y,r,kernel_function,kernel_params = NULL,alpha=0.05){
  
  # Check if X and Y are vectors, if so convert them to matrices (column vectors)
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  # Number of samples
  
  if (nrow(X) != nrow(Y)) {
    stop("Two samples must have the same number of rows.")
  }
  
  n <- nrow(X)
  
  D_size <- n*r/2
  
  D_equi <- equirep_even(n = n,r = r)[[2]]
  
  D_ind <- equirep_even(n = n,r = 1)[[2]]
  
  hS_vec_equi <- mmd_iustat(X = X,Y = Y,D = D_equi,kernel_function = kernel_function, kernel_params = kernel_params)
  
  sig2_ind <- var(mmd_iustat(X = X,Y = Y,D = D_ind,kernel_function = kernel_function, kernel_params = kernel_params)) 
  
  z_val <- (sqrt(D_size)*mean(hS_vec_equi))/sqrt(sig2_ind)
  
  p_val <- pnorm(q = z_val,mean = 0,sd = 1,lower.tail = FALSE) #one-side right tailed test
  
  z_test <- ifelse(test = p_val <= alpha,yes = 1,no = 0)
  
  return(list(z_val = z_val,
              p_val = p_val,
              z_test = z_test)
  )
  
}

MMD_rand_perm_rbf_tst <- function(X,Y,r,sig,n_perm=1000,alpha=0.05){
  
  # Check if X and Y are vectors, if so convert them to matrices (column vectors)
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  
  # Number of samples
  
  if (nrow(X) != nrow(Y)) {
    stop("Two samples must have the same number of rows.")
  }
  
  n <- nrow(X)
  
  D_size <- n*r/2
  
  #D_rand <- randes_full_k2(n = n,B = D_size) #if you want to use random designs instead
  D_rand <- equirep_even(n = n,r = r)[[2]]
  
  rand_mmd <- mean(mmd_iustat(X = X,Y = Y,D = D_rand,kernel_function = rbf_kernelpp, kernel_params = list(sigma = sig)))
  
  # Combine the two samples into one dataset
  comb_samples <- rbind(X, Y)
  tot_samples <- nrow(comb_samples)
  
  # Initialize a vector to store MMD values from permuted samples
  perm_mmds <- numeric(n_perm)
  
  # Perform the permutation test
  for (i in 1:n_perm) {
    
    # Permute the combined samples
    perm_ind <- sample(1:tot_samples, tot_samples, replace = FALSE)
    X_perm <- comb_samples[perm_ind[1:n], , drop = FALSE]
    Y_perm <- comb_samples[perm_ind[(n+1):tot_samples], , drop = FALSE]
    
    #sig <- median_heur_MMD(X = X_perm,Y = Y_perm)
    
    # Calculate MMD for the permuted samples
    perm_mmds[i] <- mean(mmd_iustat(X = X_perm,Y = Y_perm,D = D_rand,kernel_function = rbf_kernelpp, kernel_params = list(sigma = sig)))
    
  }
  
  # Calculate the p-value: the proportion of permuted MMDs that are greater than or equal to the original MMD
  p_val <- (sum(perm_mmds >= rand_mmd) +1) /(n_perm + 1)
  
  z_perm_test <- ifelse(test = p_val <= alpha,yes = 1,no = 0)
  
  # Return the p-value, the original MMD value, and the null distribution
  return(list(p_val = p_val, 
              z_perm_test = z_perm_test) 
  )
  
}
