library(PLMIX)
library(PerMallows)

# 1. Simulation settings

alpha_list <- list(
  high = 0.45,
  medium = 0.15,
  low = 0.025
)

sample_sizes <- c(20, 50, 100)

B <- 30

central_ranking <- c(1, 2, 3, 4, 5)

# 2. Generate Mallows complete ranking data

generate_mallows_data <- function(alpha, N) {
  
  complete_rankings_num <- rmm(
    n = N,
    sigma0 = central_ranking,
    theta = alpha,
    dist.name = "kendall"
  )
  
  complete_rankings_num
}

# 3. Convert complete rankings into partial rankings

make_partial_data <- function(complete_rankings_num) {
  
  partial_output <- make_partial(
    data = complete_rankings_num,
    format_input = "ordering",
    probcens = c(0.25, 0.25, 0.25, 0.25)
  )
  partial_output$partialdata
}

# 4. Fit the Bayesian Plackett-Luce model

fit_bayesian_pl <- function(partial_rankings) {
  
  # Suppress iteration information printed by gibbsPLMIX
  capture.output(bayesian_pl_fit <- gibbsPLMIX(
    pi_inv = partial_rankings,
    K = 5,
    G = 1,
    n_iter = 1000,
    n_burn = 200,
    hyper = list(
      shape0 = matrix(
        1,
        nrow = 1,
        ncol = 5
      ),
      rate0 = 0.001,
      alpha0 = 1
    )
  ))
  
  # Posterior draws of the score parameters
  posterior_scores <- bayesian_pl_fit$P
  
  # Posterior means
  bayesian_scores <- colMeans(posterior_scores)
  names(bayesian_scores) <- paste0("A",1:5)
  
  # Estimated ranking
  bayesian_ranking <- names(sort(bayesian_scores, decreasing = TRUE))
  bayesian_ranking
}

# 5. Kendall distance

kendall_distance <- function(estimated_order, true_order) {
  
  distance <- 0
  n <- length(true_order)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      # Find the positions of the two items in the estimated ranking
      pos_i <- match(true_order[i], estimated_order)
      pos_j <- match(true_order[j], estimated_order)
      
      # Add one if their relative order is reversed
      if (pos_i > pos_j) {
        distance <- distance + 1
      }
    }
  }
  distance
}

# 6. Simulation main loop

set.seed(123)

true_order <- c("A1", "A2", "A3", "A4", "A5")

for (name in names(alpha_list)) {
  
  alpha <- alpha_list[[name]]
  
  for (N in sample_sizes) {
    
    # Initialise logical vector of length B with FALSE values
    bayesian_recovery <- logical(B)
    
    # Initialise numeric vector of length B with zero values
    bayesian_distance <- numeric(B)
    
    for (b in 1:B) {
      
      # Generate complete ranking data
      complete_rankings <- generate_mallows_data(alpha, N)
      
      # Convert complete rankings into partial rankings
      partial_rankings <- make_partial_data(complete_rankings)
      
      # Fit the Bayesian Plackett-Luce model  
      bayesian_order <- fit_bayesian_pl(partial_rankings)
      
      # Calculate Kendall distance
      bayesian_distance[b] <- kendall_distance(bayesian_order, true_order)
      
      # TRUE if the estimated ranking exactly matches the true ranking
      bayesian_recovery[b] <- bayesian_distance[b] == 0
    }
    
    # Print results for this setting
    cat("\n")
    cat("Separation:", name, "\n")
    cat("Sample size:", N, "\n")
    cat("Bayesian PL recovery rate:", mean(bayesian_recovery), "\n")
    cat("Bayesian PL mean Kendall distance:", mean(bayesian_distance), "\n")
  }
}