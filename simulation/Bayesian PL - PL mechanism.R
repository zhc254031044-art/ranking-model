library(PLMIX)

# 1. Simulation settings

# True score settings
score_list <- list(
  high = c(A1 = 5, A2 = 4, A3 = 3, A4 = 2, A5 = 1),
  medium = c(A1 = 5, A2 = 4.5, A3 = 4, A4 = 3.5, A5 = 3),
  low = c(A1 = 5, A2 = 4.9, A3 = 4.8, A4 = 4.7, A5 = 4.6)
)

sample_sizes <- c(20, 50, 100)

B <- 30

# 2. Generate Plackett-Luce complete ranking data

generate_pl_data <- function(w_true, N) {
  
  items <- c(1, 2, 3, 4, 5)
  
  # Create an empty matrix to store complete rankings
  complete_rankings <- matrix(NA, nrow = N, ncol = length(items))
  
  # Generate each complete ranking
  for (n in 1:N) {
    
    remaining_items <- items
    
    for (k in 1:length(items)) {
      
      # Calculate selection probabilities for the remaining items
      probabilities <- w_true[remaining_items] / sum(w_true[remaining_items])
      
      # Select the position of one remaining item according to the PL probabilities
      selected_position  <- sample(1:length(remaining_items), size = 1, prob = probabilities)
      
      # Convert the selected position back to the actual item
      selected_item <- remaining_items[selected_position]
      
      # Store the selected item at the current rank
      complete_rankings[n, k] <- selected_item
      
      # Remove the selected item from the remaining set
      remaining_items <- setdiff(remaining_items, selected_item)
    }
  }
  complete_rankings
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

for (score_name in names(score_list)) {
  
  w_true <- score_list[[score_name]]
  
  for (N in sample_sizes) {
    
    # Initialise logical vector of length B with FALSE values
    bayesian_recovery <- logical(B)
    
    # Initialise numeric vector of length B with zero values
    bayesian_distance <- numeric(B)
    
    for (b in 1:B) {
      
      # Generate complete ranking data
      complete_rankings <- generate_pl_data(w_true, N)
      
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
    cat("Separation:", score_name, "\n")
    cat("Sample size:", N, "\n")
    cat("Bayesian PL recovery rate:", mean(bayesian_recovery), "\n")
    cat("Bayesian PL mean Kendall distance:", mean(bayesian_distance), "\n")
  }
}

# 7. Compare posterior means using 1000 and 2000 Gibbs iterations

set.seed(123)

w_true <- score_list$medium
N <- 50

# Generate one simulated dataset
complete_rankings <- generate_pl_data(w_true, N)
partial_rankings <- make_partial_data(complete_rankings)

# Fit with 1000 iterations
fit_1000  <- gibbsPLMIX(
  pi_inv = partial_rankings,
  K = 5,
  G = 1,
  n_iter = 1000,
  n_burn = 200,
  hyper = list(
    shape0 = matrix(1, nrow = 1, ncol = 5),
    rate0 = 0.001,
    alpha0 = 1
  )
)

# Fit with 2000 iterations
fit_2000  <- gibbsPLMIX(
  pi_inv = partial_rankings,
  K = 5,
  G = 1,
  n_iter = 2000,
  n_burn = 400,
  hyper = list(
    shape0 = matrix(1, nrow = 1, ncol = 5),
    rate0 = 0.001,
    alpha0 = 1
  )
)

# Posterior draws of the score parameters
posterior_scores_1000 <- fit_1000$P
posterior_scores_2000 <- fit_2000$P

# Posterior means
bayesian_scores_1000 <- colMeans(posterior_scores_1000)
names(bayesian_scores_1000) <- paste0("A",1:5)

bayesian_scores_2000 <- colMeans(posterior_scores_2000)
names(bayesian_scores_2000) <- paste0("A",1:5)

bayesian_scores_1000
bayesian_scores_2000
