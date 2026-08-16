library(PLMIX)

# =========================================================
# 1. Simulation settings
# =========================================================

score_list <- list(
  high = c(A1 = 5, A2 = 4, A3 = 3, A4 = 2, A5 = 1),
  medium = c(A1 = 5, A2 = 4.5, A3 = 4, A4 = 3.5, A5 = 3),
  low = c(A1 = 5, A2 = 4.9, A3 = 4.8, A4 = 4.7, A5 = 4.6)
)

true_order <- c("A1", "A2", "A3", "A4", "A5")


# =========================================================
# 2. Generate complete rankings under PL mechanism
# =========================================================

generate_pl_data <- function(w_true, N) {
  
  items <- names(w_true)
  
  complete_rankings <- matrix(
    NA,
    nrow = N,
    ncol = length(items),
    dimnames = list(
      NULL,
      paste0("rank_", 1:length(items))
    )
  )
  
  for (r in 1:N) {
    
    remaining_items <- items
    
    for (k in 1:length(items)) {
      
      probabilities <-
        w_true[remaining_items] /
        sum(w_true[remaining_items])
      
      selected_item <- sample(
        remaining_items,
        size = 1,
        prob = probabilities
      )
      
      complete_rankings[r, k] <- selected_item
      
      remaining_items <- setdiff(
        remaining_items,
        selected_item
      )
    }
  }
  
  as.data.frame(complete_rankings)
}


# =========================================================
# 3. Convert character rankings into numeric rankings
# =========================================================

convert_to_numeric <- function(complete_rankings) {
  
  items <- c("A1", "A2", "A3", "A4", "A5")
  
  complete_rankings_num <- as.matrix(
    data.frame(
      lapply(
        complete_rankings,
        function(x) match(x, items)
      )
    )
  )
  
  complete_rankings_num
}


# =========================================================
# 4. Convert complete rankings into partial rankings
# =========================================================

make_partial_data <- function(complete_rankings_num) {
  
  partial_output <- make_partial(
    data = complete_rankings_num,
    format_input = "ordering",
    probcens = c(
      0.25,
      0.25,
      0.25,
      0.25
    )
  )
  
  partial_output
}


# =========================================================
# 5. Fit Bayesian Plackett-Luce model
# =========================================================

fit_bayesian_pl <- function(partial_rankings_num) {
  
  partial_rankings_plmix <- as.top_ordering(
    data = partial_rankings_num,
    format_input = "ordering",
    aggr = FALSE
  )
  
  bayesian_pl_fit <- gibbsPLMIX(
    pi_inv = partial_rankings_plmix,
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
  )
  
  posterior_mean <- colMeans(
    bayesian_pl_fit$P
  )
  
  estimated_order_num <- order(
    posterior_mean,
    decreasing = TRUE
  )
  
  estimated_order <- paste0(
    "A",
    estimated_order_num
  )
  
  list(
    posterior_mean = posterior_mean,
    estimated_order = estimated_order
  )
}


# =========================================================
# 6. Kendall distance
# =========================================================

kendall_distance <- function(
    estimated_order,
    true_order) {
  
  item_pairs <- combn(
    seq_along(true_order),
    2
  )
  
  estimated_position <- match(
    true_order,
    estimated_order
  )
  
  sum(
    estimated_position[item_pairs[1, ]] >
      estimated_position[item_pairs[2, ]]
  )
}


# =========================================================
# 7. Simulation
# =========================================================

sample_sizes <- c(20, 50, 100)
B <- 30

pl_bayesian_results <- data.frame()

set.seed(123)

for (score_name in names(score_list)) {
  
  w_true <- score_list[[score_name]]
  
  for (N in sample_sizes) {
    
    recovery <- numeric(B)
    kendall <- numeric(B)
    
    cat(
      "\nPL mechanism:",
      score_name,
      "N =", N,
      "\n"
    )
    
    for (b in 1:B) {
      
      # Generate complete rankings
      
      complete_rankings <- generate_pl_data(
        w_true = w_true,
        N = N
      )
      
      complete_rankings_num <- convert_to_numeric(
        complete_rankings
      )
      
      
      # Convert to partial rankings
      
      partial_output <- make_partial_data(
        complete_rankings_num
      )
      
      partial_rankings <- partial_output$partialdata
      
      
      # Fit Bayesian PL
      
      bayesian_result <- fit_bayesian_pl(
        partial_rankings
      )
      
      estimated_order <-
        bayesian_result$estimated_order
      
      
      # Exact recovery
      
      recovery[b] <- as.integer(
        identical(
          estimated_order,
          true_order
        )
      )
      
      
      # Kendall distance
      
      kendall[b] <- kendall_distance(
        estimated_order,
        true_order
      )
      
      
      # Progress
      
      cat(
        "Replication",
        b,
        "of",
        B,
        "\r"
      )
    }
    
    
    # Store results
    
    pl_bayesian_results <- rbind(
      pl_bayesian_results,
      data.frame(
        separation = score_name,
        N = N,
        recovery_rate = mean(recovery),
        mean_kendall = mean(kendall)
      )
    )
    
    cat("\nCompleted.\n")
  }
}


# =========================================================
# 8. Results
# =========================================================

pl_bayesian_results