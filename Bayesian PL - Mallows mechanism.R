library(PLMIX)
library(PerMallows)

# =========================================================
# 1. Simulation settings
# =========================================================

alpha_list <- list(
  high = 0.45,
  medium = 0.15,
  low = 0.025
)

true_order <- c("A1", "A2", "A3", "A4", "A5")

central_ranking <- c(1, 2, 3, 4, 5)


# =========================================================
# 2. Generate complete rankings under Mallows mechanism
# =========================================================

generate_mallows_data <- function(alpha, N) {
  
  complete_rankings_num <- rmm(
    n = N,
    sigma0 = central_ranking,
    theta = alpha,
    dist.name = "kendall"
  )
  
  complete_rankings_num
}


# =========================================================
# 3. Convert complete rankings into partial rankings
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
# 4. Fit Bayesian Plackett-Luce model
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
# 5. Kendall distance
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
# 6. Simulation
# =========================================================

sample_sizes <- c(20, 50, 100)
B <- 30

mm_bayesian_results <- data.frame()

set.seed(123)

for (separation_name in names(alpha_list)) {
  
  alpha <- alpha_list[[separation_name]]
  
  for (N in sample_sizes) {
    
    recovery <- numeric(B)
    kendall <- numeric(B)
    
    cat(
      "\nMallows mechanism:",
      separation_name,
      "N =", N,
      "\n"
    )
    
    for (b in 1:B) {
      
      # Generate complete rankings
      
      complete_rankings_num <- generate_mallows_data(
        alpha = alpha,
        N = N
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
    
    mm_bayesian_results <- rbind(
      mm_bayesian_results,
      data.frame(
        separation = separation_name,
        N = N,
        recovery_rate = mean(recovery),
        mean_kendall = mean(kendall)
      )
    )
    
    cat("\nCompleted.\n")
  }
}


# =========================================================
# 7. Results
# =========================================================

mm_bayesian_results