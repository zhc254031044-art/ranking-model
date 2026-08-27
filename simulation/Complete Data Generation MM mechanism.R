library(PerMallows)
library(PlackettLuce)


# 1. Simulation settings

alpha_list <- list(
  high = 0.45,
  medium = 0.15,
  low = 0.025
)

sample_sizes <- c(100, 500, 1000)

B <- 100

central_ranking <- c(1, 2, 3, 4, 5)


# 2. Generate complete ranking data under the Mallows mechanism

generate_mallows_data <- function(alpha, N) {
  
  complete_rankings_num <- rmm(
    n = N,
    sigma0 = central_ranking,
    theta = alpha,
    dist.name = "kendall"
  )
  
  complete_rankings_num
}


# 3. Fit the Mallows model

fit_mallows <- function(complete_rankings_num) {
  
  mm_fit <- lmm(
    data = complete_rankings_num,
    dist.name = "kendall"
  )
  
  # Extract estimated consensus ranking
  mm_mode <- mm_fit$mode
  
  # Convert numeric ranking back to candidate labels
  mm_ranking <- paste0("A", mm_mode)
  mm_ranking
}


# 4. Fit the Plackett-Luce model

fit_pl <- function(complete_rankings) {
  
  # Convert numeric rankings into candidate labels required by PlackettLuce
  complete_rankings <- as.data.frame(matrix(
    paste0("A", complete_rankings),
    nrow = nrow(complete_rankings),
    ncol = ncol(complete_rankings)
  )
  )
  
  colnames(complete_rankings) <- paste0("rank_", 1:5)
  
  # Convert ordering data into Plackett-Luce ranking object
  pl_rankings <- as.rankings(complete_rankings, input = "orderings")
  pl_fit <- PlackettLuce(pl_rankings, npseudo = 0)
  
  # Extract estimated worth parameters
  pl_scores <- coef(pl_fit, log = FALSE)
  
  # Rank candidates according to estimated scores
  pl_ranking <- names(sort(pl_scores, decreasing = TRUE))
  pl_ranking
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
    
    # Initialise logical vectors of length B with FALSE values
    pl_recovery <- logical(B)
    mallows_recovery <- logical(B)
    
    # Initialise numeric vectors of length B with zero values
    pl_distance <- numeric(B)
    mallows_distance <- numeric(B)
    
    for (b in 1:B) {
      
      # Generate complete ranking data
      complete_rankings <- generate_mallows_data(alpha, N)
      
      # Fit the two models
      pl_order <- fit_pl(complete_rankings)
      mallows_order <- fit_mallows(complete_rankings)
      
      # Calculate Kendall distance
      pl_distance[b] <- kendall_distance(pl_order, true_order)
      mallows_distance[b] <- kendall_distance(mallows_order, true_order)
      
      # TRUE if the estimated ranking exactly matches the true ranking
      pl_recovery[b] <- pl_distance[b] == 0
      mallows_recovery[b] <- mallows_distance[b] == 0
    }
    
    # Print results for this setting
    cat("\n")
    cat("Separation:", name, "\n")
    cat("Sample size:", N, "\n")
    cat("PL recovery rate:", mean(pl_recovery), "\n")
    cat("Mallows recovery rate:", mean(mallows_recovery), "\n")
    cat("PL mean Kendall distance:", mean(pl_distance), "\n")
    cat("Mallows mean Kendall distance:", mean(mallows_distance), "\n")
  }
}