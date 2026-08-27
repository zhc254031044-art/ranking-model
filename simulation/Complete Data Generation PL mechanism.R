library(PerMallows)
library(PlackettLuce)

# 1. Simulation settings

# True score settings
score_list <- list(
  high = c(A1 = 5, A2 = 4, A3 = 3, A4 = 2, A5 = 1),
  medium = c(A1 = 5, A2 = 4.5, A3 = 4, A4 = 3.5, A5 = 3),
  low = c(A1 = 5, A2 = 4.9, A3 = 4.8, A4 = 4.7, A5 = 4.6)
)

sample_sizes <- c(100, 500, 1000)

B <- 100

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


# 3. Fit the Plackett-Luce model

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


# 4. Fit the Mallows model

fit_mallows <- function(complete_rankings) {
  
  mm_fit <- lmm(
    data = complete_rankings,
    dist.name = "kendall"
  )
  
  # Extract estimated consensus ranking
  mm_mode <- mm_fit$mode
  
  # Convert numeric ranking back to candidate labels
  mm_ranking <- paste0("A", mm_mode)
  mm_ranking
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
    
    # Initialise logical vectors of length B with FALSE values
    pl_recovery <- logical(B)
    mallows_recovery <- logical(B)
    
    # Initialise numeric vectors of length B with zero values
    pl_distance <- numeric(B)
    mallows_distance <- numeric(B)
    
    for (b in 1:B) {
      
      # Generate complete ranking data
      complete_rankings <- generate_pl_data(w_true, N)
      
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
    cat("Separation:", score_name, "\n")
    cat("Sample size:", N, "\n")
    cat("PL recovery rate:", mean(pl_recovery), "\n")
    cat("Mallows recovery rate:", mean(mallows_recovery), "\n")
    cat("PL mean Kendall distance:", mean(pl_distance), "\n")
    cat("Mallows mean Kendall distance:", mean(mallows_distance), "\n")
  }
}