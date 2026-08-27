library(eba)
library(BradleyTerry2)

# 1. Simulation settings

# True score settings
score_list <- list(
  high = c(A1 = 5, A2 = 4, A3 = 3, A4 = 2, A5 = 1),
  medium = c(A1 = 5, A2 = 4.5, A3 = 4, A4 = 3.5, A5 = 3),
  low = c(A1 = 5, A2 = 4.9, A3 = 4.8, A4 = 4.7, A5 = 4.6)
)

sample_sizes <- c(100, 500, 1000)

B <- 100

# 2. Generate Thurstone pairwise data

generate_thurstone_data <- function(w_true, N) {
  
  items <- names(w_true)
  
  # Extract all unique pairwise combinations
  item_pairs <- combn(items, 2)
  
  # Create a 5*5 win matrix
  pairwise_matrix <- matrix(0, nrow = length(items), ncol = length(items), dimnames = list(items, items))
  
  # Convert true scores to utilities
  utility <- log(w_true)
  
  # Loop through each pair of items
  for (k in 1:ncol(item_pairs)) {
    
    # Extract the two items in the current pair
    item_i <- item_pairs[1, k]
    item_j <- item_pairs[2, k]
    
    # Calculate the Thurstone probability that item i beats item j
    p_i_wins <- pnorm(utility[item_i] - utility[item_j])
    
    # Simulate the number of times item i wins
    win_i <- rbinom(1, size = N, prob = p_i_wins)
    
    # Store the win counts for both directions
    pairwise_matrix[item_i, item_j] <- win_i
    pairwise_matrix[item_j, item_i] <- N - win_i
  }
  pairwise_matrix
}

# 3. Fit the Thurstone model

fit_thurstone <- function(pairwise_matrix) {
  scores <- thurstone(pairwise_matrix)$estimate
  names(sort(scores, decreasing = TRUE))
}

# 4. Fit the Bradley-Terry model

fit_bt <- function(pairwise_matrix) {
  
  n <- nrow(pairwise_matrix)
  bt_data <- data.frame()
  
  # Convert the win matrix into pairwise comparison data required by the BT model
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      bt_data <- rbind(bt_data, data.frame(
        player1 = paste0("A", i),
        player2 = paste0("A", j),
        win1 = pairwise_matrix[i, j],
        win2 = pairwise_matrix[j, i]))
    }
  }
  
  # Define common candidate levels for both comparison variables
  item_levels <- paste0("A", 1:n)
  
  bt_data$player1 <- factor(
    bt_data$player1,
    levels = item_levels
  )
  
  bt_data$player2 <- factor(
    bt_data$player2,
    levels = item_levels
  )
  
  bt_fit <- BTm(
    cbind(win1, win2),
    player1,
    player2,
    data = bt_data
  )
  
  bt_scores <- BTabilities(bt_fit)
  bt_ranking <- rownames(bt_scores)[order(bt_scores[, "ability"], decreasing = TRUE)]
  bt_ranking
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
    thurstone_recovery <- logical(B)
    bt_recovery <- logical(B)

    # Initialise numeric vectors of length B with zero values
    thurstone_distance <- numeric(B)
    bt_distance <- numeric(B)
    
    for (b in 1:B) {
      
      # Generate pairwise data
      pairwise_matrix <- generate_thurstone_data(w_true, N)
      
      # Fit the two models
      thurstone_order <- fit_thurstone(pairwise_matrix)
      bt_order <- fit_bt(pairwise_matrix)
      
      # Calculate Kendall distance
      thurstone_distance[b] <- kendall_distance(thurstone_order, true_order)
      bt_distance[b] <- kendall_distance(bt_order, true_order)
      
      # TRUE if the estimated ranking exactly matches the true ranking
      thurstone_recovery[b] <- thurstone_distance[b] == 0
      bt_recovery[b] <- bt_distance[b] == 0
    }
    
    # Print results for this setting
    cat("\n")
    cat("Separation:", score_name, "\n")
    cat("Sample size:", N, "\n")
    cat("Thurstone recovery rate:", mean(thurstone_recovery), "\n" )
    cat("BT recovery rate:", mean(bt_recovery), "\n" )
    cat("Thurstone mean Kendall distance:", mean(thurstone_distance), "\n" )
    cat("BT mean Kendall distance:", mean(bt_distance), "\n")
  }
}
