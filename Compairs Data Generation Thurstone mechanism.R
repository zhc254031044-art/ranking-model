library(eba)
library(BradleyTerry2)


# 1. Simulation settings

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
  item_pairs <- combn(items, 2)
  
  pairwise_matrix <- matrix(
    0,
    nrow = length(items),
    ncol = length(items),
    dimnames = list(items, items)
  )
  
  utility <- log(w_true)
  
  for (k in 1:ncol(item_pairs)) {
    
    item_i <- item_pairs[1, k]
    item_j <- item_pairs[2, k]
    
    p_i_wins <- pnorm(
      utility[item_i] - utility[item_j]
    )
    
    win_i <- rbinom(1, size = N, prob = p_i_wins)
    
    pairwise_matrix[item_i, item_j] <- win_i
    pairwise_matrix[item_j, item_i] <- N - win_i
  }
  
  pairwise_matrix
}


# 3. Fit the Thurstone model

fit_thurstone <- function(pairwise_matrix) {
  
  thurstone_fit <- thurstone(pairwise_matrix)
  thurstone_scores <- thurstone_fit$estimate
  
  thurstone_result <- data.frame(
    item = names(thurstone_scores),
    score = as.numeric(thurstone_scores)
  )
  
  thurstone_result <- thurstone_result[
    order(thurstone_result$score, decreasing = TRUE),
  ]
  
  rownames(thurstone_result) <- NULL
  
  thurstone_result
}


# 4. Fit the Bradley-Terry model

fit_bt <- function(pairwise_matrix) {
  
  items <- rownames(pairwise_matrix)
  item_pairs <- combn(items, 2)
  
  pairwise_data <- data.frame(
    item1 = item_pairs[1, ],
    item2 = item_pairs[2, ],
    win1 = pairwise_matrix[
      cbind(item_pairs[1, ], item_pairs[2, ])
    ],
    win2 = pairwise_matrix[
      cbind(item_pairs[2, ], item_pairs[1, ])
    ]
  )
  
  pairwise_data$item1 <- factor(
    pairwise_data$item1,
    levels = items
  )
  
  pairwise_data$item2 <- factor(
    pairwise_data$item2,
    levels = items
  )
  
  bt_fit <- BTm(
    outcome = cbind(win1, win2),
    player1 = item1,
    player2 = item2,
    formula = ~ item,
    id = "item",
    data = pairwise_data
  )
  
  bt_ability <- BTabilities(bt_fit)
  
  bt_result <- data.frame(
    item = rownames(bt_ability),
    ability = bt_ability[, "ability"]
  )
  
  bt_result <- bt_result[
    order(bt_result$ability, decreasing = TRUE),
  ]
  
  rownames(bt_result) <- NULL
  
  bt_result
}


# 5. Kendall distance

kendall_distance <- function(estimated_order, true_order) {
  
  item_pairs <- combn(seq_along(true_order), 2)
  estimated_position <- match(true_order, estimated_order)
  
  sum(
    estimated_position[item_pairs[1, ]] >
      estimated_position[item_pairs[2, ]]
  )
}


# 6. Repeat one simulation setting

run_thurstone_simulation <- function(w_true, N, B) {
  
  true_order <- names(sort(w_true, decreasing = TRUE))
  
  simulation_results <- data.frame(
    thurstone_recovery = logical(B),
    bt_recovery = logical(B),
    thurstone_distance = numeric(B),
    bt_distance = numeric(B)
  )
  
  for (b in 1:B) {
    
    pairwise_matrix <- generate_thurstone_data(
      w_true,
      N
    )
    
    thurstone_order <- fit_thurstone(
      pairwise_matrix
    )$item
    
    bt_order <- fit_bt(
      pairwise_matrix
    )$item
    
    simulation_results$thurstone_recovery[b] <-
      identical(thurstone_order, true_order)
    
    simulation_results$bt_recovery[b] <-
      identical(bt_order, true_order)
    
    simulation_results$thurstone_distance[b] <-
      kendall_distance(
        thurstone_order,
        true_order
      )
    
    simulation_results$bt_distance[b] <-
      kendall_distance(
        bt_order,
        true_order
      )
  }
  
  simulation_results
}


# 7. Run all Thurstone-generated settings

set.seed(123)

thurstone_simulation_summary <- data.frame()

for (score_name in names(score_list)) {
  
  for (N in sample_sizes) {
    
    results <- run_thurstone_simulation(
      w_true = score_list[[score_name]],
      N = N,
      B = B
    )
    
    setting_summary <- data.frame(
      separation = score_name,
      sample_size = N,
      thurstone_recovery_rate =
        mean(results$thurstone_recovery),
      bt_recovery_rate =
        mean(results$bt_recovery),
      thurstone_mean_kendall =
        mean(results$thurstone_distance),
      bt_mean_kendall =
        mean(results$bt_distance)
    )
    
    thurstone_simulation_summary <- rbind(
      thurstone_simulation_summary,
      setting_summary
    )
  }
}

thurstone_simulation_summary