library(PlackettLuce)
library(PerMallows)


# 1. Simulation settings

score_list <- list(
  high = c(A1 = 5, A2 = 4, A3 = 3, A4 = 2, A5 = 1),
  medium = c(A1 = 5, A2 = 4.5, A3 = 4, A4 = 3.5, A5 = 3),
  low = c(A1 = 5, A2 = 4.9, A3 = 4.8, A4 = 4.7, A5 = 4.6)
)

sample_sizes <- c(100, 500, 1000)

B <- 100


# 2. Generate Plackett-Luce complete ranking data

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
      
      probabilities <- w_true[remaining_items] /
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
  
  complete_rankings <- as.data.frame(
    complete_rankings
  )
  
  complete_rankings
}


# 3. Fit the Plackett-Luce model

fit_pl <- function(complete_rankings) {
  
  pl_rankings <- as.rankings(
    complete_rankings,
    input = "orderings"
  )
  
  pl_fit <- PlackettLuce(
    pl_rankings,
    npseudo = 0
  )
  
  pl_worth <- coef(
    pl_fit,
    log = FALSE
  )
  
  pl_result <- data.frame(
    item = names(pl_worth),
    worth = as.numeric(pl_worth)
  )
  
  pl_result <- pl_result[
    order(pl_result$worth, decreasing = TRUE),
  ]
  
  rownames(pl_result) <- NULL
  
  pl_result
}


# 4. Fit the Mallows model

fit_mallows <- function(complete_rankings) {
  
  items <- names(w_true)
  
  complete_rankings_num <- as.matrix(
    data.frame(
      lapply(
        complete_rankings,
        function(x) match(x, items)
      )
    )
  )
  
  
  colnames(complete_rankings_num) <- paste0(
    "rank_",
    1:length(items)
  )
  
  mm_fit <- lmm(
    data = complete_rankings_num,
    dist.name = "kendall"
  )
  
  mm_mode <- mm_fit$mode
  
  mm_result <- data.frame(
    rank = 1:length(mm_mode),
    item = items[mm_mode]
  )
  
  mm_result
}


# 5. Kendall distance

kendall_distance <- function(estimated_order, true_order) {
  
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


# 6. Repeat one simulation setting

run_pl_simulation <- function(w_true, N, B) {
  
  true_order <- names(
    sort(w_true, decreasing = TRUE)
  )
  
  simulation_results <- data.frame(
    pl_recovery = logical(B),
    mallows_recovery = logical(B),
    pl_distance = numeric(B),
    mallows_distance = numeric(B)
  )
  
  for (b in 1:B) {
    
    complete_rankings <- generate_pl_data(
      w_true,
      N
    )
    
    pl_order <- fit_pl(
      complete_rankings
    )$item
    
    mallows_order <- fit_mallows(
      complete_rankings
    )$item
    
    simulation_results$pl_recovery[b] <-
      identical(pl_order, true_order)
    
    simulation_results$mallows_recovery[b] <-
      identical(mallows_order, true_order)
    
    simulation_results$pl_distance[b] <-
      kendall_distance(
        pl_order,
        true_order
      )
    
    simulation_results$mallows_distance[b] <-
      kendall_distance(
        mallows_order,
        true_order
      )
  }
  
  simulation_results
}


# 7. Run all Plackett-Luce-generated settings

set.seed(123)

pl_simulation_summary <- data.frame()

for (score_name in names(score_list)) {
  
  for (N in sample_sizes) {
    
    results <- run_pl_simulation(
      w_true = score_list[[score_name]],
      N = N,
      B = B
    )
    
    setting_summary <- data.frame(
      separation = score_name,
      sample_size = N,
      pl_recovery_rate =
        mean(results$pl_recovery),
      mallows_recovery_rate =
        mean(results$mallows_recovery),
      pl_mean_kendall =
        mean(results$pl_distance),
      mallows_mean_kendall =
        mean(results$mallows_distance)
    )
    
    pl_simulation_summary <- rbind(
      pl_simulation_summary,
      setting_summary
    )
  }
}

pl_simulation_summary