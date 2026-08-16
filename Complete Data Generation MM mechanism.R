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

true_order <- c("A1", "A2", "A3", "A4", "A5")

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
  
  mm_mode <- mm_fit$mode
  
  mallows_order <- paste0(
    "A",
    mm_mode
  )
  
  return(mallows_order)
}


# 4. Fit the Plackett-Luce model

fit_pl <- function(complete_rankings_num) {
  
  complete_rankings <- as.data.frame(
    matrix(
      paste0("A", complete_rankings_num),
      nrow = nrow(complete_rankings_num),
      ncol = ncol(complete_rankings_num)
    )
  )
  
  colnames(complete_rankings) <- paste0(
    "rank_",
    1:5
  )
  
  pl_rankings <- as.rankings(
    complete_rankings,
    input = "orderings"
  )
  
  pl_fit <- PlackettLuce(
    pl_rankings,
    npseudo = 0
  )
  
  pl_score <- coef(
    pl_fit,
    log = FALSE
  )
  
  pl_order <- names(
    sort(pl_score, decreasing = TRUE)
  )
  
  return(pl_order)
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

run_mallows_simulation <- function(alpha, N, B) {
  
  simulation_results <- data.frame(
    mallows_recovery = logical(B),
    pl_recovery = logical(B),
    mallows_distance = numeric(B),
    pl_distance = numeric(B)
  )
  
  for (b in 1:B) {
    
    complete_rankings_num <- generate_mallows_data(
      alpha,
      N
    )
    
    mallows_order <- fit_mallows(
      complete_rankings_num
    )
    
    pl_order <- fit_pl(
      complete_rankings_num
    )
    
    simulation_results$mallows_recovery[b] <-
      identical(
        mallows_order,
        true_order
      )
    
    simulation_results$pl_recovery[b] <-
      identical(
        pl_order,
        true_order
      )
    
    simulation_results$mallows_distance[b] <-
      kendall_distance(
        mallows_order,
        true_order
      )
    
    simulation_results$pl_distance[b] <-
      kendall_distance(
        pl_order,
        true_order
      )
  }
  
  simulation_results
}


# 7. Run all Mallows-generated settings

set.seed(123)

mallows_simulation_summary <- data.frame()

for (setting_name in names(alpha_list)) {
  
  for (N in sample_sizes) {
    
    results <- run_mallows_simulation(
      alpha = alpha_list[[setting_name]],
      N = N,
      B = B
    )
    
    setting_summary <- data.frame(
      separation = setting_name,
      sample_size = N,
      mallows_recovery_rate =
        mean(results$mallows_recovery),
      pl_recovery_rate =
        mean(results$pl_recovery),
      mallows_mean_kendall =
        mean(results$mallows_distance),
      pl_mean_kendall =
        mean(results$pl_distance)
    )
    
    mallows_simulation_summary <- rbind(
      mallows_simulation_summary,
      setting_summary
    )
  }
}

mallows_simulation_summary