library(PlackettLuce)
library(PLMIX)

# 1. Read APA data

lines <- readLines("D:/KCL/ranking/code/Empirical Analysis/APA.soi")

# Remove comments and empty lines
data_lines <- lines[!grepl("^#", lines)]
data_lines <- data_lines[nchar(data_lines) > 0]

# Extract the number of voters and ranking information
counts <- as.integer(sub(":.*", "", data_lines))
orders_text <- trimws(sub("^[0-9]+:", "", data_lines))

# Split each ranking into individual items and convert them to integers
orders <- strsplit(orders_text, ",")
orders <- lapply(orders, function(x) as.integer(trimws(x)))

# 2. Expand rankings according to voter counts

# Create an empty matrix to store individual complete rankings
apa_rankings_num <- matrix(nrow = sum(counts),ncol = 5)
row_index <- 1

# Add each ranking pattern according to its corresponding voter count
for (i in seq_along(orders)) {
  ranking <- orders[[i]]
  count <- counts[i]
  
  # Add NA for unranked candidates to make all rankings have length 5
  ranking_full <- c(ranking, rep(NA, 5 - length(ranking)))
  
  # Repeat the same ranking pattern for all voters with this ranking
  apa_rankings_num[row_index:(row_index + count - 1), ] <- matrix(
    rep(ranking_full, times = count),
    nrow = count,
    ncol = 5,
    byrow = TRUE
  )
  
  # Update the starting row for the next ranking pattern
  row_index <- row_index + count
}

# 3. Plackett-Luce model

# Convert candidate numbers into candidate labels required by PlackettLuce
apa_rankings <- apa_rankings_num
for (i in 1:nrow(apa_rankings)) {
  for (j in 1:ncol(apa_rankings)) {
    if (!is.na(apa_rankings[i,j])) {
      apa_rankings[i,j] <- paste0("A", apa_rankings[i,j])
    }
  }
}

apa_rankings <- as.data.frame(apa_rankings)
colnames(apa_rankings) <- paste0("rank_",1:5)

# Convert ordering data into Plackett-Luce ranking object
pl_rankings <- as.rankings(apa_rankings, input = "orderings")
pl_fit <- PlackettLuce(pl_rankings, npseudo = 0)

# Extract estimated worth parameters
pl_scores <- coef(pl_fit, log = FALSE)

# Rank candidates according to estimated scores
pl_ranking <- names(sort(pl_scores, decreasing = TRUE))

# 4. Bayesian Plackett-Luce model

# Convert NA values to 0 for unranked candidates
bayesian_rankings_num <- apa_rankings_num
bayesian_rankings_num[is.na(bayesian_rankings_num)] <- 0

# Convert ranking data into PLMIX top-ordering format
bayesian_rankings <- as.top_ordering(
  data = bayesian_rankings_num,
  format_input = "ordering",
  aggr = FALSE
)

# MCMC settings
n_iter <- 2000
n_burn <- 400

# Fit Bayesian Plackett--Luce model
bayesian_pl_fit <- gibbsPLMIX(
  pi_inv = bayesian_rankings,
  K = 5,
  G = 1,
  n_iter = n_iter,
  n_burn = n_burn,
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

# 5. Bayesian posterior results

# Posterior draws of the score parameters
posterior_scores <- bayesian_pl_fit$P

# Posterior means
bayesian_scores <- colMeans( posterior_scores)

names(bayesian_scores) <- paste0("A",1:5)

# Estimated ranking
bayesian_ranking <- names(sort( bayesian_scores, decreasing = TRUE))

# 6. Results

pl_scores
pl_ranking

bayesian_scores
bayesian_ranking

