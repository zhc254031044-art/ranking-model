library(PerMallows)
library(PlackettLuce)

# 1. Read APA data

lines <- readLines("data_analysis/APA.soi")

# Remove comments and empty lines
data_lines <- lines[!grepl("^#", lines)]
data_lines <- data_lines[nchar(data_lines) > 0]

# Extract the number of voters and ranking information
counts <- as.integer(sub(":.*", "", data_lines))
orders_text <- trimws(sub("^[0-9]+:", "", data_lines))

# Split each ranking into individual items and convert them to integers
orders <- strsplit(orders_text, ",")
orders <- lapply(orders, function(x) as.integer(trimws(x)))

# 2. Keep complete rankings only, Top-4 rankings are also treated as complete rankings

ranking_lengths <- sapply(orders, length)
complete_index <- ranking_lengths >= 4

complete_orders <- orders[complete_index]
complete_counts <- counts[complete_index]

# Add missing candidate for top-4 rankings

for (i in 1:length(complete_orders)) {
  ranking <- complete_orders[[i]]
  # If only four candidates are ranked, add the remaining candidate as the last position
  if (length(ranking) == 4) {
    missing_candidate <- setdiff(1:5, ranking)
    complete_orders[[i]] <- c(ranking, missing_candidate)
  }
}

# Expand ranking patterns according to voter counts

# Create an empty matrix to store individual complete rankings
complete_rankings_num <- matrix(nrow = sum(complete_counts),ncol = 5)
row_index <- 1

# Add each ranking pattern according to its corresponding voter count
for (i in 1:length(complete_orders)) {
  
  ranking <- complete_orders[[i]]
  count <- complete_counts[i]
  
  # Repeat the same ranking pattern for all voters with this ranking
  complete_rankings_num[row_index:(row_index + count - 1),] <- matrix(
    rep(ranking, times = count),
    nrow = count,
    ncol = 5,
    byrow = TRUE
  )
  
  # Update the starting row for the next ranking pattern
  row_index <- row_index + count
}

# 3. Plackett-Luce model

# Convert numeric rankings into candidate labels required by PlackettLuce
complete_rankings <- as.data.frame(matrix(
    paste0("A", complete_rankings_num),
    nrow = nrow(complete_rankings_num),
    ncol = ncol(complete_rankings_num)
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

# 4. Mallows model

mm_fit <- lmm(
  data = complete_rankings_num,
  dist.name = "kendall"
  #dist.name = "cayley"
)

# Extract estimated consensus ranking and concentration parameter
mm_mode <- mm_fit$mode
mm_theta <- mm_fit$theta

# 5. Results

pl_scores
pl_ranking

mm_mode
mm_theta

# Count how many rankings place A3 before A4
a3_before_a4 <- 0

for (i in 1:nrow(complete_rankings_num)) {
  
  ranking <- complete_rankings_num[i, ]
  
  if (which(ranking == 3) < which(ranking == 4)) {
    a3_before_a4 <- a3_before_a4 + 1
  }
}

a4_before_a3 <- nrow(complete_rankings_num) - a3_before_a4

a3_before_a4
a4_before_a3
