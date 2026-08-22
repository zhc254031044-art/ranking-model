library(eba)
library(BradleyTerry2)

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

# 2. Construct pairwise win matrix

n <- 5

# Create a 5*5 win matrix
win_mat <- matrix(0, nrow = n, ncol = n, dimnames = list( paste0("A", 1:n), paste0("A", 1:n)))

# Convert ranking data into pairwise win counts
for (s in seq_along(orders)) {
  
  ranking <- orders[[s]]
  count <- counts[s]
  
  # Record the ranking position of each candidate
  pos <- match(1:n, ranking)
  
  # Compare every pair of candidates
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      
      # Check whether each candidate is ranked in this vote
      i_ranked <- !is.na(pos[i])
      j_ranked <- !is.na(pos[j])
      
      # No comparison if both candidates are unranked
      if (!i_ranked && !j_ranked) {
        next
      }
      
      # Both candidates are ranked
      if (i_ranked && j_ranked) {
        if (pos[i] < pos[j]) {
          win_mat[i, j] <- win_mat[i, j] + count
        } else {
          win_mat[j, i] <- win_mat[j, i] + count
        }
      }
      
      # Only i is ranked
      else if (i_ranked && !j_ranked) {
        win_mat[i, j] <- win_mat[i, j] + count
      }
      
      # Only j is ranked
      else if (!i_ranked && j_ranked) {
        win_mat[j, i] <- win_mat[j, i] + count
      }
    }
  }
}

# 3. Thurstone model fitting

thurstone_fit <- thurstone(win_mat)
thurstone_scores <- thurstone_fit$estimate
thurstone_ranking <- names(sort(thurstone_scores, decreasing = TRUE))

# 4. Bradley-Terry model fitting

bt_data <- data.frame()

# Convert the win matrix into pairwise comparison data required by the BT model
for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    
    bt_data <- rbind(bt_data, data.frame(
        player1 = paste0("A", i),
        player2 = paste0("A", j),
        win1 = win_mat[i, j],
        win2 = win_mat[j, i]))
  }
}

# Define common candidate levels for both comparison variables
item_levels <- paste0("A", 1:n)

# Convert candidate labels into factors
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

# 5. Results

win_mat
thurstone_scores
bt_scores

thurstone_ranking
bt_ranking