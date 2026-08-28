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
  
  # Calculate the length of each ranking
  ranking_lengths <- sapply(orders, length)
  
  # Sum the number of voters for each ranking length
  length_counts <- tapply(counts, ranking_lengths, sum)
  
  # Convert the results into a data frame
  length_table <- data.frame(
    RankingLength = as.integer(names(length_counts)),
    NumberOfVoters = as.integer(length_counts)
  )
  
  # Plot
  bp <- barplot(
    height = length_table$NumberOfVoters,
    names.arg = length_table$RankingLength,
    xlab = "Ranking length",
    ylab = "Number of voters",
    ylim = c(0, max(length_table$NumberOfVoters) * 1.1) 
  )
  
  # Add the number of voters above each bar
  text(
    x = bp,
    y = length_table$NumberOfVoters,
    labels = length_table$NumberOfVoters,
    pos = 3,
    cex = 0.8, 
    col = "black"
  )
  
  
  # Pairwise win matrix
  n <- 5
  win_mat <- matrix(0, nrow = n, ncol = n)
  comp_mat <- matrix(0, nrow = n, ncol = n)
  
  # Convert ranking data into pairwise win counts
  for (s in 1:length(orders)) {
    
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
        
        # Total number of comparisons between i and j
        comp_mat[i,j] <- comp_mat[i,j] + count
        comp_mat[j,i] <- comp_mat[j,i] + count
         
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
  
  # Pairwise win probabilities
  prob_mat <- win_mat / comp_mat
  library(pheatmap)
  
  # Diagonal set to 0.5, because self-comparison no meaning
  diag(prob_mat) <- 0.5
  
  rownames(prob_mat) <- paste("Cand", 1:n)
  colnames(prob_mat) <- paste("Cand", 1:n)

  pheatmap(
    prob_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("firebrick3", "white", "steelblue4"))(100),
    breaks = seq(0, 1, length.out = 101),
    display_numbers = FALSE,
    legend = TRUE,
    border_color = "grey80"
  )
