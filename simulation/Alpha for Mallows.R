# PL expected Kendall distance
pl_expected_kendall <- function(w) {
  total <- 0
  for (i in 1:(length(w) - 1)) {
    for (j in (i + 1):length(w)) {
      total <- total + w[j] / (w[i] + w[j])
    }
  }
  total
}

# Mallows expected Kendall distance
mallows_expected_kendall <- function(alpha, n = 5) {
  q <- exp(-alpha)
  j <- 1:n
  
  n * q / (1 - q) - sum(j * q^j / (1 - q^j))
}


# Find alpha by matching expected Kendall distances
find_alpha <- function(w) {
  target <- pl_expected_kendall(w)
  
  alpha <- uniroot(
    function(a) mallows_expected_kendall(a) - target,
    c(0.000001, 10)
  )$root
  alpha
}


# Three score settings
w_high <- c(5, 4, 3, 2, 1)
w_medium <- c(5, 4.5, 4, 3.5, 3)
w_low <- c(5, 4.9, 4.8, 4.7, 4.6)

find_alpha(w_high)
find_alpha(w_medium)
find_alpha(w_low)
