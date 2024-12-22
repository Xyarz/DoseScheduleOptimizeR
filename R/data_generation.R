
#' @title generate_full_factorial_design
#'
#' @param doses the dose levels
#' @param schedules the schedules
#' @param n_pat sample size to be allocated across the factorial space
#'
#' @return matrix
#' @export
generate_full_factorial_design <- function(
    doses,
    schedules,
    n_pat
    ) {
  create_matrix <- function(
    nrow = doses,
    ncol = schedules,
    n = n_pat
    ) {
    # Create an empty matrix with the given dimensions
    matrix <- matrix(0, length(nrow), length(ncol))

    # Calculate the base value for each cell and the remainder
    base_value <- floor(n / (length(nrow) * length(ncol)))
    remainder <- n - (base_value * length(nrow) * length(ncol))

    # Fill the matrix with the base value
    matrix[] <- base_value

    # Distribute the remainder randomly across the cells
    if (remainder > 0) {
      # Create a vector of indices for the cells in the matrix
      indices <- 1:(length(nrow) * length(ncol))

      # Randomly select 'remainder' number of indices
      selected_indices <- sample(indices, remainder)

      # Add 1 to the selected cells
      matrix[selected_indices] <- matrix[selected_indices] + 1
    }

    colnames(matrix) <- ncol
    rownames(matrix) <- nrow
    return(matrix)
  }
  return(create_matrix())
}

#' @title generate_optimal_design_linear
#'
#' @param doses doses
#' @param schedules schedules
#' @param n_pat sample size to be allocated across the factorial space
#' @param criterion optimality criteria, default "D"
#' @param random_order boolean, default TRUE
#' @param interaction boolean, default TRUE
#'
#' @return matrix
#' @export
generate_optimal_design_linear <- function(
    doses,
    schedules,
    n_pat,
    criterion = "D",
    random_order = TRUE,
    interaction = TRUE
    ) {

  full.factorial <- expand.grid(doses = doses, schedules = schedules)
  full.factorial$linear_interaction <- doses * schedules
  full.factorial$linear_total <- NULL
  if(interaction) {
    full.factorial$linear_total <- ifelse(doses == 0 | schedules == 0, 0,
                                          doses + schedules)
  } else {
    full.factorial$linear_total <- ifelse(doses == 0 | schedules == 0, 0,
                                        doses + schedules + full.factorial$linear_interaction)
  }
  full.factorial_candidates <- do.call(bind_rows, replicate(n_pat, full.factorial, simplify = FALSE))

  d.optimal <- AlgDesign::optFederov(~ linear_total,
                          data = full.factorial_candidates, nTrials = n_pat, criterion = criterion)
  pat_alloc_freq <- table(d.optimal$design$doses, d.optimal$design$schedules)

  # Convert the table into a matrix
  m <- as.matrix(pat_alloc_freq)
  allocation_matrix <-  matrix(0, nrow = length(doses), ncol = length(schedules))
  rownames(allocation_matrix) <- doses
  colnames(allocation_matrix) <- schedules

  common_rows <- intersect(rownames(m), rownames(allocation_matrix))
  common_cols <- intersect(colnames(m), colnames(allocation_matrix))

  for (row in common_rows) {
    for (col in common_cols) {
      allocation_matrix[row, col] <- as.numeric(m[row, col])
    }
  }

  return(
    list(
      design = d.optimal,
      pat_alloc = pat_alloc_freq,
      mat = allocation_matrix
    )
  )
}


#' @title generate_optimal_design_emax
#'
#' @param doses  doses
#' @param schedules  schedules
#' @param n_pat  sample size to be allocated across the factorial space
#' @param criterion  optimality criteria, default "D"
#' @param random_order  boolean, default TRUE
#' @param emax_d  EMax of the doses
#' @param ed50_d  ED50 of the doses
#' @param emax_s  EMax of the schedules
#' @param ed50_s  ED50 of the schedules
#' @param interaction  boolean, default TRUE
#'
#' @return matrix
#' @export
generate_optimal_design_emax <- function(
    doses,
    schedules,
    n_pat,
    criterion = "D",
    random_order = TRUE,
    emax_d,
    ed50_d,
    emax_s,
    ed50_s,
    interaction = TRUE
) {
  full.factorial <- expand.grid(doses = doses, schedules = schedules)
  full.factorial$emax_doses <- full.factorial$doses*emax_d / (full.factorial$doses + ed50_d)
  full.factorial$emax_schedules <- full.factorial$schedules*emax_s / (full.factorial$schedules + ed50_s)
  full.factorial$emax_interaction <- full.factorial$emax_doses * full.factorial$emax_schedules
  full.factorial$emax_total <- NULL
  if(interaction) {
    full.factorial$emax_total <- ifelse(full.factorial$emax_doses == 0 | full.factorial$emax_schedules == 0, 0,
                                        full.factorial$emax_doses + full.factorial$emax_schedules)
  } else {
    full.factorial$emax_total <- ifelse(full.factorial$emax_doses == 0 | full.factorial$emax_schedules == 0, 0,
                                        full.factorial$emax_doses + full.factorial$emax_schedules + full.factorial$emax_interaction)
  }
  full.factorial_candidates <- do.call(bind_rows, replicate(n_pat, full.factorial, simplify = FALSE))

  d.optimal <- optFederov(~ emax_total,
                          data = full.factorial_candidates, nTrials = n_pat, criterion = criterion)
  pat_alloc_freq <- table(d.optimal$design$doses, d.optimal$design$schedules)

  m <- as.matrix(pat_alloc_freq)
  allocation_matrix <-  matrix(0, nrow = length(doses), ncol = length(schedules))
  rownames(allocation_matrix) <- doses
  colnames(allocation_matrix) <- schedules

  common_rows <- intersect(rownames(m), rownames(allocation_matrix))
  common_cols <- intersect(colnames(m), colnames(allocation_matrix))

  for (row in common_rows) {
    for (col in common_cols) {
      allocation_matrix[row, col] <- as.numeric(m[row, col])
    }
  }

  return(
    list(
      design = d.optimal,
      pat_alloc = pat_alloc_freq,
      mat = allocation_matrix
    )
  )
}


#' @title generate_custom_corner_mid
#'
#' @param doses  dose levels
#' @param schedules  schedule levels
#' @param n_pat  sample size to be allocated across the factorial space
#'
#' @return matrix
#' @export
generate_custom_corner_mid <- function(
    doses,
    schedules,
    n_pat
  ) {
  mat <- matrix(0, nrow = length(doses), ncol = length(schedules))
  rownames(mat) <- doses
  colnames(mat) <- schedules

  ncol <- length(schedules)
  nrow <- length(doses)

  size <- floor(n_pat / 5)

  remainder <- n_pat - size * 5

  # Assign a different value to the four corners
  mat[1, 1] <- size
  mat[1, ncol] <- size
  mat[nrow, 1] <- size
  mat[nrow, ncol] <- size

  # If the number of rows/columns is even, there will be 4 middle cells
  if (nrow %% 2 == 0 && ncol %% 2 == 0) {
    mid_row <- nrow / 2
    mid_col <- ncol / 2
    idx <- sample(1:4, 1)
    switch(idx,
           '1' = mat[mid_row, mid_col] <- size,
           '2' = mat[mid_row + 1, mid_col] <- size,
           '3' = mat[mid_row, mid_col + 1] <- size,
           '4' = mat[mid_row + 1, mid_col + 1] <- size,
    )
  } else {
    mid_row <- ceiling(nrow / 2)
    mid_col <- ceiling(ncol / 2)
    mat[mid_row, mid_col] <- size
  }

  non_zero_cells <- which(mat != 0)

  while(remainder > 0) {
    index <- sample(non_zero_cells, 1)
    mat[index] <- mat[index] + 1
    remainder <- remainder - 1
  }

  return(mat)
}

#' @title generate_custom
#'
#' @param doses  dose levels
#' @param schedules  schedule levels
#' @param n_pat  sample size to be allocated across the factorial space
#'
#' @return matrix
#' @export
generate_custom <- function(
    doses,
    schedules,
    n_pat
) {
  mat <- matrix(0, nrow = length(doses), ncol = length(schedules))
  rownames(mat) <- doses
  colnames(mat) <- schedules

  ncol <- length(schedules)
  nrow <- length(doses)

  size <- floor(n_pat / 5)

  remainder <- n_pat - size * 5

  # Assign a different value to the four corners
  mat[1, 1] <- size
  mat[1, ncol] <- size
  mat[nrow, 1] <- size
  mat[nrow, ncol] <- size

  # If the number of rows/columns is even, there will be 4 middle cells
  if (nrow %% 2 == 0 && ncol %% 2 == 0) {
    mid_row <- nrow / 2
    mid_col <- ncol / 2
    idx <- sample(1:4, 1)
    switch(idx,
           '1' = mat[mid_row, mid_col] <- size,
           '2' = mat[mid_row + 1, mid_col] <- size,
           '3' = mat[mid_row, mid_col + 1] <- size,
           '4' = mat[mid_row + 1, mid_col + 1] <- size,
    )
  } else {
    mid_row <- ceiling(nrow / 2)
    mid_col <- ceiling(ncol / 2)
    mat[mid_row, mid_col] <- size
  }

  non_zero_cells <- which(mat != 0)

  while(remainder > 0) {
    index <- sample(non_zero_cells, 1)
    mat[index] <- mat[index] + 1
    remainder <- remainder - 1
  }

  return(mat)
}

#' @title generate_custom_design
#'
#' @param n_pat  number of patients allocated across the factorial space
#' @param post_probs  posterior probabilities of the grid the patients should be allocated on
#'
#' @return matrix
#' @export
generate_custom_design <- function(
    n_pat,
    post_probs
) {
  probabilities <- as.vector(post_probs)
  probabilities <- probabilities / sum(probabilities)


  # Initialize the allocation matrix
  allocation_grid <- matrix(0, nrow = nrow(post_probs), ncol = ncol(post_probs))

  # Nested loop to assign each cell with the probabilities times 200
  for (i in 1:nrow(post_probs)) {
    for (j in 1:ncol(post_probs)) {
      allocation_grid[i, j] <- floor(probabilities[(i - 1) * ncol(post_probs) + j] * n_pat)
    }
  }
  # Ensure there are no NA values in allocation_grid and post_probs
  allocation_grid[is.na(allocation_grid)] <- 0
  post_probs[is.na(post_probs)] <- 0

  # Calculate the difference
  difference <- n_pat - sum(allocation_grid, na.rm = TRUE)


  if(difference > 0) {
    for(i in 1:difference){
      x <- sample(1:nrow(post_probs), 1)
      y <- sample(1:ncol(post_probs), 1)
      allocation_grid[x, y] <- allocation_grid[x, y] + 1
    }
  }


  # Print the allocation grid
  # print(allocation_grid)


  return(allocation_grid)
}

