compress_paths_multi <- function(paths) {
  
  # Sort once
  paths <- paths[order(paths$id, paths$time), ]
  
  # Extract state number
  paths$state <- as.integer(
    sub(".*?(\\d+).*", "\\1", paths$state)
  )
  
  # Extract cluster number
  paths$cluster <- as.integer(
    sub(".*?(\\d+).*", "\\1", paths$cluster)
  )
  
  n <- nrow(paths)
  
  # Identify consecutive duplicates within subject
  same_as_prev <- c(
    FALSE,
    paths$id[-1] == paths$id[-n] &
      paths$state[-1] == paths$state[-n]
  )
  
  # Keep rows that are NOT duplicates
  keep <- !same_as_prev
  
  # Always keep last row of each subject
  last_rows <- c(
    which(paths$id[-1] != paths$id[-n]),
    n
  )

  keep[last_rows] <- TRUE
  
  compressed <- paths[keep, , drop = FALSE]
  
  rownames(compressed) <- NULL
  
  compressed
}