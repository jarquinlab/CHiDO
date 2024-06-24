# List all R script files in the specified directory
list_r_files <- function() {
  list.files("code", full.names = TRUE, recursive = TRUE)
}

# Extract library calls from a file
extract_libraries <- function(file) {
  lines <- readLines(file, warn = FALSE)
  libraries <- unique(c(
    unlist(regmatches(lines, gregexpr("library\\(([^)]+)\\)", lines))),
    unlist(regmatches(lines, gregexpr("require\\(([^)]+)\\)", lines))),
    unlist(regmatches(lines, gregexpr("[a-zA-Z0-9.]+::", lines)))
  ))
  libraries <- gsub("library\\(|require\\(|\\)", "", libraries)
  libraries <- gsub("::.*", "", libraries)
  return(libraries)
}

# Get unique libraries used in all files
get_all_libraries <- function(directory = "R") {
  r_files <- c(list_r_files(), "app.R")
  all_libraries <- unlist(lapply(r_files, extract_libraries))
  unique(all_libraries)
}

# Print the list of libraries
print(get_all_libraries())
