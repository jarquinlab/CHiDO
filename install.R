packages <- readLines("requirements.txt")

installed <- character()
failed <- character()

for (pkg in packages) {
  if (require(pkg, character.only = TRUE, quietly = TRUE)) {
    next
  }
  tryCatch(
    {
      install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE, quiet = TRUE)
      if (require(pkg, character.only = TRUE, quietly = TRUE)) {
        installed <- c(installed, pkg)
      } else {
        failed <- c(failed, pkg)
      }
    },
    error = function(e) {
      failed <<- c(failed, pkg)
    }
  )
}

if (length(installed) > 0) cat("Installed:", paste(installed, collapse = ", "), "\n")
if (length(failed)    > 0) cat("Failed:   ", paste(failed,    collapse = ", "), "\n")
if (length(failed) == 0 && length(installed) == 0) cat("All packages already installed.\n")
