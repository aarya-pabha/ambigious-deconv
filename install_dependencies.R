#!/usr/bin/env Rscript
required_packages <- c("quadprog", "MASS", "Matrix", "clue", "lpSolve")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) {
  install.packages(new_packages, repos = "https://cloud.r-project.org")
  cat("All dependencies successfully installed.\n")
} else {
  cat("All dependencies are already installed.\n")
}
