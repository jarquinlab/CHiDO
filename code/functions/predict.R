extract_omics <- function(formula_str) {
  clean_formula <- gsub("\\s+", "", formula_str)
  terms <- unlist(strsplit(clean_formula, "\\+"))
  vars <- unique(unlist(strsplit(terms, "\\*|\\+")))
  return(setdiff(vars, c("E", "L")))
}
