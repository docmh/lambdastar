library(lambdastar)

# Deliberately base R: the tests need no optional testing installation.
.count <- 0L
.failures <- character()
test <- function(name, code) {
  tryCatch({ force(code); .count <<- .count + 1L; cat("PASS:", name, "\n") },
    error = function(e) {
      .failures <<- c(.failures, paste(name, conditionMessage(e), sep = ": "))
      cat("FAIL:", tail(.failures, 1), "\n")
    })
}
equal <- function(actual, expected, tol = 1e-10) {
  result <- all.equal(actual, expected, tolerance = tol, check.attributes = TRUE)
  if (!isTRUE(result)) stop(paste(result, collapse = "; "))
}
relative <- function(actual, expected, tol = 1e-8) {
  stopifnot(length(actual) == length(expected), all(is.finite(actual)),
            all(is.finite(expected)), all(expected != 0),
            all(abs(actual / expected - 1) <= tol))
}
fails <- function(code, pattern) {
  err <- tryCatch({ force(code); NULL }, error = identity)
  if (is.null(err) || !grepl(pattern, conditionMessage(err)))
    stop("Expected error matching: ", pattern)
}
test_root <- if (dir.exists("cases")) "." else "tests"
fixture <- function(name) file.path(test_root, "fixtures", name)
for (file in list.files(file.path(test_root, "cases"), full.names = TRUE)) source(file)
cat(sprintf("\nLambda Star: %d passed, %d failed.\n", .count, length(.failures)))
if (length(.failures)) stop(paste(.failures, collapse = "\n"))
