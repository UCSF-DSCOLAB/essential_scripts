#' Print a timestamped log message
#' @param ... any number of strings (or other types that can be converted to a string)
#' @return None
#' @details This function prints a timestamped (ts) log message where all elements
#' to the function are concatenated with no spacer and printed after a timestamp
#' in the form of: \code{[2025-03-11T14:05:47] <msg>}
#' @examples
#' ts_log("work starting")
#' 
#' # Also works easily with variable data too:
#' x <- "XXX-HS25"
#' ts_log("Working on sample: ", x)

ts_log <- function(...) {
    cat("[", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), "] ", ..., "\n", sep = "")
}
