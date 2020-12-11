#' Generate a random alphanumeric string
#'
#' @param length Length of string
#'
#' @return A character string
#' @export
#'
#' @examples
#' random_alphanumeric(20)
random_alphanumeric <- function(length) {
  paste(sample(c(LETTERS, letters, 0:9),
               size = length,
               replace = TRUE),
        collapse = "")
}
