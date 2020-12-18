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

#' Median
#'
#' See \code{stats::\link[stats:median]{median}} for details.
#'
#' @name median
#' @rdname median
#' @keywords internal
#' @export
#' @importFrom stats median
NULL

#' Data
#'
#' See \code{utils::\link[utils:data]{data}} for details.
#'
#' @name data
#' @rdname data
#' @keywords internal
#' @export
#' @importFrom utils data
NULL

#' head
#'
#' See \code{utils::\link[utils:head]{head}} for details.
#'
#' @name head
#' @rdname head
#' @keywords internal
#' @export
#' @importFrom utils head
NULL

#' read.delim
#'
#' See \code{utils::\link[utils:read.table]{read.table}} for details.
#'
#' @name read.delim
#' @rdname read.delim
#' @keywords internal
#' @export
#' @importFrom utils read.delim
NULL

#' stack
#'
#' See \code{utils::\link[utils:stack]{stack}} for details.
#'
#' @name stack
#' @rdname stack
#' @keywords internal
#' @export
#' @importFrom utils stack
NULL
