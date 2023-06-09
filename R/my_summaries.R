
#' Summary for a numeric vector
#'
#' \code{numeric_summary} calculates extended summary for a numeric vector
#'
#' @param x numeric vector
#' @param na.rm whether to include NA values in the calculations
#'
#' @return a named numeric vector with 6 values
#' @export
#' @import dplyr
#'
#' @examples
#'
#' x <- rnorm(100)
#' numeric_summary(x)
#'
numeric_summary <- function(x, na.rm=FALSE){

  if (!is.null(x) && !is.numeric(x))
    stop("'x' must be a numeric vector")

  min = min(x, na.rm=na.rm)
  max = max(x, na.rm=na.rm)
  mean = mean(x, na.rm=na.rm)
  sd = sd(x, na.rm=na.rm)
  length = length(x)
  Nmiss = sum(is.na(x))

  c(min=min, max=max, mean=mean, sd=sd, length=length, Nmiss=Nmiss)

}

# This function creates a summary for a character vector
#' Summary for character vector
#'
#' Description blah blah
#'
#' @param x character vector
#' @param na.rm whether to remove NAs
#'
#' @return a named list of 3
#' @export
#'
#' @examples
#'
#' x <- c("Boston", "NY", "Boston", "New Haven")
#' char_summary(x)
#'
char_summary <- function(x, na.rm=FALSE){

  if (!is.null(x) && !is.character(x))
    stop("'x' must be a character vector")

  length = length(x)
  Nmiss = sum(is.na(x))
  Nunique = length(unique(x))

  c(length = length,
    Nmiss = Nmiss,
    Nunique = Nunique )

}
