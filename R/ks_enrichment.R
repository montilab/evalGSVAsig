## NOTE: copied from hypeR/R/ks_enrichment.R
##
#' An empty ggplot
#'
#' @return A ggplot object
#'
#' @importFrom ggplot2 qplot geom_rug ggplot geom_hline geom_vline annotate theme theme_void element_text element_line element_rect element_blank
#'
.ggempty <- function() {
  ggplot2::ggplot() +
    ggplot2::theme_void()
}
#' Enrichment plot implemented in ggplot
#'
#' @param n The length of a ranked list
#' @param positions A vector of positions in the ranked list
#' @param x_axis The x-axis of a running enrichment score
#' @param y_axis The y-axis of a running enrichment score
#' @param title Plot title
#' @return A ggplot object
#'
#' @importFrom ComplexHeatmap anno_lines
#'
.ggeplot <- function(n, positions, x_axis, y_axis, title = "") {
  score <- which.max(abs(y_axis))
  ggplot2::qplot(x_axis,
    y_axis,
    main = title,
    ylab = "Running Enrichment Score",
    xlab = "Position in Ranked List of Genes",
    geom = "line"
  ) +
    ggplot2::geom_rug(data = data.frame(positions), aes(x = positions), inherit.aes = FALSE) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = n / 2, linetype = "dotted") +
    ggplot2::annotate("point", x = x_axis[score], y = y_axis[score], color = "red") +
    ggplot2::annotate("text", x = x_axis[score] + n / 20, y = y_axis[score], label = round(y_axis[score], 2)) +
    ggplot2::annotate("point", x = x_axis[score], y = y_axis[score], color = "red") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, size = 1)
    )
}
## this needs to extrapolate so that x and y are same lenght as n
.ks_anno_lines <- function( n, positions, x_axis, y_axis, title="" )
{
  score <- which.max(abs(y_axis))
  ComplexHeatmap::anno_lines(
    x = data.frame(x = x_axis, y = y_axis) |> as.matrix(),
    which = "row", smooth = TRUE)
}
#' One-sided Kolmogorov–Smirnov test
#'
#' @param n.x The length of a ranked list
#' @param y A vector of positions in the ranked list
#' @param weights Weights for weighted score (Subramanian et al.)
#' @param weights.pwr Exponent for weights (Subramanian et al.)
#' @param absolute Takes max-min score rather than the max deviation from null
#' @param plotting Use true to generate plot
#' @param plot.title Plot title
#' @return A list of data and plots
#'
#' @importFrom stats ks.test
#'
#' @keywords internal
.kstest <- function(n.x,
                    y,
                    weights = NULL,
                    weights.pwr = 1,
                    absolute = FALSE, # this is not really implemented, should be removed
                    plotting = FALSE,
                    plot.title = "") {
  n.y <- length(y)
  err <- list(score = 0, pval = 1, leading_edge = NULL, leading_hits = NULL, plot = .ggempty())

  if ( n.y < 1 || any(y > n.x) || any(y < 1) ) {
    return(err)
  }
  x.axis <- y.axis <- NULL
  leading_edge <- NULL # recording the x corresponding to the highest y value

  ## If weights are provided
  if ( !is.null(weights) ) {
    weights <- abs(weights[y])^weights.pwr

    Pmis <- rep(1, n.x)
    Pmis[y] <- 0
    Pmis <- cumsum(Pmis)
    Pmis <- Pmis / (n.x - n.y)
    Phit <- rep(0, n.x)
    Phit[y] <- weights
    Phit <- cumsum(Phit)
    Phit <- Phit / Phit[n.x]
    z <- Phit - Pmis

    score <- if (absolute) max(z) - min(z) else z[leading_edge <- which.max(abs(z))]

    x.axis <- 1:n.x
    y.axis <- z
  }
  ## Without weights
  else {
    y <- sort(y)
    n <- n.x * n.y / (n.x + n.y)
    hit <- 1 / n.y
    mis <- 1 / n.x

    Y <- sort(c(y - 1, y)) # append the positions preceding hits
    Y <- Y[diff(Y) != 0] # remove repeated position
    y.match <- match(y, Y) # find the hits' positions
    D <- rep(0, length(Y))
    D[y.match] <- (1:n.y)
    zero <- which(D == 0)[-1]
    D[zero] <- D[zero - 1]

    z <- D * hit - Y * mis

    score <- if (absolute) max(z) - min(z) else z[leading_edge <- which.max(abs(z))]

    x.axis <- Y
    y.axis <- z

    if (Y[1] > 0) {
      x.axis <- c(0, x.axis)
      y.axis <- c(0, y.axis)
    }
    if (max(Y) < n.x) {
      x.axis <- c(x.axis, n.x)
      y.axis <- c(y.axis, 0)
    }
  }
  leading_edge <- x.axis[leading_edge]

  # One-sided Kolmogorov–Smirnov test
  results <- suppressWarnings(ks.test(1:n.x, y, alternative = "less"))
  results$statistic <- score # Use the signed statistic

  # Enrichment plot
  p <- if (plotting) {
    .ggeplot(n.x, y, x.axis, y.axis, plot.title) +
      geom_vline(xintercept = leading_edge, linetype = "dotted", color = "red", size = 0.25)
  } else {
    .ggempty()
  }
  ##p <- .ks_anno_lines(n.x, y, x.axis, y.axis)
  return(list(
    score = as.numeric(results$statistic),
    pval = results$p.value,
    leading_edge = leading_edge,
    x_axis = x.axis,
    y_axis = y.axis,
    plot = p
  ))
}

