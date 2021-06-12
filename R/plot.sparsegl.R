#' Plot solution paths from a `sparsegl` object
#'
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{\link{sparsegl}} object. The result is a `ggplot`. Additional user
#' modifications can be added as desired.
#'
#' @param x fitted \code{\link{sparsegl}} model
#' @param y_axis variable is on the y_axis. Either the coefficients (default)
#'   or the group norm.
#' @param x_axis variable on the x-axis. Either the (log)-lambda
#'   sequence (default) or value of the penalty. The penalty is scaled by its
#'   maximum along the path.
#' @param add_legend show the legend. Often, with many groups/predictors, this
#'   can become overwhelming.
#' @param \dots other graphical parameters to plot
#' @method plot sparsegl
#' @export
plot.sparsegl <- function(x,
                          y_axis = c("coef", "group"),
                          x_axis = c("lambda", "penalty"),
                          add_legend = TRUE,
                          ...) {

    y_axis <- match.arg(y_axis)
    x_axis <- match.arg(x_axis)

    xb <- x$beta
    nonzeros <- sort(unique(xb@i)) + 1
    assertthat::assert_that(
        length(nonzeros) > 0,
        msg = "No nonzero betas / groups are available to plot")

    xb <- xb[nonzeros, , drop = FALSE]
    g <- x$group[nonzeros]
    uni_group <- unique(g)
    n_groups <- length(uni_group)
    sgnorm <- apply(xb, 2, sp_group_norm, gr = g, asparse = x$asparse)

    if (y_axis == "group") {
        xb <- apply(xb, 2, grouped_sp_norm, gr = g, asparse = x$asparse)
        rownames(xb) <- paste0("group", uni_group)
    }
    df <- as.data.frame(t(as.matrix(xb)))
    df$lambda <- x$lambda
    df$penalty <- sgnorm / max(sgnorm)
    df <- df %>%
        tidyr::pivot_longer(!c(.data$lambda, .data$penalty), names_to = y_axis)

    plot_layer <- df %>%
        ggplot2::ggplot(
            ggplot2::aes(x = !!rlang::sym(x_axis),
                         y = .data$value,
                         color = !!rlang::sym(y_axis))) +
        ggplot2::geom_hline(yintercept = 0)

    if (x_axis == "norm") {
        xlab_layer <- ggplot2::xlab("penalty / max (penalty)")
    } else {
        xlab_layer <- ggplot2::xlab("lambda") + ggplot2::scale_x_log10()
    }

    if (y_axis == "group") {
        plot_layer <- plot_layer +
            ggplot2::geom_line(ggplot2::aes(group = .data$group)) +
            ggplot2::ylab("group norm")
    } else {
        plot_layer <- plot_layer +
            ggplot2::geom_line() +
            ggplot2::ylab("coefficients")
    }

    legend_layer <- ggplot2::scale_color_viridis_d()
    theme_layer <- ggplot2::theme_bw()
    if (!add_legend)
        theme_layer <- theme_layer + ggplot2::theme(legend.position = "none")

    p <- plot_layer + xlab_layer + legend_layer + theme_layer
    return(p)
}

