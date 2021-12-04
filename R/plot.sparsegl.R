#' Plot solution paths from a `sparsegl` object.
#'
#' Produces a coefficient profile plot of a fitted
#' [sparsegl()] object. The result is a [ggplot2::ggplot()]. Additional user
#' modifications can be added as desired.
#'
#' @param x Fitted [sparsegl()] object.
#' @param y_axis Variable on the y_axis. Either the coefficients (default)
#'   or the group norm.
#' @param x_axis Variable on the x-axis. Either the (log)-lambda
#'   sequence (default) or the value of the penalty. In the second case,
#'   the penalty is scaled by its maximum along the path.
#' @param add_legend Show the legend. Often, with many groups/predictors, this
#'   can become overwhelming.
#' @param ... Not used.
#' @seealso [sparsegl()].
#' @method plot sparsegl
#' @export
#' @examples
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' beta_star <- c(rep(5, 5), c(5, -5, 2, 0, 0), rep(-5, 5), rep(0, (p - 15)))
#' y <- X %*% beta_star + eps
#' groups <- rep(1:(p / 5), each = 5)
#' fit1 <- sparsegl(X, y, group = groups)
#' plot(fit1, y_axis = "coef", x_axis = "penalty")
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
        ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(x_axis),
                                     y = .data$value,
                                     color = !!rlang::sym(y_axis))) +
        ggplot2::geom_hline(yintercept = 0)

    if (x_axis == "penalty") xlab_layer <- ggplot2::xlab("penalty / max (penalty)")
    else xlab_layer <- ggplot2::xlab("lambda") + ggplot2::scale_x_log10()

    if (y_axis == "group") {
        plot_layer <- plot_layer +
            ggplot2::geom_line(ggplot2::aes(group = .data$group)) +
            ggplot2::ylab("group norm")
    } else {
        plot_layer <- plot_layer +
            ggplot2::geom_line() +
            ggplot2::ylab("coefficients")
    }

    legend_layer <- ggplot2::scale_color_viridis_d(name = "")

    theme_layer <- ggplot2::theme_bw()
    if (!add_legend)
        theme_layer <- theme_layer + ggplot2::theme(legend.position = "none")

    p <- plot_layer + xlab_layer + legend_layer + theme_layer
    return(p)
}

