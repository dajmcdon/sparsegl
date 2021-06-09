#' Plot solution paths from a `sparsegl` object
#' 
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{\link{sparsegl}} object.
#' 
#' A coefficient profile plot is produced.
#' 
#' @param x fitted \code{\link{sparsegl}} model
#' @param grouped if applying the user-specified grouping of features. Default
#' is FALSE.
#' @param y_axis what is on the Y_axis. Plot \code{"coef"} the coefficients or 
#' \code{"norm"} group norm. Default is \code{"coef"}.
#' @param x_axis what is on the X-axis. Plot against the \code{"log_lambda"} log-lambda
#' sequence or \code{"norm"} a scaled norm vector. Default is \code{"log_lambda"}.
#' @param asparse the weight to put on the ell1 norm in sparse group lasso. Default
#' is 0.05
#' @param \dots other graphical parameters to plot
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression
#' @method plot sparsegl
#' @export
plot.sparsegl <- function(x, grouped = FALSE, asparse = 0.05, 
                          y_axis = c("coef", "group_norm"),
                          x_axis = c("log_lambda", "norm"), ...) {
    
    y_axis <- match.arg(y_axis)
    x_axis <- match.arg(x_axis)
    
    xb <- x$beta  
    nonzeros <- sort(unique(xb@i)) + 1
    assertthat::assert_that(length(nonzeros) > 0, 
                            msg = "No nonzero betas / groups are available to plot")
    
    tmp <- xb[nonzeros, , drop = FALSE]  
    g <- as.numeric(x$group[nonzeros])
    uni_group <- unique(g)
    n.g <- length(uni_group)
    sgnorm <- apply(xb, 2, function(y) sp_group_norm(y, x$group))
    
    if (grouped && y_axis == "group_norm") {
        beta <- matrix(NA, n.g, length(x$lambda))  
        j <- 1
        for (i in uni_group) {
            beta[j, ] <- apply(tmp[g == i, ], 2, function(x) {
                asparse * sum(abs(x)) + (1 - asparse) * sqrt(sum(x^2))})
            j <- j + 1
            }
        } else {
            beta <- tmp
    }
    
    transform_outputs <- function(df) {
        if (x_axis == "norm") {
            df <- df %>% dplyr::mutate(norm = sgnorm / max(sgnorm))
        } else {
            df <- df %>% dplyr::mutate(lambda = x$lambda) 
            }
        
        if (grouped && y_axis == "group_norm") {
            if (x_axis == "norm") {
                df <- df %>% 
                    tidyr::pivot_longer(!.data$norm, names_to = "variable")
            } else {
            df <- df %>% 
                tidyr::pivot_longer(!.data$lambda, names_to = "variable")
            }
            
        } else {
            if (x_axis == "norm") {
                df <- df %>% 
                    tidyr::pivot_longer(!c(.data$norm, .data$group), names_to = "variable")
            } else {
            df <- df %>% 
                tidyr::pivot_longer(!c(.data$lambda, .data$group), names_to = "variable")
            }
        }
        
        df <- df %>% dplyr::mutate(variable = factor(
            .data$variable, levels = sort(as.numeric(unique(.data$variable)))))
        
        return(df)
    }
    
    
    outputs <- tibble::as_tibble(t(as.matrix(beta)))
    if (grouped && y_axis == "group_norm")  {
        colnames(outputs) <- uni_group 
        outputs <- transform_outputs(outputs)
    } else {
        colnames(outputs) <- nonzeros
        outputs <- transform_outputs(outputs %>% dplyr::mutate(group = 0))
        
        j <- 1
        for (i in uni_group) {
            outputs <- outputs %>%
                dplyr::mutate(group = replace(.data$group, .data$variable %in%
                                                  nonzeros[g == i], uni_group[j]))
            j <- j + 1
        }
        outputs <- outputs %>% dplyr::mutate(group = factor(.data$group))
    }
    
    if (!grouped || (grouped && y_axis == "group_norm"))  {
        if (x_axis == "norm") {
            plot_layer <- ggplot2::ggplot(outputs, ggplot2::aes(x = .data$norm, 
                                                                y = .data$value, 
                                                                color = .data$variable))
        } else {
            plot_layer <- ggplot2::ggplot(outputs, ggplot2::aes(x = .data$lambda, 
                                                                y = .data$value, 
                                                                color = .data$variable))
        }
    } else {
        if (x_axis == "norm") {
            plot_layer <- ggplot2::ggplot(outputs, ggplot2::aes(x = .data$norm,
                                                                y = .data$value,
                                                                group = .data$variable,
                                                                color = .data$group))
        } else {
            plot_layer <- ggplot2::ggplot(outputs, ggplot2::aes(x = .data$lambda,
                                                                y = .data$value,
                                                                group = .data$variable,
                                                                color = .data$group))
        }
    }
        
    plot_layer <- plot_layer +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0)
    
    if (x_axis == "norm") {
        xlab_layer <- ggplot2::xlab("beta norm / max (beta norm)")
    } else {
        xlab_layer <- ggplot2::xlab("Log Lambda") + ggplot2::scale_x_log10() 
    }
    
    if (grouped && y_axis == "group_norm") {
        ylab_layer <- ggplot2::ylab("Group norm")
    } else {
        ylab_layer <- ggplot2::ylab("Coefficients")
    }
    
    if (grouped) {
        legend_layer <- ggplot2::scale_color_discrete(
            name = "Group", 
            labels = paste0("group", uni_group))
    } else {
        legend_layer <- ggplot2::scale_color_discrete(
            name = "Variable", 
            labels = paste0("variable", nonzeros))
    }
    
    theme_layer <- ggplot2::theme_bw()
    
    p <- plot_layer +
        xlab_layer +
        ylab_layer +
        legend_layer +
        theme_layer
    
    return(p)
    
}

