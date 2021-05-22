utils::globalVariables(c("starts_with", "lambda", "value","variable"))
#' Plot solution paths from a "sparsegl" object
#' 
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{\link{sparsegl}} object.
#' 
#' A coefficient profile plot is produced.
#' 
#' @param x fitted \code{\link{sparsegl}} model
#' @param group what is on the Y-axis. Plot the norm of each group if
#' \code{TRUE}. Plot each coefficient if \code{FALSE}.
#' @param log.l what is on the X-axis. Plot against the log-lambda sequence if
#' \code{TRUE}. Plot against the lambda sequence if \code{FALSE}.
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
plot.sparsegl <- function(x, group = TRUE, log.l = TRUE, asparse = 0.05, ...) {
    
    xb <- x$beta  
    
    if (nrow(xb) == 1) {
        if (any(abs(xb) > 0)) {
            nonzeros <- 1                             
        } else nonzeros <- NULL
    } else {
        nonzeros <- which(apply(abs(xb), 1, sum) > 0)
    }
    
    tmp <- xb[nonzeros, , drop = FALSE]  
    g <- as.numeric(as.factor(x$group[nonzeros])) 
    l <- x$lambda  
    n.g <- max(g) 
    
    if (group) {
        bs <- as.integer(as.numeric(table(g)))  
        ix <- rep(NA, n.g)
        iy <- rep(NA, n.g)
        j <- 1
        for (i in 1:n.g) {
            ix[i] <- j
            iy[i] <- j + bs[i] - 1
            j <- j + bs[i]
        }
        beta <- matrix(NA, n.g, length(l))  
        for (i in 1:n.g) {
            crossp <- apply(tmp[ix[i]:iy[i], ], 2, function(x) {asparse * sum(abs(x)) + (1 - asparse) * sqrt(crossprod(x))})
            beta[i,] <- crossp
        }
    } else beta <- tmp

    if (log.l) l <- log(l)

    outputs <- as.data.frame(t(as.matrix(beta)))
    colnames(outputs) <- c(1:n.g)
    outputs <- outputs %>% dplyr::mutate(lambda = l)
    outputs1 <- as.data.frame(outputs)
    outputs1 <- outputs1 %>% tidyr::pivot_longer(cols = !starts_with("lambda"), names_to = "variable")
    
    # ggplot
    p1 <- outputs1 %>%
        ggplot2::ggplot(ggplot2::aes(x = lambda, y = value, col = factor(variable))) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) + 
        ggplot2::ylab("Coefficients")

        
    if (group) {
        p1 <- p1 + ggplot2::scale_color_discrete(name = "Group", 
                                                 labels = paste0("group", 1:n.g), 
                                                 breaks = sort(as.numeric(levels(factor(outputs1$variable))))); 
    } else {
        p1 <- p1 + ggplot2::scale_color_discrete(name = "Variable", 
                                                 labels = paste0("variable", 1:n.g),
                                                 breaks = sort(as.numeric(levels(factor(outputs1$variable)))))
        }

    if (log.l) p1 <- p1 + ggplot2::xlab("Log Lambda") else p1 <- p1 + ggplot2::xlab("Lambda")

    # standardized group norm on x-axis
    sgnorm <- apply(xb, 2, function(y) sp_group_norm(y, x$group))
    outputs2 <- outputs
    outputs2 <- outputs2 %>% 
        dplyr::mutate(lambda = sgnorm / max(sgnorm)) %>% 
        tidyr::pivot_longer(!starts_with("lambda"), names_to = "variable")

    p2 <- outputs2 %>%
        ggplot2::ggplot(ggplot2::aes(x = lambda, y = value, col = factor(variable))) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::ylab("Coefficients") +
        ggplot2::xlab("Standardized Lambda") +
        ggplot2::theme(legend.position = "none")
    
    outputs <- tibble::as_tibble(t(as.matrix(tmp))) 
    names <- c(1:dim(tmp)[1])
    colnames(outputs) <- names
    outputs <- outputs %>% dplyr::mutate(lambda = l)
    outputs3 <- outputs %>% tidyr::pivot_longer(cols = !starts_with("lambda"), names_to = "variable")
    outputs4 <- outputs %>% dplyr::mutate(lambda = sgnorm / max(sgnorm)) %>% tidyr::pivot_longer(cols = !starts_with("lambda"), names_to = "variable")
    if (group) {
        outputs3 <- outputs3 %>% dplyr::mutate(groups = 0)
        outputs4 <- outputs4 %>% dplyr::mutate(groups = 0)
        for (i in 1:n.g) {
           outputs3 <- outputs3 %>% dplyr::mutate(groups = replace(groups, variable %in% names[ix[i]:iy[i]], i))
           outputs4 <- outputs4 %>% dplyr::mutate(groups = replace(groups, variable %in% names[ix[i]:iy[i]], i))
        }
        p3 <- outputs3 %>% 
            ggplot2::ggplot(ggplot2::aes(x = lambda, y = value, group = variable, color = factor(groups))) +
            ggplot2::geom_line() +
            ggplot2::geom_hline(yintercept = 0) +
            ggplot2::ylab("Coefficients") +
            ggplot2::scale_color_discrete(name = "Group", labels = paste0("group", 1:n.g)) 
        
        p4 <- outputs4 %>% 
            ggplot2::ggplot(ggplot2::aes(x = lambda, y = value, group = variable, col = factor(groups))) +
            ggplot2::geom_line() +
            ggplot2::geom_hline(yintercept = 0) +
            ggplot2::scale_color_discrete(name = "Group", labels = paste0("group", 1:n.g)) +
            ggplot2::xlab("Standardized Lambda") +
            ggplot2::ylab("Coefficients")
            
    }
    
    if (group) {
        p_left <- ggpubr::ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "left") 
        p_right <- ggpubr::ggarrange(p3, p4, nrow = 2, common.legend = TRUE, legend = "right")
        p <- ggpubr::ggarrange(p_left, p_right, nrow = 1, ncol = 2)
        return(p)
    } else {
        p <- ggpubr::ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "right")
        return(p)
    }
}
    

