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
plot.sparsegl <- function(x, group = FALSE, log.l = TRUE, asparse = 0.05,
                          ...) {
    
    xb <- x$beta  
    
    if (nrow(xb) == 1) {
        if (any(abs(xb) > 0)) {
            nonzeros <- 1                             
        } else nonzeros <- NULL
    } else {
        nonzeros <- which(apply(abs(xb), 1, sum) > 0)
    }
    
    tmp <- xb[nonzeros, , drop = FALSE]  
    g <- as.numeric(x$group[nonzeros])
    l <- x$lambda  
    uni_group <- unique(g)
    n.g <- length(uni_group)
    
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
            crossp <- apply(tmp[ix[i]:iy[i], ], 2, function(x) {
                asparse * sum(abs(x)) + (1 - asparse) * sqrt(crossprod(x))})
            beta[i,] <- crossp
        }
    } else beta <- tmp
    
    if (log.l) l <- log(l)
    
    outputs <- tibble::as_tibble(t(as.matrix(beta)), .name_repair = "unique")
    if (group) colnames(outputs) <- uni_group else colnames(outputs) <- c(1:dim(outputs)[2])
    outputs <- outputs %>% dplyr::mutate(lambda = l)
    outputs1 <- outputs %>%
        tidyr::pivot_longer(cols = !.data$lambda, names_to = "variable") %>% 
        dplyr::mutate(variable = factor(
            .data$variable, levels = sort(as.numeric(unique(.data$variable)))))
    
    # ggplot
    p1 <- outputs1 %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$lambda, 
                                     y = .data$value, 
                                     col = .data$variable)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) + 
        ggplot2::theme_bw()
    
    
    if (group) {
        p1 <- p1 + 
            ggplot2::ylab("Group norm") +
            ggplot2::scale_color_discrete(name = "Group", 
                                          labels = paste0("group", uni_group))
    } else {
        p1 <- p1 + 
            ggplot2::ylab("Coefficients") +
            ggplot2::scale_color_discrete(name = "Variable", 
                                          labels = paste0("variable", as.numeric(nonzeros)))
    }
    
    if (log.l) p1 <- p1 + ggplot2::xlab("Log Lambda") else p1 <- p1 + ggplot2::xlab("Lambda")
    
    # standardized group norm on x-axis
    sgnorm <- apply(xb, 2, function(y) sp_group_norm(y, x$group))
    outputs2 <- outputs
    outputs2 <- outputs2 %>% 
        dplyr::mutate(lambda = sgnorm / max(sgnorm)) %>% 
        tidyr::pivot_longer(!.data$lambda, names_to = "variable") %>% 
        dplyr::mutate(variable = factor(.data$variable, 
                                        levels = sort(as.numeric(unique(.data$variable)))))
    
    p2 <- outputs2 %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$lambda,
                                     y = .data$value, 
                                     col = .data$variable)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::xlab("beta norm / max(beta norm)") +
        ggplot2::theme(legend.position = "none") + 
        ggplot2::theme_bw()
    
    if (group) {
        p2 <- p2 + 
            ggplot2::ylab("Group norm") +
            ggplot2::scale_color_discrete(name = "Group", 
                                          labels = paste0("group", uni_group))
    } else {
        p2 <- p2 + 
            ggplot2::ylab("Coefficients") +
            ggplot2::scale_color_discrete(name = "Variable", 
                                          labels = paste0("variable", as.numeric(nonzeros)))
    }
    
    if (!group) {
        return(list(p1, p2))
    }
    
    outputs <- tibble::as_tibble(t(as.matrix(tmp)), .name_repair = "unique") 
    names <- c(1:dim(tmp)[1])
    colnames(outputs) <- names
    outputs3 <- outputs %>% 
        dplyr::mutate(lambda = l) %>% 
        dplyr::mutate(group = 0) %>% 
        tidyr::pivot_longer(!c(.data$lambda, .data$group), names_to = "variable") %>% 
        dplyr::mutate(variable = factor(
            .data$variable, levels = sort(as.numeric(unique(.data$variable)))))
    outputs4 <- outputs %>%
        dplyr::mutate(lambda = sgnorm / max(sgnorm)) %>%
        dplyr::mutate(group = 0) %>% 
        tidyr::pivot_longer(!c(.data$lambda, .data$group), names_to = "variable") %>% 
        dplyr::mutate(variable = factor(
            .data$variable, levels = sort(as.numeric(unique(.data$variable)))))
    
    j <- 1
    for (i in 1:n.g) {
        outputs3 <- outputs3 %>%
            dplyr::mutate(group = replace(group, .data$variable %in%
                                              names[ix[i]:iy[i]], uni_group[j]))
        outputs4 <- outputs4 %>%
            dplyr::mutate(group = replace(group, .data$variable %in%
                                              names[ix[i]:iy[i]], uni_group[j]))
        j <- j + 1
    }
    
    outputs3 <- outputs3 %>% dplyr::mutate(group = factor(group))
    outputs4 <- outputs4 %>% dplyr::mutate(group = factor(group))
    p3 <- outputs3 %>% 
        ggplot2::ggplot(ggplot2::aes(x = .data$lambda, 
                                     y = .data$value, 
                                     group = .data$variable, 
                                     color = .data$group)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::ylab("Coefficients") +
        ggplot2::scale_color_discrete(name = "Group", labels = paste0("group", uni_group)) +
        ggplot2::theme_bw()
    if (log.l) p3 <- p3 + ggplot2::xlab("Log Lambda") else p3 <- p3 + ggplot2::xlab("Lambda")
    
    
    p4 <- outputs4 %>% 
        ggplot2::ggplot(ggplot2::aes(x = .data$lambda,
                                     y = .data$value, 
                                     group = .data$variable, 
                                     col = .data$group)) +
        ggplot2::geom_line() +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::scale_color_discrete(name = "Group", labels = paste0("group", uni_group)) +
        ggplot2::xlab("beta norm / (beta norm") +
        ggplot2::ylab("Coefficients") +
        ggplot2::theme_bw()
        
    
    return(list(p1, p3, p2, p4))
}
    

