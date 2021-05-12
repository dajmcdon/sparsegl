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
#' @param \dots other graphical parameters to plot
#' @author Yi Yang and Hui Zou\cr Maintainer: Yi Yang <yi.yang6@@mcgill.ca>
#' @references Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for
#' Computing Group-Lasso Penalized Learning Problems,'' \emph{Statistics and
#' Computing}. 25(6), 1129-1141.\cr BugReport:
#' \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords models regression
#' @method plot sparsegl
#' @export
plot.sparsegl <- function(x, group = FALSE, log.l = TRUE, alpha = 0.05, ...) {
    
    xb <- x$beta   # xb is the betas of output
    
    # only have one variable b1 or not
    # 'nonzeros': a vector recording the variables which are not zero w.r.t at least one lambda
    if (nrow(xb) == 1) {
        if (any(abs(xb) > 0)) {
            nonzeros <- 1
        } else nonzeros <- NULL
    } else {
        nonzeros <- which(apply(abs(xb), 1, sum) > 0)
    }
    
    tmp <- xb[nonzeros, , drop = FALSE]  # 'tmp' is the matrix filtered by 'nonzeros'
    g <- as.numeric(as.factor(x$group[nonzeros]))  # 'g': the group numbers of nonzero variables 
    p <- nrow(tmp)  # 'p': the number of nonzero variables
    l <- x$lambda  # 'l': the vectors of all lambdas
    n.g <- max(g)  # 'n.g': the largest group number
    
    if(group){
        bs <- as.integer(as.numeric(table(g)))  # 'bs': the number of nonzero variables within each group
        ix <- rep(NA, n.g)
        iy <- rep(NA, n.g)
        j <- 1
        for (g in 1:n.g) {
            ix[g] <- j
            iy[g] <- j + bs[g] - 1
            j <- j + bs[g]
        }
        # ix, iy: vectors indicating the position of each group
        
        beta <- matrix(NA, n.g, length(l))   # beta: #groups-by-#lambdas matrix
        for (g in 1:n.g) {
            # crossp <- apply(tmp[ix[g]:iy[g], ], 2, crossprod)
            # beta[g, ] <- sqrt(crossp)
            crossp <- apply(tmp[ix[g]:iy[g], ], 2, function(x) {alpha * sum(abs(x)) + (1 - alpha) * sqrt(crossprod(x))})
            beta[g,] <- crossp
        }
    } else beta <- tmp

    if (log.l) {
        l <- log(l)
        xlab <- "Log Lambda"
    } else xlab <- "Lambda"

    outputs <- as.data.frame(t(as.matrix(beta)))
    outputs$lambda <- l
    if(group){
        names <- c(1:(n.g+1))
        for (i in 1:n.g){
            names[i] <- paste("group", i)
        }
        names[n.g + 1] <- "lambda"
    }else{
        names <- c(1: (dim(outputs)[2]))
        for (i in 1:(dim(outputs)[2] - 1)){
            names[i] <- paste("variable", i)
        }
        names[dim(outputs)[2]] <- "lambda"
    }
    colnames(outputs) <- names
    outputs1 <- reshape2::melt(outputs, id.vars = 'lambda')

    # ggplot
    p1 <- outputs1 %>% ggplot2::ggplot(aes(x = lambda, y = value, col = variable)) +
        geom_line()+
        geom_hline(yintercept = 0)

    if (group){
        p1 <- p1 + xlab("Log Lambda")
    }else{
        p1 <- p1 + xlab("Lambda")
    }


    # standardized group norm on x-axis
    sgnorm <- apply(xb, 2, function(x) sp_group_norm(x, gr))
    outputs2 <- outputs
    outputs2$lambda <- sgnorm / max(sgnorm)
    outputs2 <- reshape2::melt(outputs2, id.vars = 'lambda')

    p2 <- outputs2 %>% ggplot2::ggplot(aes(x = lambda, y = value, col = variable)) +
        geom_line() +
        geom_hline(yintercept = 0) +
        xlab("Standardized Lambda")
    
    
    
    betas <- as.data.frame(t(as.matrix(xb)))
    betas$lambda <- l
    outputs3 <- reshape2::melt(betas, id.vars = "lambda")

    p3 <- outputs3 %>% ggplot2::ggplot(aes(x = lambda, y = value, col = variable)) +
        geom_line() +
        geom_hline(yintercept = 0)
    
    if (group){
        p3 <- p3 + xlab("Log Lambda")
    }else{
        p3 <- p3 + xlab("Lambda")
    }
    
    sgnorm <- apply(xb, 2, function(x) sp_group_norm(x, gr))
    outputs4 <- betas
    outputs4$lambda <- sgnorm / max(sgnorm)
    outputs4 <- reshape2::melt(outputs4, id.vars = 'lambda')

    p4 <- outputs4 %>% ggplot2::ggplot(aes(x = lambda, y = value, col = variable)) +
        geom_line() +
        geom_hline(yintercept = 0) +
        xlab("Standardized Lambda")
    
    names <- c(1:(dim(betas)[2]-1))
    if (group){
        bs <- as.integer(as.numeric(table(x$group)))  # 'bs': the number of nonzero variables within each group
        ix <- rep(NA, n.g)
        iy <- rep(NA, n.g)
        for (g in 1:length(bs)) {
            ix[g] <- j
            iy[g] <- j + bs[g] - 1
            j <- j + bs[g]
        }
        j <- 1
        for (g in 1:length(bs)){
            for (i in ix[g] : iy[g]){
                names[j] <- paste("group", g)
                j <- j + 1
            }
        }
    }else{
        for (j in 1:(dim(betas)[2]-1)){
            names[j] <- paste("variable", j)
        }
    }
    p3 <- p3 + scale_colour_hue(name = "Group", breaks=colnames(betas)[-length(colnames(betas))],labels=names)
    p4 <- p4 + scale_colour_hue(name = "Group", breaks=colnames(betas)[-length(colnames(betas))],labels=names)

    return(list(p1, p2, p3, p4))
    
} 

# plot.args <- list(x = l, y = 1:length(l), ylim = range(beta), xlab = xlab,
#     ylab = "Coefficients", type = "n", xlim = range(l))
# new.args <- list(...)
# if (length(new.args)) {
#     new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
#     plot.args[names(new.plot.args)] <- new.plot.args
# }
# do.call("plot", plot.args)
# line.args <- list(col = rainbow(n.g + 1, start = 0.7, end = 0.95)[1:n.g],
#     lwd = 1 + 1.2^(-p/20), lty = 1)
# 
# if (length(new.args))
#     line.args[names(new.args)] <- new.args
# line.args$x <- l
# line.args$y <- t(beta)
# line.args$col <- rep(line.args$col, table(g))
# do.call("matlines", line.args)
# 
# abline(h = 0, lwd = line.args$lwd)

