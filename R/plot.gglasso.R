#' Plot solution paths from a "gglasso" object
#' 
#' Produces a coefficient profile plot of the coefficient paths for a fitted
#' \code{\link{gglasso}} object.
#' 
#' A coefficient profile plot is produced.
#' 
#' @param x fitted \code{\link{gglasso}} model
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
#' @examples
#' 
#' # load gglasso library
#' library(gglasso)
#' 
#' # load data set
#' data(bardet)
#' 
#' # define group index
#' group <- rep(1:20,each=5)
#' 
#' # fit group lasso
#' m1 <- gglasso(x=bardet$x,y=bardet$y,group=group,loss="ls")
#' 
#' # make plots
#' par(mfrow=c(1,3))
#' plot(m1) # plots the coefficients against the log-lambda sequence 
#' plot(m1,group=TRUE) # plots group norm against the log-lambda sequence 
#' plot(m1,log.l=FALSE) # plots against the lambda sequence
#' @method plot gglasso
#' @export
plot.gglasso <- function(x, group = FALSE, log.l = TRUE, ...) {
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
    p <- nrow(tmp)
    l <- x$lambda
    n.g <- max(g)
    
    if(group){
        bs <- as.integer(as.numeric(table(g)))
        ix <- rep(NA, n.g)
        iy <- rep(NA, n.g)
        j <- 1
        for (g in 1:n.g) {
            ix[g] <- j
            iy[g] <- j + bs[g] - 1
            j <- j + bs[g]
        }
        beta <- matrix(NA, n.g, length(l))
        for (g in 1:n.g) {
            crossp <- apply(tmp[ix[g]:iy[g], ], 2, crossprod)
            beta[g, ] <- sqrt(crossp)
        }
    } else beta <- tmp

    if (log.l) {
        l <- log(l)
        xlab <- "Log Lambda"
    } else xlab <- "Lambda"
    
    plot.args <- list(x = l, y = 1:length(l), ylim = range(beta), xlab = xlab, 
        ylab = "Coefficients", type = "n", xlim = range(l))
    new.args <- list(...)
    if (length(new.args)) {
        new.plot.args <- new.args[names(new.args) %in% c(names(par()), names(formals(plot.default)))]
        plot.args[names(new.plot.args)] <- new.plot.args
    }
    do.call("plot", plot.args)
    line.args <- list(col = rainbow(n.g + 1, start = 0.7, end = 0.95)[1:n.g], 
        lwd = 1 + 1.2^(-p/20), lty = 1)
    
    if (length(new.args)) 
        line.args[names(new.args)] <- new.args
    line.args$x <- l
    line.args$y <- t(beta)
    line.args$col <- rep(line.args$col, table(g))
    do.call("matlines", line.args)
    
    abline(h = 0, lwd = line.args$lwd)
} 
