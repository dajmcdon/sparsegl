#' Simplified gene expression data from Scheetz et al. (2006)
#' 
#' Gene expression data (20 genes for 120 samples) from the microarray
#' experiments of mammalian eye tissue samples of Scheetz et al. (2006).
#' 
#' This data set contains 120 samples with 100 predictors (expanded from 20
#' genes using 5 basis B-splines, as described in Yang, Y. and Zou, H. (2015)).
#' 
#' @return A list with the following elements: \item{x}{a [120 x 100] matrix
#' (expanded from a [120 x 20] matrix) giving the expression levels of 20
#' filtered genes for the 120 samples. Each row corresponds to a subject, each
#' 5 consecutive columns to a grouped gene.} \item{y}{a numeric vector of
#' length 120 giving expression level of gene TRIM32, which causes Bardet-Biedl
#' syndrome.}
#' @references Scheetz, T., Kim, K., Swiderski, R., Philp, A., Braun, T.,
#' Knudtson, K., Dorrance, A., DiBona, G., Huang, J., Casavant, T. et al.
#' (2006), ``Regulation of gene expression in the mammalian eye and its
#' relevance to eye disease'', \emph{Proceedings of the National Academy of
#' Sciences} \bold{103}(39), 14429-14434. \cr
#' 
#' Huang, J., S. Ma, and C.-H. Zhang (2008). ``Adaptive Lasso for sparse
#' high-dimensional regression models''. \emph{Statistica Sinica} 18,
#' 1603-1618.\cr
#' 
#' Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for Computing
#' Group-Lasso Penalized Learning Problems,'' \emph{Statistics and Computing}.
#' 25(6), 1129-1141.\cr BugReport: \url{https://github.com/emeryyi/gglasso}\cr
#' @keywords datasets
#' @examples
#' 
#' # load gglasso library
#' library(gglasso)
#' 
#' # load data set
#' data(bardet)
#' 
#' # how many samples and how many predictors ?
#' dim(bardet$x)
#' 
#' # repsonse y
#' bardet$y
#' 
"bardet"