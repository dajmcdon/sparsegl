#' Simplified gene expression data from Alon et al. (1999)
#' 
#' Gene expression data (20 genes for 62 samples) from the microarray
#' experiments of colon tissue samples of Alon et al. (1999).
#' 
#' This data set contains 62 samples with 100 predictors (expanded from 20
#' genes using 5 basis B-splines, as described in Yang, Y. and Zou, H. (2015)):
#' 40 tumor tissues, coded 1 and 22 normal tissues, coded -1.
#' 
#' @return A list with the following elements: \item{x}{a [62 x 100] matrix
#' (expanded from a [62 x 20] matrix) giving the expression levels of 20 genes
#' for the 62 colon tissue samples. Each row corresponds to a patient, each 5
#' consecutive columns to a grouped gene.} \item{y}{a numeric vector of length
#' 62 giving the type of tissue sample (tumor or normal).}
#' @references Alon, U. and Barkai, N. and Notterman, D.A. and Gish, K. and
#' Ybarra, S. and Mack, D. and Levine, A.J. (1999). ``Broad patterns of gene
#' expression revealed by clustering analysis of tumor and normal colon tissues
#' probed by oligonucleotide arrays'', \emph{Proc. Natl. Acad. Sci. USA},
#' \bold{96}(12), 6745--6750.\cr
#' 
#' Yang, Y. and Zou, H. (2015), ``A Fast Unified Algorithm for Computing
#' Group-Lasso Penalized Learning Problems,'' \emph{Statistics and Computing}.
#' 25(6), 1129-1141.\cr BugReport: \url{https://github.com/emeryyi/gglasso}\cr
#' @source The data are described in Alon et al. (1999) and can be freely
#' downloaded from
#' \url{http://microarray.princeton.edu/oncology/affydata/index.html}.
#' @keywords datasets
#' @examples
#' 
#' # load gglasso library
#' library(gglasso)
#' 
#' # load data set
#' data(colon)
#' 
#' # how many samples and how many predictors ?
#' dim(colon$x)
#' 
#' # how many samples of class -1 and 1 respectively ?
#' sum(colon$y==-1)
#' sum(colon$y==1)
#' 
"colon"