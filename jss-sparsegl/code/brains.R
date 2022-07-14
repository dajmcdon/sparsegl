
# Load brain data ---------------------------------------------------------

library(sparsegl)
library(Matrix)

A <- Matrix::readMM('jss-sparsegl/large-data/A.mtx.gz')
Y <- readRDS('jss-sparsegl/large-data/Y.rds')
G <- readRDS('jss-sparsegl/large-data/G.rds')
gpf <- readRDS('jss-sparsegl/large-data/Gpf.rds')

# Fit and estimate risk ---------------------------------------------------

system.time(fit <- sparsegl(A, Y, group=G, pf_group=gpf, asparse=0.0)) # ~ 1.5 minutes

df <- estimate_risk(fit, A, approx_df=TRUE)
write.csv(df, "jss-sparsegl/large-data/fit.csv", row.names=FALSE, quote=FALSE)
saveRDS(fit, file = "jss-sparsegl/large-data/brain-fit.rds")
