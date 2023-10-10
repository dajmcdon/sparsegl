library("sparsegl")
library("Matrix")

# Download and load brain data ---------------------------------------------------------
download <- function(url, path) {
    file.exists(path) || download.file(url, path, mode = 'wb')
}

options(timeout = max(15 * 60)) # 15 minutes, assuming 1MB/s download speed

# https://doi.org/10.6084/m9.figshare.20314917
download("https://figshare.com/ndownloader/files/36288819", "jss-sparsegl/large-data/A.mtx.gz")
download("https://figshare.com/ndownloader/files/36288825", "jss-sparsegl/large-data/Y.rds")
download("https://figshare.com/ndownloader/files/36288822", "jss-sparsegl/large-data/G.rds")
download("https://figshare.com/ndownloader/files/36294165", "jss-sparsegl/large-data/Gpf.rds")
download("https://figshare.com/ndownloader/files/36288828", "jss-sparsegl/large-data/brain-fit.rds")

A <- Matrix::readMM("jss-sparsegl/large-data/A.mtx.gz")
Y <- readRDS("jss-sparsegl/large-data/Y.rds")
G <- readRDS("jss-sparsegl/large-data/G.rds")
gpf <- readRDS("jss-sparsegl/large-data/Gpf.rds")

# Fit and estimate risk ---------------------------------------------------

system.time(fit <- sparsegl(A, Y, group = G, pf_group = gpf, asparse = 0.0)) # ~ 1.5 minutes

df <- estimate_risk(fit, A, approx_df = TRUE)
# write.csv(df, "jss-sparsegl/large-data/fit.csv", row.names = FALSE, quote = FALSE)
saveRDS(fit, file = "jss-sparsegl/large-data/brain-fit.rds")
