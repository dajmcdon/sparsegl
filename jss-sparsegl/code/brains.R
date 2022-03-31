
# Load brain data ---------------------------------------------------------

library(here)
library(sparsegl)
library(Matrix)
source(here("..", "Pestilli","brain-tensor-data","tensor-processing.R"))

phi <- readRDS(here("..", "Pestilli", "data", "phi-tensor.rds"))
dict <- readRDS(here("..","Pestilli", "data", "dictionary.rds"))
y <- readRDS(here("..","Pestilli", "data", "y-vector.rds"))
gr <- readRDS(here("..","Pestilli", "data", "groups.rds"))
meta <- readRDS(here("..","Pestilli", "data", "meta-info.rds"))


# Build up the model components -------------------------------------------

ijk <- as.matrix(phi[,1:3])
val <- dplyr::pull(phi[,4])
rm(phi)
dims <- c(ncol(dict), meta$nvoxels, meta$nstreamlines)
Xunfold <- unfold(ijk, val, dims)
p <- meta$nstreamlines
n <- meta$nangles * meta$nvoxels
X <- dense_sp_mult(dict, Xunfold, n, p)
rm(Xunfold, dict)

# Fit and estimate risk ---------------------------------------------------

system.time(fit <- sparsegl(X, y, gr, asparse = 0)) # ~ two minutes

saveRDS(fit, file = "jss-sparsegl/large-data/brain-fit.rds")
saveRDS(meta, file = "jss-sparsegl/large-data/brain-meta.rds")
