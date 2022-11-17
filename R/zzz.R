## Funs loaded before the package

.onLoad <- function(libname, pkgname) {

  ## Propagate all warnings from the Matrix package to the checks.
  ## Coercions like `as(<matrix>, "dgCmatrix")` no longer supported.
  options(Matrix.warnDeprecatedCoerce = 2)
}
