# sparsegl 1.1.1

* Add CITATION and links to JSS article

# sparsegl 1.1.0

* Address some CRAN issues for FORTRAN builds on GCC 14.1

# sparsegl 1.0.2

* Corrects a bug / adds support for `weights` and `offset` in cv
* Some minor tweaks to error and warning messages

# sparsegl 1.0.1

* Remove unnecessary attributes from the included data

# sparsegl 1.0.0

* Improved documentation
* Added functionality to implement arbitrary `stats::family()` fitting.
* Enhanced (and corrected) S3 interface for `predict` and `coef` methods.

# sparsegl 0.5.0

* Minor internal updates to pass additional CRAN checks
* Refactor some documentation
* Intercept calculation is internal to Fortran source in all cases.


# sparsegl 0.4.0

* Add the option to weight individual coefficients in the l1 penalty.
* Remove coercions of the type `as(<matrix>, "dgCMatrix")`. These are deprecated in `{Matrix}`>=1.4-2 and will Warn on CRAN checks.
* Compute MSE internally in Fortran for `family = "Gaussian"`. Avoids the creation of a potentially large matrix of predicted values for the purposes
of risk estimation. 
* Revise `estimate_risk()` signature. Now `x` is optional and `y` is not required.

# sparsegl 0.3.0

* Initial version on CRAN.
