# sparsegl (development version)

* Add the option to weight individual coefficients in the l1 penalty.
* Compute MSE internally in Fortran for `family = "Gaussian"`. Avoids the creation of a potentially large matrix of predicted values for the purposes
of risk estimation. 
* Revise `estimate_risk()` signature. Now `x` is optional and `y` is not required.

# sparsegl 0.3.0

* Initial version on CRAN.
