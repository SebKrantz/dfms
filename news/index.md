# Changelog

## dfms 0.3.2

CRAN release: 2025-09-24

- Minor internal C++ changes to ensure compatibility with
  *RcppArmadillo* 15.0.2.

## dfms 0.3.1

CRAN release: 2025-08-20

- Fixed bug which occurred with only one quarterly variable
  ([\#73](https://github.com/SebKrantz/dfms/issues/73)). Thanks
  [@SantiagoD999](https://github.com/SantiagoD999) for reporting.

## dfms 0.3.0

CRAN release: 2025-05-18

- Added argument `quarterly.vars`, enabling mixed-frequency estimation
  with monthly and quarterly data following Banbura and Modugno (2014).
  The data matrix should contain the quarterly variables at the end
  (after the monthly ones).

## dfms 0.2.2

CRAN release: 2024-06-09

- Replace Armadillo `inv_sympd()` by Armadillo `inv()` in C++ Kalman
  Filter to improve numerical robustness at a minor performance cost.

## dfms 0.2.1

CRAN release: 2023-04-03

- Fixed print bug in `summary.dfm`: print method showed that model had
  AR(1) errors even though `idio.ar1 = FALSE` by default.

## dfms 0.2.0

CRAN release: 2023-03-31

- Added argument `idio.ar1 = TRUE` allowing estimation of approximate
  DFM’s with AR(1) observation errors.

- Added a small theoretical vignette entitled ‘Dynamic Factor Models: A
  Very Short Introduction’. This vignette lays a foundation for the
  present and future functionality of *dfms*. I plan to implement all
  features described in this vignette until summer 2023.

## dfms 0.1.5

- Added argument `na.keep = TRUE` to `fitted.dfm`. Setting
  `na.keep = FALSE` allows interpolation of data based on the DFM.
  Thanks [@apoorvalal](https://github.com/apoorvalal)
  ([\#45](https://github.com/SebKrantz/dfms/issues/45)).

## dfms 0.1.4

CRAN release: 2023-01-12

- Fixed minor bug in `summary.dfm` occurring if only one factor was
  estimated (basically an issue with dropping matrix dimensions which
  lead the factor summary statistics to be displayed without names).

## dfms 0.1.3

CRAN release: 2022-10-12

- Implemented some minor CRAN comments, no changes to functionality.

## dfms 0.1.2

- New default `em.method = "auto"`, which uses `"BM"` if the data has
  any missing values and `"DGR"` otherwise.

- Added vignette providing a walkthrough of the main features.

## dfms 0.1.1

- Renamed package from *DFM* to *dfms*. Lowercase names are preferred by
  [rOpenSci](https://devguide.ropensci.org/building.html?q=package%20name#package-name-and-metadata),
  and this also helps distinguish the package name from the main
  function [`DFM()`](https://sebkrantz.github.io/dfms/reference/DFM.md).
  The new name was inspired by the *vars* package.
