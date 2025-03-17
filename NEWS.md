# dfms 0.3.0

* Added argument `quarterly.vars`, enabling mixed-frequency estimation with monthly and quarterly data following Banbura and Modugno (2014). The data matrix should contain the quarterly variables at the end (after the monthly ones).

# dfms 0.2.2

* Replace Armadillo `inv_sympd()` by Armadillo `inv()` in C++ Kalman Filter to improve numerical robustness at a minor performance cost.

# dfms 0.2.1

* Fixed print bug in `summary.dfm`: print method showed that model had AR(1) errors even though `idio.ar1 = FALSE` by default. 

# dfms 0.2.0

* Added argument `idio.ar1 = TRUE` allowing estimation of approximate DFM's with AR(1) observation errors. 

* Added a small theoretical vignette entitled 'Dynamic Factor Models: A Very Short Introduction'. This vignette lays a foundation for the present and future functionality of *dfms*. I plan to implement all features described in this vignette until summer 2023. 

# dfms 0.1.5

* Added argument `na.keep = TRUE` to `fitted.dfm`. Setting `na.keep = FALSE` allows interpolation of data based on the DFM. Thanks @apoorvalal (#45).

# dfms 0.1.4

* Fixed minor bug in `summary.dfm` occurring if only one factor was estimated (basically an issue with dropping matrix dimensions which lead the factor summary statistics to be displayed without names).

# dfms 0.1.3

* Implemented some minor CRAN comments, no changes to functionality. 

# dfms 0.1.2

* New default `em.method = "auto"`, which uses `"BM"` if the data has any missing values and `"DGR"` otherwise. 

* Added vignette providing a walkthrough of the main features. 

# dfms 0.1.1

* Renamed package from *DFM* to *dfms*. Lowercase names are preferred by [rOpenSci](<https://devguide.ropensci.org/building.html?q=package%20name#package-name-and-metadata>), and this also helps distinguish the package name from the main function `DFM()`. The new name was inspired by the *vars* package. 
