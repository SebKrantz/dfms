# Armadillo's Inverse Functions

Matrix inverse and pseudo-inverse by the Armadillo C++ library.

## Usage

``` r
ainv(x)

apinv(x)
```

## Arguments

- x:

  a numeric matrix, must be square for `ainv`.

## Value

The matrix-inverse or pseudo-inverse.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
ainv(crossprod(diff(EuStockMarkets)))
#>                DAX           SMI           CAC          FTSE
#> DAX   1.671464e-06 -5.411415e-07 -7.533943e-07 -3.323765e-07
#> SMI  -5.411415e-07  8.147114e-07 -1.596419e-07 -1.770041e-07
#> CAC  -7.533943e-07 -1.596419e-07  2.013654e-06 -4.872422e-07
#> FTSE -3.323765e-07 -1.770041e-07 -4.872422e-07  1.234721e-06
```
