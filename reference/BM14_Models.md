# Euro Area Macroeconomic Data from Banbura and Modugno 2014

A data extract from BM 2014 replication files. Some proprietary series
(mostly PMI's) are excluded. The dataset `BM14_Models` provides
information about all series and their inclusion in the 'small',
'medium' and 'large' sized dynamic factor models estimated by BM 2014.
The actual data is contained in *xts* format in `BM14_M` for monthly
data and `BM14_Q` for quarterly data.

## Usage

``` r
BM14_Models
BM14_M
BM14_Q
```

## Format

`BM14_Models` is a data frame with 101 obs. (series) and 8 columns:

- series:

  BM14 series code (converted to snake case for R)

- label:

  BM14 series label

- code:

  original series code from data source

- freq:

  series frequency

- log_trans:

  logical indicating whether the series was transformed by the natural
  log before differencing. Note that all data are provided in
  untransformed levels, and all data was (log-)differenced by BM14
  before estimation.

- small:

  logical indicating series included in the 'small' model of BM14.
  Proprietary series are excluded.

- medium:

  logical indicating series included in the 'medium' model of BM14.
  Proprietary series are excluded.

- large:

  logical indicating series included in the 'large' model of BM14. This
  comprises all series, thus the variable is redundant but included for
  completeness. Proprietary series are excluded.

## Source

Banbura, M., & Modugno, M. (2014). Maximum likelihood estimation of
factor models on datasets with arbitrary pattern of missing data.
*Journal of Applied Econometrics, 29*(1), 133-160.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
library(magrittr)
library(xts)
#> Loading required package: zoo
#> 
#> Attaching package: ‘zoo’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.Date, as.Date.numeric

# Constructing the database for the large model
BM14 = merge(BM14_M, BM14_Q)
BM14[, BM14_Models$log_trans] %<>% log()
BM14[, BM14_Models$freq == "M"] %<>% diff()
BM14[, BM14_Models$freq == "Q"] %<>% diff(3)

# Small Model Database
head(BM14[, BM14_Models$small])
#>            ip_tot_cstr new_cars orders ret_turnover_defl ecs_ec_sent_ind
#> 1980-01-31          NA       NA     NA                NA              NA
#> 1980-02-29          NA       NA     NA      -0.045867752              NA
#> 1980-03-31          NA       NA     NA       0.025129863              NA
#> 1980-04-30          NA       NA     NA      -0.014409195              NA
#> 1980-05-31          NA       NA     NA      -0.002828186              NA
#> 1980-06-30          NA       NA     NA      -0.007937021              NA
#>            pms_pmi urx extra_ea_trade_exp_val euro325      raw_mat          gdp
#> 1980-01-31      NA  NA                     NA      NA           NA           NA
#> 1980-02-29      NA  NA           -0.006628411      NA -0.017094433           NA
#> 1980-03-31      NA  NA           -0.010029199      NA  0.041649238           NA
#> 1980-04-30      NA  NA           -0.009230433      NA  0.007688121           NA
#> 1980-05-31      NA  NA           -0.009886927      NA -0.035638515           NA
#> 1980-06-30      NA  NA           -0.030326833      NA -0.031090587 -0.004706623
#>                    empl capacity      gdp_us
#> 1980-01-31           NA       NA          NA
#> 1980-02-29           NA       NA          NA
#> 1980-03-31           NA       NA          NA
#> 1980-04-30           NA       NA          NA
#> 1980-05-31           NA       NA          NA
#> 1980-06-30 0.0005918654       NA -0.02070885

# Medium-Sized Model Database
head(BM14[, BM14_Models$medium])
#>            ip_tot_cstr ip_constr ip_im_goods ip_capital ip_d_cstr ip_nd_cons
#> 1980-01-31          NA        NA          NA         NA        NA         NA
#> 1980-02-29          NA        NA          NA         NA        NA         NA
#> 1980-03-31          NA        NA          NA         NA        NA         NA
#> 1980-04-30          NA        NA          NA         NA        NA         NA
#> 1980-05-31          NA        NA          NA         NA        NA         NA
#> 1980-06-30          NA        NA          NA         NA        NA         NA
#>            ip_en new_cars orders ret_turnover_defl ecs_ec_sent_ind ecs_ind_conf
#> 1980-01-31    NA       NA     NA                NA              NA           NA
#> 1980-02-29    NA       NA     NA      -0.045867752              NA           NA
#> 1980-03-31    NA       NA     NA       0.025129863              NA           NA
#> 1980-04-30    NA       NA     NA      -0.014409195              NA           NA
#> 1980-05-31    NA       NA     NA      -0.002828186              NA           NA
#> 1980-06-30    NA       NA     NA      -0.007937021              NA           NA
#>            ecs_ind_prod_exp ecs_ind_x_orders ecs_ind_empl_exp ecs_cons_conf
#> 1980-01-31               NA               NA               NA            NA
#> 1980-02-29               NA               NA               NA            NA
#> 1980-03-31               NA               NA               NA            NA
#> 1980-04-30               NA               NA               NA            NA
#> 1980-05-31               NA               NA               NA            NA
#> 1980-06-30               NA               NA               NA            NA
#>            ecs_cstr_conf ecs_ret_tr_conf ecs_serv_conf ecs_serv_empl_exp
#> 1980-01-31            NA              NA            NA                NA
#> 1980-02-29            NA              NA            NA                NA
#> 1980-03-31            NA              NA            NA                NA
#> 1980-04-30            NA              NA            NA                NA
#> 1980-05-31            NA              NA            NA                NA
#> 1980-06-30            NA              NA            NA                NA
#>            pms_pmi pms_serv_out urx empl_total extra_ea_trade_exp_val
#> 1980-01-31      NA           NA  NA         NA                     NA
#> 1980-02-29      NA           NA  NA         NA           -0.006628411
#> 1980-03-31      NA           NA  NA         NA           -0.010029199
#> 1980-04-30      NA           NA  NA         NA           -0.009230433
#> 1980-05-31      NA           NA  NA         NA           -0.009886927
#> 1980-06-30      NA           NA  NA         NA           -0.030326833
#>            extra_ea_trade_imp_val         us_ip us_ip_manuf_exp us_cons_exp
#> 1980-01-31                     NA            NA              NA          NA
#> 1980-02-29            0.033595314  0.0008954556              10         0.8
#> 1980-03-31            0.113723637 -0.0029476935              -2       -10.6
#> 1980-04-30           -0.096103676 -0.0203344816             -22         0.1
#> 1980-05-31            0.044716794 -0.0254443045             -24         0.9
#> 1980-06-30           -0.003696892 -0.0125115289               4         7.7
#>                     m3 loans ir_long ir_short eer     exr_usd euro325
#> 1980-01-31          NA    NA      NA       NA  NA          NA      NA
#> 1980-02-29 0.007535875    NA      NA       NA  NA -0.01314636      NA
#> 1980-03-31 0.010894232    NA      NA       NA  NA -0.07252762      NA
#> 1980-04-30 0.004549395    NA      NA       NA  NA -0.01100038      NA
#> 1980-05-31 0.007720339    NA      NA       NA  NA  0.05541221      NA
#> 1980-06-30 0.007523561    NA      NA       NA  NA  0.02075198      NA
#>                  dow_j  raw_mat_en raw_mat_oil_fwd          gdp  priv_cons
#> 1980-01-31          NA          NA              NA           NA         NA
#> 1980-02-29  0.02009644  0.03875620              NA           NA         NA
#> 1980-03-31 -0.08884307  0.01467532              NA           NA         NA
#> 1980-04-30 -0.02167110 -0.02849797              NA           NA         NA
#> 1980-05-31  0.05186822 -0.03043713              NA           NA         NA
#> 1980-06-30  0.04908563 -0.03595893              NA -0.004706623 -0.0061348
#>                 invest      export      import         empl prductivity
#> 1980-01-31          NA          NA          NA           NA          NA
#> 1980-02-29          NA          NA          NA           NA          NA
#> 1980-03-31          NA          NA          NA           NA          NA
#> 1980-04-30          NA          NA          NA           NA          NA
#> 1980-05-31          NA          NA          NA           NA          NA
#> 1980-06-30 -0.01693449 -0.04461216 -0.02727455 0.0005918654          NA
#>            capacity      gdp_us
#> 1980-01-31       NA          NA
#> 1980-02-29       NA          NA
#> 1980-03-31       NA          NA
#> 1980-04-30       NA          NA
#> 1980-05-31       NA          NA
#> 1980-06-30       NA -0.02070885
```
