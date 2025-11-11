# DFM Residuals and Fitted Values

The residuals \\\textbf{e}\_t = \textbf{x}\_t - \textbf{C}
\textbf{F}\_t\\ or fitted values \\\textbf{C} \textbf{F}\_t\\ of the DFM
observation equation.

## Usage

``` r
# S3 method for class 'dfm'
residuals(
  object,
  method = switch(object$em.method, none = "2s", "qml"),
  orig.format = FALSE,
  standardized = FALSE,
  na.keep = TRUE,
  ...
)

# S3 method for class 'dfm'
fitted(
  object,
  method = switch(object$em.method, none = "2s", "qml"),
  orig.format = FALSE,
  standardized = FALSE,
  na.keep = TRUE,
  ...
)
```

## Arguments

- object:

  an object of class 'dfm'.

- method:

  character. The factor estimates to use: one of `"qml"`, `"2s"` or
  `"pca"`.

- orig.format:

  logical. `TRUE` returns residuals/fitted values in a data format
  similar to `X`.

- standardized:

  logical. `FALSE` will put residuals/fitted values on the original data
  scale.

- na.keep:

  logical. `TRUE` inserts missing values where `X` is missing (default
  `TRUE` as residuals/fitted values are only defined for observed data).
  `FALSE` returns the raw prediction, which can be used to interpolate
  data based on the DFM. For residuals, `FALSE` returns the difference
  between the prediction and the initial imputed version of `X` use for
  PCA to initialize the Kalman Filter.

- ...:

  not used.

## Value

A matrix of DFM residuals or fitted values. If `orig.format = TRUE` the
format may be different, e.g. a data frame.

## See also

[dfms-package](https://sebkrantz.github.io/dfms/reference/dfms-package.md)

## Examples

``` r
# \donttest{
library(xts)
# Fit DFM with 3 factors and 3 lags in the transition equation
mod <- DFM(diff(BM14_M), r = 3, p = 3)
#> Converged after 26 iterations.

# Residuals
head(resid(mod))
#>      ip_total ip_tot_cstr ip_tot_cstr_en ip_constr ip_im_goods ip_capital
#> [1,]       NA          NA             NA        NA          NA         NA
#> [2,]       NA          NA             NA        NA          NA         NA
#> [3,]       NA          NA             NA        NA          NA         NA
#> [4,]       NA          NA             NA        NA          NA         NA
#> [5,]       NA          NA             NA        NA          NA         NA
#> [6,]       NA          NA             NA        NA          NA         NA
#>      ip_d_cstr ip_nd_cons ip_en ip_en_2 ip_manuf ip_metals ip_chemicals
#> [1,]        NA         NA    NA      NA       NA        NA           NA
#> [2,]        NA         NA    NA      NA       NA        NA           NA
#> [3,]        NA         NA    NA      NA       NA        NA           NA
#> [4,]        NA         NA    NA      NA       NA        NA           NA
#> [5,]        NA         NA    NA      NA       NA        NA           NA
#> [6,]        NA         NA    NA      NA       NA        NA           NA
#>      ip_electric ip_machinery ip_paper ip_plastic new_cars orders
#> [1,]          NA           NA       NA         NA       NA     NA
#> [2,]          NA           NA       NA         NA       NA     NA
#> [3,]          NA           NA       NA         NA       NA     NA
#> [4,]          NA           NA       NA         NA       NA     NA
#> [5,]          NA           NA       NA         NA       NA     NA
#> [6,]          NA           NA       NA         NA       NA     NA
#>      ret_turnover_defl ecs_ec_sent_ind ecs_ind_conf ecs_ind_order_book
#> [1,]        -3.6576879              NA           NA                 NA
#> [2,]         1.7727190              NA           NA                 NA
#> [3,]        -1.1002519              NA           NA                 NA
#> [4,]        -0.1454957              NA           NA                 NA
#> [5,]        -0.5699934              NA           NA                 NA
#> [6,]         1.6003246              NA           NA                 NA
#>      ecs_ind_stocks ecs_ind_prod_exp ecs_ind_prod_rec_m ecs_ind_x_orders
#> [1,]             NA               NA                 NA               NA
#> [2,]             NA               NA                 NA               NA
#> [3,]             NA               NA                 NA               NA
#> [4,]             NA               NA                 NA               NA
#> [5,]             NA               NA                 NA               NA
#> [6,]             NA               NA                 NA               NA
#>      ecs_ind_empl_exp ecs_cons_conf ecs_cons_sit_over_next_12
#> [1,]               NA            NA                        NA
#> [2,]               NA            NA                        NA
#> [3,]               NA            NA                        NA
#> [4,]               NA            NA                        NA
#> [5,]               NA            NA                        NA
#> [6,]               NA            NA                        NA
#>      ecs_cons_exp_unempl ecs_cons_gen_last_12m ecs_cstr_conf
#> [1,]                  NA                    NA            NA
#> [2,]                  NA                    NA            NA
#> [3,]                  NA                    NA            NA
#> [4,]                  NA                    NA            NA
#> [5,]                  NA                    NA            NA
#> [6,]                  NA                    NA            NA
#>      ecs_cstr_order_books ecs_cstr_empl_exp ecs_cstr_prod_recent
#> [1,]                   NA                NA                   NA
#> [2,]                   NA                NA                   NA
#> [3,]                   NA                NA                   NA
#> [4,]                   NA                NA                   NA
#> [5,]                   NA                NA                   NA
#> [6,]                   NA                NA                   NA
#>      ecs_ret_tr_conf ecs_ret_tr_bus_sit ecs_ret_tr_stocks ecs_ret_tr_exp_bus
#> [1,]              NA                 NA                NA                 NA
#> [2,]              NA                 NA                NA                 NA
#> [3,]              NA                 NA                NA                 NA
#> [4,]              NA                 NA                NA                 NA
#> [5,]              NA                 NA                NA                 NA
#> [6,]              NA                 NA                NA                 NA
#>      ecs_ret_tr_empl ecs_serv_conf ecs_serv_empl_exp pms_comp_output
#> [1,]              NA            NA                NA              NA
#> [2,]              NA            NA                NA              NA
#> [3,]              NA            NA                NA              NA
#> [4,]              NA            NA                NA              NA
#> [5,]              NA            NA                NA              NA
#> [6,]              NA            NA                NA              NA
#>      pms_comp_empl pms_pmi pms_manuf_empl pms_manuf_output pms_manuf_product
#> [1,]            NA      NA             NA               NA                NA
#> [2,]            NA      NA             NA               NA                NA
#> [3,]            NA      NA             NA               NA                NA
#> [4,]            NA      NA             NA               NA                NA
#> [5,]            NA      NA             NA               NA                NA
#> [6,]            NA      NA             NA               NA                NA
#>      pms_serv_out pms_serv_empl pms_serv_new_bus pms_serv_product urx
#> [1,]           NA            NA               NA               NA  NA
#> [2,]           NA            NA               NA               NA  NA
#> [3,]           NA            NA               NA               NA  NA
#> [4,]           NA            NA               NA               NA  NA
#> [5,]           NA            NA               NA               NA  NA
#> [6,]           NA            NA               NA               NA  NA
#>      empl_total empl_tot_xc empl_cstr empl_manuf extra_ea_trade_exp_val
#> [1,]         NA          NA        NA         NA             190582.858
#> [2,]         NA          NA        NA         NA            -370236.171
#> [3,]         NA          NA        NA         NA             318427.157
#> [4,]         NA          NA        NA         NA               5589.842
#> [5,]         NA          NA        NA         NA             250820.296
#> [6,]         NA          NA        NA         NA             806903.878
#>      intra_ea_trade_exp_val extra_ea_trade_imp_val intra_ea_trade_imp_val
#> [1,]               98385.12               594946.9               339686.3
#> [2,]              395717.71              1991220.5             -1795697.5
#> [3,]              972838.03             -1359308.8              1292021.9
#> [4,]              411894.61              1573422.8               527296.4
#> [5,]             1286616.09               674687.8               874156.6
#> [6,]              -19528.59               336060.0               813911.4
#>           us_ip      us_urx     us_empl us_retail_sales us_ip_manuf_exp
#> [1,] -0.2772139  0.04836678 -112.456851              NA        5.324118
#> [2,] -0.2644710 -0.02617550 -365.912104              NA       -3.224801
#> [3,] -0.6596791  0.41765710 -299.274986              NA      -15.436624
#> [4,] -0.8140125  0.39418780  -56.356502              NA      -14.602936
#> [5,] -0.4127089 -0.01475864 -209.665861              NA        7.100297
#> [6,] -0.4506433  0.17780271   -0.211103              NA        1.243730
#>      us_cons_exp    us_r3_m us_r10_year          m3 loans ir_long ir_short
#> [1,]   -1.122283  0.4861193   1.3943184 -0.08865737    NA      NA       NA
#> [2,]  -11.783175  2.1103821   0.1684485 -0.08831266    NA      NA       NA
#> [3,]    2.069787 -1.3568310  -0.9003209 -0.31414432    NA      NA       NA
#> [4,]    4.281491 -3.6692419  -0.7184414 -0.32873688    NA      NA       NA
#> [5,]    8.629724 -1.0978836  -0.1540432 -0.21035360    NA      NA       NA
#> [6,]   -0.395273  0.9553125   0.4662679 -0.10306526    NA      NA       NA
#>      ir_1_year ir_2_year ir_5_year eer eer_cpi eer_ppi       exr_usd
#> [1,]        NA        NA        NA  NA      NA      NA -0.0007111725
#> [2,]        NA        NA        NA  NA      NA      NA  0.0028584538
#> [3,]        NA        NA        NA  NA      NA      NA  0.0016526492
#> [4,]        NA        NA        NA  NA      NA      NA  0.0311893194
#> [5,]        NA        NA        NA  NA      NA      NA  0.0079394990
#> [6,]        NA        NA        NA  NA      NA      NA  0.0058606558
#>           exr_gbp    rxr_yen euro50 euro325         sp500       dow_j
#> [1,] -0.005337100   6.959517     NA      NA -11.574282006 -106.400843
#> [2,]  0.010539456  -3.996799     NA      NA  -9.330390459  -67.400779
#> [3,] -0.015277038   2.059095     NA      NA  25.722292255  195.046508
#> [4,] -0.018721359 -15.982046     NA      NA  37.754219317  292.312532
#> [5,] -0.009910865  -9.425426     NA      NA  19.593190651  143.727912
#> [6,] -0.005135135   7.589525     NA      NA  -0.007684972    3.769579
#>      raw_mat_en raw_mat_oil raw_mat_gold raw_mat_oil_fwd   raw_mat
#> [1,]  1.3983662  -12.206736    -5.636493              NA -4.886608
#> [2,] -2.1803877    8.016667   -91.962489              NA  1.512282
#> [3,]  0.3466068   12.951518   -36.349745              NA  6.457894
#> [4,]  2.7357011   11.407453   -11.553640              NA  5.039393
#> [5,] -0.4720561    2.670849    81.338684              NA  1.007499
#> [6,]  0.9367513   -5.982752    40.826881              NA -2.203898
plot(resid(mod, orig.format = TRUE)) # this is an xts object


# Fitted values
head(fitted(mod))
#>      ip_total ip_tot_cstr ip_tot_cstr_en ip_constr ip_im_goods ip_capital
#> [1,]       NA          NA             NA        NA          NA         NA
#> [2,]       NA          NA             NA        NA          NA         NA
#> [3,]       NA          NA             NA        NA          NA         NA
#> [4,]       NA          NA             NA        NA          NA         NA
#> [5,]       NA          NA             NA        NA          NA         NA
#> [6,]       NA          NA             NA        NA          NA         NA
#>      ip_d_cstr ip_nd_cons ip_en ip_en_2 ip_manuf ip_metals ip_chemicals
#> [1,]        NA         NA    NA      NA       NA        NA           NA
#> [2,]        NA         NA    NA      NA       NA        NA           NA
#> [3,]        NA         NA    NA      NA       NA        NA           NA
#> [4,]        NA         NA    NA      NA       NA        NA           NA
#> [5,]        NA         NA    NA      NA       NA        NA           NA
#> [6,]        NA         NA    NA      NA       NA        NA           NA
#>      ip_electric ip_machinery ip_paper ip_plastic new_cars orders
#> [1,]          NA           NA       NA         NA       NA     NA
#> [2,]          NA           NA       NA         NA       NA     NA
#> [3,]          NA           NA       NA         NA       NA     NA
#> [4,]          NA           NA       NA         NA       NA     NA
#> [5,]          NA           NA       NA         NA       NA     NA
#> [6,]          NA           NA       NA         NA       NA     NA
#>      ret_turnover_defl ecs_ec_sent_ind ecs_ind_conf ecs_ind_order_book
#> [1,]      0.0836697879              NA           NA                 NA
#> [2,]      0.1650839195              NA           NA                 NA
#> [3,]     -0.0168157840              NA           NA                 NA
#> [4,]     -0.0718749623              NA           NA                 NA
#> [5,]     -0.0367624807              NA           NA                 NA
#> [6,]      0.0007295562              NA           NA                 NA
#>      ecs_ind_stocks ecs_ind_prod_exp ecs_ind_prod_rec_m ecs_ind_x_orders
#> [1,]             NA               NA                 NA               NA
#> [2,]             NA               NA                 NA               NA
#> [3,]             NA               NA                 NA               NA
#> [4,]             NA               NA                 NA               NA
#> [5,]             NA               NA                 NA               NA
#> [6,]             NA               NA                 NA               NA
#>      ecs_ind_empl_exp ecs_cons_conf ecs_cons_sit_over_next_12
#> [1,]               NA            NA                        NA
#> [2,]               NA            NA                        NA
#> [3,]               NA            NA                        NA
#> [4,]               NA            NA                        NA
#> [5,]               NA            NA                        NA
#> [6,]               NA            NA                        NA
#>      ecs_cons_exp_unempl ecs_cons_gen_last_12m ecs_cstr_conf
#> [1,]                  NA                    NA            NA
#> [2,]                  NA                    NA            NA
#> [3,]                  NA                    NA            NA
#> [4,]                  NA                    NA            NA
#> [5,]                  NA                    NA            NA
#> [6,]                  NA                    NA            NA
#>      ecs_cstr_order_books ecs_cstr_empl_exp ecs_cstr_prod_recent
#> [1,]                   NA                NA                   NA
#> [2,]                   NA                NA                   NA
#> [3,]                   NA                NA                   NA
#> [4,]                   NA                NA                   NA
#> [5,]                   NA                NA                   NA
#> [6,]                   NA                NA                   NA
#>      ecs_ret_tr_conf ecs_ret_tr_bus_sit ecs_ret_tr_stocks ecs_ret_tr_exp_bus
#> [1,]              NA                 NA                NA                 NA
#> [2,]              NA                 NA                NA                 NA
#> [3,]              NA                 NA                NA                 NA
#> [4,]              NA                 NA                NA                 NA
#> [5,]              NA                 NA                NA                 NA
#> [6,]              NA                 NA                NA                 NA
#>      ecs_ret_tr_empl ecs_serv_conf ecs_serv_empl_exp pms_comp_output
#> [1,]              NA            NA                NA              NA
#> [2,]              NA            NA                NA              NA
#> [3,]              NA            NA                NA              NA
#> [4,]              NA            NA                NA              NA
#> [5,]              NA            NA                NA              NA
#> [6,]              NA            NA                NA              NA
#>      pms_comp_empl pms_pmi pms_manuf_empl pms_manuf_output pms_manuf_product
#> [1,]            NA      NA             NA               NA                NA
#> [2,]            NA      NA             NA               NA                NA
#> [3,]            NA      NA             NA               NA                NA
#> [4,]            NA      NA             NA               NA                NA
#> [5,]            NA      NA             NA               NA                NA
#> [6,]            NA      NA             NA               NA                NA
#>      pms_serv_out pms_serv_empl pms_serv_new_bus pms_serv_product urx
#> [1,]           NA            NA               NA               NA  NA
#> [2,]           NA            NA               NA               NA  NA
#> [3,]           NA            NA               NA               NA  NA
#> [4,]           NA            NA               NA               NA  NA
#> [5,]           NA            NA               NA               NA  NA
#> [6,]           NA            NA               NA               NA  NA
#>      empl_total empl_tot_xc empl_cstr empl_manuf extra_ea_trade_exp_val
#> [1,]         NA          NA        NA         NA              -305115.0
#> [2,]         NA          NA        NA         NA               198378.8
#> [3,]         NA          NA        NA         NA              -475081.2
#> [4,]         NA          NA        NA         NA              -171789.3
#> [5,]         NA          NA        NA         NA              -750484.2
#> [6,]         NA          NA        NA         NA              -813549.6
#>      intra_ea_trade_exp_val extra_ea_trade_imp_val intra_ea_trade_imp_val
#> [1,]               44046.37               58442.46             -152565.57
#> [2,]             -130482.95              390811.63               91613.75
#> [3,]             -974592.57             -671158.98             -732883.48
#> [4,]             -713775.15             -652894.14             -476253.32
#> [5,]             -808583.16             -752361.67             -824700.12
#> [6,]             -496135.16             -542413.81             -708081.06
#>            us_ip      us_urx    us_empl us_retail_sales us_ip_manuf_exp
#> [1,]  0.32921387 -0.04836678  228.45685              NA        4.675882
#> [2,]  0.09347103  0.02617550   83.91210              NA        1.224801
#> [3,] -0.50632094  0.18234290 -180.72501              NA       -6.563376
#> [4,] -0.61198747  0.20581220 -231.64350              NA       -9.397064
#> [5,] -0.27529107  0.11475864  -53.33414              NA       -3.100297
#> [6,]  0.06064331  0.02219729  114.21110              NA        1.756270
#>      us_cons_exp     us_r3_m us_r10_year        m3 loans ir_long ir_short
#> [1,]   1.9222826  0.37201468  0.21814568 0.2019563    NA      NA       NA
#> [2,]   1.1831749  0.23425448  0.16586228 0.2536198    NA      NA       NA
#> [3,]  -1.9697874 -0.64078803 -0.37539336 0.3837109    NA      NA       NA
#> [4,]  -3.3814910 -0.95694856 -0.57346338 0.4475182    NA      NA       NA
#> [5,]  -0.9297239 -0.40497354 -0.24929015 0.3269930    NA      NA       NA
#> [6,]   0.7952730  0.02676539  0.00336416 0.2450060    NA      NA       NA
#>      ir_1_year ir_2_year ir_5_year eer eer_cpi eer_ppi     exr_usd      exr_gbp
#> [1,]        NA        NA        NA  NA      NA      NA -0.01243519 -0.007376398
#> [2,]        NA        NA        NA  NA      NA      NA -0.07538607 -0.019812599
#> [3,]        NA        NA        NA  NA      NA      NA -0.01265303  0.007510884
#> [4,]        NA        NA        NA  NA      NA      NA  0.02422289  0.019962715
#> [5,]        NA        NA        NA  NA      NA      NA  0.01281248  0.009991140
#> [6,]        NA        NA        NA  NA      NA      NA  0.01182098  0.003544841
#>         rxr_yen euro50 euro325      sp500      dow_j raw_mat_en raw_mat_oil
#> [1,] -0.9472834     NA      NA  16.176270  123.87371  2.2016338    4.306736
#> [2,] -7.7484387     NA      NA  -1.323895   -7.25715  3.5803877   -7.416667
#> [3,] -3.3943497     NA      NA -27.482682 -212.27318 -3.0466068  -12.251518
#> [4,] -0.3183934     NA      NA -32.856037 -250.45063 -5.5357011  -11.207453
#> [5,] -0.2184031     NA      NA -12.872888 -102.06125 -2.7279439   -4.670849
#> [6,]  0.7041374     NA      NA   5.220149   36.15176 -0.5367513    2.082752
#>      raw_mat_gold raw_mat_oil_fwd    raw_mat
#> [1,]    -3.050195              NA  3.3856082
#> [2,]   -19.654177              NA  2.1901851
#> [3,]    -1.221445              NA -5.7574272
#> [4,]     9.116140              NA -8.2415268
#> [5,]     5.179530              NA -3.7092994
#> [6,]     4.137405              NA  0.1024977
head(fitted(mod, orig.format = TRUE)) # this is an xts object
#>            ip_total ip_tot_cstr ip_tot_cstr_en ip_constr ip_im_goods ip_capital
#> 1980-01-31       NA          NA             NA        NA          NA         NA
#> 1980-02-29       NA          NA             NA        NA          NA         NA
#> 1980-03-31       NA          NA             NA        NA          NA         NA
#> 1980-04-30       NA          NA             NA        NA          NA         NA
#> 1980-05-31       NA          NA             NA        NA          NA         NA
#> 1980-06-30       NA          NA             NA        NA          NA         NA
#>            ip_d_cstr ip_nd_cons ip_en ip_en_2 ip_manuf ip_metals ip_chemicals
#> 1980-01-31        NA         NA    NA      NA       NA        NA           NA
#> 1980-02-29        NA         NA    NA      NA       NA        NA           NA
#> 1980-03-31        NA         NA    NA      NA       NA        NA           NA
#> 1980-04-30        NA         NA    NA      NA       NA        NA           NA
#> 1980-05-31        NA         NA    NA      NA       NA        NA           NA
#> 1980-06-30        NA         NA    NA      NA       NA        NA           NA
#>            ip_electric ip_machinery ip_paper ip_plastic new_cars orders
#> 1980-01-31          NA           NA       NA         NA       NA     NA
#> 1980-02-29          NA           NA       NA         NA       NA     NA
#> 1980-03-31          NA           NA       NA         NA       NA     NA
#> 1980-04-30          NA           NA       NA         NA       NA     NA
#> 1980-05-31          NA           NA       NA         NA       NA     NA
#> 1980-06-30          NA           NA       NA         NA       NA     NA
#>            ret_turnover_defl ecs_ec_sent_ind ecs_ind_conf ecs_ind_order_book
#> 1980-01-31                NA              NA           NA                 NA
#> 1980-02-29        0.08366979              NA           NA                 NA
#> 1980-03-31        0.16508392              NA           NA                 NA
#> 1980-04-30       -0.01681578              NA           NA                 NA
#> 1980-05-31       -0.07187496              NA           NA                 NA
#> 1980-06-30       -0.03676248              NA           NA                 NA
#>            ecs_ind_stocks ecs_ind_prod_exp ecs_ind_prod_rec_m ecs_ind_x_orders
#> 1980-01-31             NA               NA                 NA               NA
#> 1980-02-29             NA               NA                 NA               NA
#> 1980-03-31             NA               NA                 NA               NA
#> 1980-04-30             NA               NA                 NA               NA
#> 1980-05-31             NA               NA                 NA               NA
#> 1980-06-30             NA               NA                 NA               NA
#>            ecs_ind_empl_exp ecs_cons_conf ecs_cons_sit_over_next_12
#> 1980-01-31               NA            NA                        NA
#> 1980-02-29               NA            NA                        NA
#> 1980-03-31               NA            NA                        NA
#> 1980-04-30               NA            NA                        NA
#> 1980-05-31               NA            NA                        NA
#> 1980-06-30               NA            NA                        NA
#>            ecs_cons_exp_unempl ecs_cons_gen_last_12m ecs_cstr_conf
#> 1980-01-31                  NA                    NA            NA
#> 1980-02-29                  NA                    NA            NA
#> 1980-03-31                  NA                    NA            NA
#> 1980-04-30                  NA                    NA            NA
#> 1980-05-31                  NA                    NA            NA
#> 1980-06-30                  NA                    NA            NA
#>            ecs_cstr_order_books ecs_cstr_empl_exp ecs_cstr_prod_recent
#> 1980-01-31                   NA                NA                   NA
#> 1980-02-29                   NA                NA                   NA
#> 1980-03-31                   NA                NA                   NA
#> 1980-04-30                   NA                NA                   NA
#> 1980-05-31                   NA                NA                   NA
#> 1980-06-30                   NA                NA                   NA
#>            ecs_ret_tr_conf ecs_ret_tr_bus_sit ecs_ret_tr_stocks
#> 1980-01-31              NA                 NA                NA
#> 1980-02-29              NA                 NA                NA
#> 1980-03-31              NA                 NA                NA
#> 1980-04-30              NA                 NA                NA
#> 1980-05-31              NA                 NA                NA
#> 1980-06-30              NA                 NA                NA
#>            ecs_ret_tr_exp_bus ecs_ret_tr_empl ecs_serv_conf ecs_serv_empl_exp
#> 1980-01-31                 NA              NA            NA                NA
#> 1980-02-29                 NA              NA            NA                NA
#> 1980-03-31                 NA              NA            NA                NA
#> 1980-04-30                 NA              NA            NA                NA
#> 1980-05-31                 NA              NA            NA                NA
#> 1980-06-30                 NA              NA            NA                NA
#>            pms_comp_output pms_comp_empl pms_pmi pms_manuf_empl
#> 1980-01-31              NA            NA      NA             NA
#> 1980-02-29              NA            NA      NA             NA
#> 1980-03-31              NA            NA      NA             NA
#> 1980-04-30              NA            NA      NA             NA
#> 1980-05-31              NA            NA      NA             NA
#> 1980-06-30              NA            NA      NA             NA
#>            pms_manuf_output pms_manuf_product pms_serv_out pms_serv_empl
#> 1980-01-31               NA                NA           NA            NA
#> 1980-02-29               NA                NA           NA            NA
#> 1980-03-31               NA                NA           NA            NA
#> 1980-04-30               NA                NA           NA            NA
#> 1980-05-31               NA                NA           NA            NA
#> 1980-06-30               NA                NA           NA            NA
#>            pms_serv_new_bus pms_serv_product urx empl_total empl_tot_xc
#> 1980-01-31               NA               NA  NA         NA          NA
#> 1980-02-29               NA               NA  NA         NA          NA
#> 1980-03-31               NA               NA  NA         NA          NA
#> 1980-04-30               NA               NA  NA         NA          NA
#> 1980-05-31               NA               NA  NA         NA          NA
#> 1980-06-30               NA               NA  NA         NA          NA
#>            empl_cstr empl_manuf extra_ea_trade_exp_val intra_ea_trade_exp_val
#> 1980-01-31        NA         NA                     NA                     NA
#> 1980-02-29        NA         NA              -305115.0               44046.37
#> 1980-03-31        NA         NA               198378.8             -130482.95
#> 1980-04-30        NA         NA              -475081.2             -974592.57
#> 1980-05-31        NA         NA              -171789.3             -713775.15
#> 1980-06-30        NA         NA              -750484.2             -808583.16
#>            extra_ea_trade_imp_val intra_ea_trade_imp_val       us_ip
#> 1980-01-31                     NA                     NA          NA
#> 1980-02-29               58442.46             -152565.57  0.32921387
#> 1980-03-31              390811.63               91613.75  0.09347103
#> 1980-04-30             -671158.98             -732883.48 -0.50632094
#> 1980-05-31             -652894.14             -476253.32 -0.61198747
#> 1980-06-30             -752361.67             -824700.12 -0.27529107
#>                 us_urx    us_empl us_retail_sales us_ip_manuf_exp us_cons_exp
#> 1980-01-31          NA         NA              NA              NA          NA
#> 1980-02-29 -0.04836678  228.45685              NA        4.675882   1.9222826
#> 1980-03-31  0.02617550   83.91210              NA        1.224801   1.1831749
#> 1980-04-30  0.18234290 -180.72501              NA       -6.563376  -1.9697874
#> 1980-05-31  0.20581220 -231.64350              NA       -9.397064  -3.3814910
#> 1980-06-30  0.11475864  -53.33414              NA       -3.100297  -0.9297239
#>               us_r3_m us_r10_year        m3 loans ir_long ir_short ir_1_year
#> 1980-01-31         NA          NA        NA    NA      NA       NA        NA
#> 1980-02-29  0.3720147   0.2181457 0.2019563    NA      NA       NA        NA
#> 1980-03-31  0.2342545   0.1658623 0.2536198    NA      NA       NA        NA
#> 1980-04-30 -0.6407880  -0.3753934 0.3837109    NA      NA       NA        NA
#> 1980-05-31 -0.9569486  -0.5734634 0.4475182    NA      NA       NA        NA
#> 1980-06-30 -0.4049735  -0.2492902 0.3269930    NA      NA       NA        NA
#>            ir_2_year ir_5_year eer eer_cpi eer_ppi     exr_usd      exr_gbp
#> 1980-01-31        NA        NA  NA      NA      NA          NA           NA
#> 1980-02-29        NA        NA  NA      NA      NA -0.01243519 -0.007376398
#> 1980-03-31        NA        NA  NA      NA      NA -0.07538607 -0.019812599
#> 1980-04-30        NA        NA  NA      NA      NA -0.01265303  0.007510884
#> 1980-05-31        NA        NA  NA      NA      NA  0.02422289  0.019962715
#> 1980-06-30        NA        NA  NA      NA      NA  0.01281248  0.009991140
#>               rxr_yen euro50 euro325      sp500      dow_j raw_mat_en
#> 1980-01-31         NA     NA      NA         NA         NA         NA
#> 1980-02-29 -0.9472834     NA      NA  16.176270  123.87371   2.201634
#> 1980-03-31 -7.7484387     NA      NA  -1.323895   -7.25715   3.580388
#> 1980-04-30 -3.3943497     NA      NA -27.482682 -212.27318  -3.046607
#> 1980-05-31 -0.3183934     NA      NA -32.856037 -250.45063  -5.535701
#> 1980-06-30 -0.2184031     NA      NA -12.872888 -102.06125  -2.727944
#>            raw_mat_oil raw_mat_gold raw_mat_oil_fwd   raw_mat
#> 1980-01-31          NA           NA              NA        NA
#> 1980-02-29    4.306736    -3.050195              NA  3.385608
#> 1980-03-31   -7.416667   -19.654177              NA  2.190185
#> 1980-04-30  -12.251518    -1.221445              NA -5.757427
#> 1980-05-31  -11.207453     9.116140              NA -8.241527
#> 1980-06-30   -4.670849     5.179530              NA -3.709299
# }
```
