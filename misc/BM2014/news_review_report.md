## BM2014 News Comparison Report

This report compares the current R `news()` implementation and documentation to:
- the BM2014 theory (Section 2.3 + Appendix D), and
- the MATLAB reference (`News_DFM_ML.m`).

### Summary

The R implementation is theoretically aligned with BM2014 and the MATLAB logic. Differences are limited to numerical stabilization and output semantics. These differences are justified and do not alter the decomposition formula or interpretation.

### Theory Alignment (Paper)

BM2014 defines the news decomposition:
- Innovation: \u03bd_i = x_i^{new} - E[x_i | \u03a9_old]
- Revision: y_t^{new} - y_t^{old} = \u2211_i g_i \u03bd_i
- Gain: g = \u03c3_y C_y P1 P2^{-1}

R implements these components:
- Innovations computed from old-vintage smoothing.
- P1/P2 built using lag-augmented covariance blocks and measurement error R.
- Correct handling of lag direction (transpose when t_miss > t_fcst).
- Special cases handled: target observed, no releases.

### MATLAB vs R: Algorithm Comparison

Common elements:
- Same innovation definition and standardization.
- Same k selection: k = max(|lag|, max(lag)-min(lag)).
- Same P1/P2 covariance block indexing and WW adjustment.
- Same use of measurement error in P2.

Differences:
- R pre-computes P1/P2 for multiple targets (efficiency).
- R uses explicit re-scaling of old data to new-vintage stats.
- R supports MQ/AR1 via full state matrices (`ss_full`).

These differences are either efficiency improvements or required due to R’s compact dfm output.

### Output Semantics

R returns `news_df` with columns:
- `actual`: actual release (if any)
- `forecast`: old-vintage forecast of the release
- `news`: total innovation per series on output scale (equals actual-forecast if single release)
- `gain`: effective weight on `news` (output scale)
- `gain_std`: effective weight on standardized innovations
- `impact`: contribution to target revision (impact = news * gain)

This is a clarity improvement over MATLAB’s vector outputs and is consistent with BM2014’s formula.

### Numerical and Stability Differences

1) **Q stabilization in KFS**  
R nudges zero diagonal entries in Q to 1e-8 before KFS to avoid singularity under ragged edge patterns.  
This matches BM2014 practice in mixed-frequency initialization and is a justified numerical safeguard.

2) **Augmented state initialization**  
R uses stationary covariance for augmented states (plus small diagonal), while MATLAB uses zeros.  
This yields more stable smoothing without changing the theoretical decomposition.

### MQ / AR(1) Extensions

MATLAB’s `News_DFM_ML.m` assumes baseline DFM (no MQ, no AR1).  
R extends the decomposition to MQ and AR(1) by requiring full state matrices in `dfm$ss_full`.  
This is necessary and consistent with the BM2014 augmented state formulation.

### Practical Assessment

The R implementation is theoretically sound and matches BM2014’s decomposition logic.  
Differences are justified and improve stability and clarity.  
No substantive deviations were found that change the meaning of the news decomposition.

