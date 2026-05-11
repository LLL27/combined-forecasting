# Final results — 2026-05-11

Pipeline run from msliu1996 local, delivered via DM in 4 waves 2026-05-11 ~12:10-12:20 UTC.

## Files (20 total)

### Summary tables (full h=1..12 coverage)
- `MSPE_summary.xlsx` — validation-set MSPE, 5 base models × 12 horizons + best model per horizon
- `evaluation.xlsx` — test-set metrics (RMSE/MAE/SMAPE/Theil's U1/MDAPE), 6 models (5 base + Combined) × 12 horizons
- `weights_summary.xlsx` — combination weights by horizon

### Per-forecast predictions (geometric inverse-MSPE, α=2)
- `predictions_geo_A_h{1..12}.xlsx` — 12 files, each contains Forecast / Actual / AR / FAVAR / DFM / XGBoost / NNAR / Combined

### Metrics summary (redundant with evaluation.xlsx; preserved for completeness)
- `predictions_A_h{4..8}.xlsx` — 5 files with per-horizon model metrics

## Configuration (per `main_horizon_weights_new_parallel_for_mac.R`)
- Val window: H_val = 24 months (rolling)
- Test: N_test = 12 months (h=1..12)
- Bootstrap: B = 100 replications (MBB)
- Weighting: `compute_geometric_mspe_weights` with α=2 (equivalent to inverse-MSPE²)
- Base models: AR / FAVAR (K=4, plag=2, nrep=3000) / DFM (r=5, p=1) / XGBoost (embed_lag=3, nrounds=1000, eta=0.01, max_depth=4) / NNAR (size=24, repeats=50)

## Known issue
- `mspe_val_h{h}.csv` not output locally due to newline bug in `main_horizon_weights_new_parallel_for_mac.R` file path. Fix pending.

## Next
- Step 0 (diagnostic): residual correlation + error covariance, reads `predictions_geo_A_h{1..12}.xlsx`
- Step 1 (weighting schemes): 8 alternative weighting methods, reads `MSPE_summary.xlsx` + `predictions_geo_A_h{1..12}.xlsx`
