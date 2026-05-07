# Final pipeline results — 2026-05-07

Source: msliu1996 local pipeline run, delivered via DM 03:46:55 UTC.

## Variant & weighting

- **Variant A**: inverse-MSPE^k weighting with **k=2** (w_i ∝ 1/MSPE_i²)
- **Variant B (geo)**: geometric combination variant; only h=11 shown in this delivery

## Files

- `weights_summary.xlsx` — variant A combination weights, horizons 1-4 and 11-12 (h=5..10 missing from this delivery; to be confirmed with msliu1996)
- `mspe_data_summary.csv` — validation-set MSPE by horizon, 5 base models. Covers h=11,12 only.
- `predictions_A_h11.csv` — test-set predictions at h=11 (5 base + Combined), RMSE / MAE / SMAPE / Theil's U1 / MDAPE. 2 evaluation points.
- `predictions_A_h12.csv` — test-set predictions at h=12. 1 evaluation point.
- `predictions_geo_A_h11.csv` — variant B (geometric) per-forecast predictions at h=11.

## Headline numbers

Combined RMSE vs best single model:

| h | Combined | Best single | Winner |
|---|---|---|---|
| 1 | 0.234 | XGBoost 0.217 | XGBoost |
| 2 | 0.211 | Combined 0.211 | Combined |
| 3 | 0.184 | FAVAR 0.174 | FAVAR |
| 4 | 0.169 | DFM 0.154 | DFM |
| 11 | 0.059 | FAVAR 0.108 | **Combined (-46%)** |
| 12 | 0.076 | DFM 0.142 | **Combined (-47%)** |

Long-horizon combination gains are substantial.

## Caveats

- Test set has 12 monthly points total; at h=11 only 2 evaluation points exist, at h=12 only 1 point. Small-sample noise dominates these estimates.
- MDAPE values in predictions files are inflated by small denominators (see paper §4 caption footnote).
- No RW baseline / no Theil's U2 by msliu1996 decision (variant A evaluation uses U1).

## Open questions

- Full h=5..10 coverage: confirm with msliu1996 whether pipeline ran or results are separately delivered.
- Variant B (geo): only h=11 predictions shown; full h range + MSPE for geo to be confirmed.
