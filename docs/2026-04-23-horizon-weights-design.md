# CPI Inflation Forecasting — Horizon-Specific Validation Weights

## Status
Draft — pending user approval

---

## Problem Statement

Current `main_horizon_weights_mbb_v4_ridge.R` uses a single-phase validation approach:
- **Phase 1**: NPREV=12 fixed rolling windows, each predicting 1 step → one unified weight vector
- **Phase 2**: Apply that single weight vector to all horizons h ∈ {1, 3, 6}

**Issue**: 1-step validation MSPE may not reflect 3-step or 6-step forecast ability. A model that is best at 1-step may be inferior at longer horizons. Applying 1-step weights to h=3,6 likely produces suboptimal combination forecasts.

**Example**: AR dominates 1-step but XGBoost dominates 3-step → using AR's high weight at h=3 hurts accuracy.

---

## Solution: Per-Horizon Validation Weights

For each h ∈ {1, 3, 6}:
1. Run rolling window validation where each window predicts **h steps directly**
2. Compute MSPE across all validation points for each model
3. Derive inverse-MSP weights specific to that horizon
4. Apply those weights at test time when forecasting h steps

---

## Architecture

### Data Split

```
Train:   1 ~ train_end
Val:     val_start ~ val_end   (H_val = 24 rows)
Test:    test_start ~ test_end  (N_test = 12 rows)
```

- `split_data()` computes boundaries unchanged
- Rolling window validation only touches the val window

### Phase 1: Per-Horizon Validation

**For each h in horizons:**

```
n_val_points = H_val - h + 1

for v in 1:n_val_points:
    train_start = train_end - NPREV + v
    train_end_w = train_end + v - 1        ← note: window slides WITH val point
    train_window = Y[train_start : train_end_w]

    actual_t = Y[train_end_w + h]           ← actual value h steps ahead

    for b in 1:B:
        mbb_sample = get_mbb_sample(train_window, l)
        for model in models:
            pred_h[b] = direct_predict_h(model, mbb_sample, train_window, h)

    MSPE_v[model] = mean((pred_h - actual_t)^2)

MSPE_h[model] = mean(MSPE_v[model])
weights_h = inverse_mspe_weights(MSPE_h)
```

**Validation point count by horizon:**
| h | n_val_points |
|---|-------------|
| 1 | 24 - 1 + 1 = 24 |
| 3 | 24 - 3 + 1 = 22 |
| 6 | 24 - 6 + 1 = 19 |

### Phase 2: Test — Unchanged Except Weight Selection

- Loop over horizons unchanged
- At each h: select `weights = weights_by_horizon[[h]]` instead of global weights
- All else identical to current v4_ridge

---

## Changes

### `functions_mbb_v4.R`

1. **`calc_block_length(n, h, embed_lag=3)`** — current signature already takes `h` but ignores it. Currently hardcoded to `2L * ceiling(n^(1/3))`. No change needed since block length formula is horizon-agnostic.

2. **No new functions needed** — existing `fit_*_direct` / `predict_*_direct_from_real` already handle arbitrary h via the direct prediction wrapper. Each already takes an `h` parameter.

### `main_horizon_weights_mbb_v4_ridge.R`

1. **Phase 1 outer loop**: Replace single 1-step MSPE computation with `for (h in horizons)` computing per-horizon weights.

2. **Phase 1 inner loop**: For each validation point `v`:
   - Compute `train_end_w = train_end + v - 1` (sliding)
   - Actual = `Y[train_end_w + h]`
   - Each model's h-step direct prediction uses current h (not hardcoded 1)

3. **Phase 2 weight selection**: `weights <- weights_by_horizon[[as.character(h)]]` — already correct in current v2, confirm v4_ridge uses it.

---

## Key Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Validation window slides? | Yes (`train_end + v - 1`) | Matches current v4_ridge behavior |
| Block length formula | Unchanged (horizon-agnostic) | MBB block length doesn't depend on h |
| Weight combination | Inverse MSPE per horizon | Keeps existing method, just applied at correct granularity |
| Embed lag | Unchanged (embed_lag=3) | No reason to change for this fix |
| B (bootstrap count) | B=5 (v4_ridge smoke test value) | No change from current |

---

## Edge Cases

1. **n_val_points too small**: For h=6, n_val_points=19. Acceptable (NPREV=12 is larger, gives 12 windows).

2. **Window size minimum**: `train_end_w - train_start + 1 = NPREV = 12` always. Each model requires at least `embed_lag + h` rows for direct features. For h=6: 12 - 6 - 3 = 3 rows minimum; if insufficient, skip that model for that validation point.

3. **Missing predictions**: If any model fails in MBB bootstrap, that model's MSPE for that validation point is set to NA and excluded from the mean (handled by `mean(..., na.rm=TRUE)`).

4. **Horizon 1 consistency**: When h=1, this reduces to the current behavior (n_val_points=24, identical structure). Results should match current v4_ridge.

---

## Testing Plan

1. **Smoke test**: Run with B=5, all 3 horizons. Confirm no crashes.
2. **Horizon 1 regression**: Compare h=1 weights with current v4_ridge — should be identical (same validation structure, same 1-step).
3. **Weight plausibility**: Inspect weights for h=3 and h=6 — confirm they differ from h=1 and are not trivially uniform.
4. **Evaluation comparison**: Compare test RMSE for h=3 and h=6 between current (single weights) and new (per-horizon weights) approaches.

---

## File Changes Summary

| File | Change |
|------|--------|
| `main_horizon_weights_mbb_v4_ridge.R` | Phase 1: outer loop over horizons; compute h-step MSPE per horizon; store `weights_by_horizon` |
| `functions_mbb_v4.R` | No changes required |
