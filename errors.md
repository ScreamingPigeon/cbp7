# Potential Bugs / Issues

## TageAheadHC_IR.hpp

### 1. Overly conservative pred_ram update (line 1058)

```cpp
do_pred_update = t_m1[I] & t_phw & table_wrong;
```

Only updates provider's prediction when hyst is weak (`t_phw`) AND table was wrong.
Standard TAGE updates provider on any misprediction regardless of hyst strength.
This may miss training opportunities on strong-but-wrong entries.

### 2. Missing `.fo1()` on tag_hit pipeline shift (line 514)

```cpp
current_tag_hit[I] = prefetch_tag_hit[I];  // no .fo1()
```

Other pipeline shifts (current_idx, current_ctag, etc.) use `.fo1()`. This one
doesn't, risking unintended fanout propagation from the prefetch stage.

### 3. FB bimodal hyst one-way ratchet (line ~960)

`bh_gate = do_train & t_m1[NT] & mispredict & bh_changed`

Bimodal hysteresis is only updated on mispredict cycles. So hyst can go
strong(1)→weak(0) when wrong, but never weak(0)→strong(1) on correct predictions.
Once weak, the FB prediction bit becomes permanently flippable.

Same issue exists in template TageAhead.hpp (line 2249).

## TageAhead.hpp (Template)

### 4. FB bimodal hyst one-way ratchet (line 2249)

Same as #3 above. `fb_bh_read` is only read inside `execute_if(mispredict, ...)`
(line 1715), so the hyst RAM is never accessed on correct predictions, preventing
any strong-transition update.

### 5. fb_gate fanout over-declaration (line 2237)

`fb_gate.fanout(hard<5>{})` assumes FB_BANKS=2 write structure. For FB_BANKS=1
(single execute_if), only 1 fanout needed. Over-declaration is harmless but
misleading for audits.
