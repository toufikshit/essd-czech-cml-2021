# QC processing parameters: revision 1 (`QC_CORRECTED_CANDIDATE_REV1.nc`)

**Script:** `scripts\s3_04_qc_corrected_rev1.py`, sha256 3765739a7319c3cad7d73509ce97bc888bb827400c3e0f3c6e69c4d246fe1004
(stored in the NetCDF). It differs from `s3_01_qc_corrected.py` only in the blackout eligibility, two added counters
and the metadata text (`scripts\s3_04_diff_vs_s3_01.txt`).

**Decision:** the author's, of 2026-10-28 (`PROCESSING_LOG.md`).

**Everything else** is exactly as in `QC_PROCESSING_PARAMETERS.md`, which is preserved unchanged:
* period;
* population;
* P5 and P6;
* invalid values;
* plateau;
* interpolation;
* fluctuation;
* detection limit;
* sublink and CML rules.

## Blackout reconstruction (revised): adaptation of Blettner et al. (2023)

| Element | Definition |
|---|---|
| Blackout gap | A whole run of missing RX of ≤ 60 minutes, after invalid-value, P6 and plateau masking, with RX < −65 dBm on at least one edge |
| **Eligibility (revision 1)** | The gap is reconstructed only if it contains **no** minute masked by P6 (`qc_flag` bit 2) and **no** minute removed by the plateau filter (bit 4). Otherwise it is rejected: `trsl` stays missing and the minutes keep their flags |
| Reconstructed value | `trsl = max(tx) − min(rx)`, where max(tx) and min(rx) are taken over the same sublink's observations in the QC reference period 2021-04-01 00:00 – 2021-10-31 23:59 UTC, after masking. One value per sublink |
| What is written | Only `trsl` (bit 8). TX and RX are never filled. **RX is missing at every bit-8 minute; TX may be observed** (the gap is defined on RX) |
| Use | Excluded from the fluctuation statistics; never used as an interpolation end point |
| Difference from Blettner et al. (2023) and standalone v3 | They write max TX and min RX into the TX and RX series and fill any qualifying gap, including plateau-removed minutes. Here, no TX or RX value is created, and gaps with plateau-removed or P6-masked minutes are left missing |

## Flag bits (`qc_flag`)

| Bit | Meaning |
|---|---|
| 1 | Invalid value masked |
| 2 | Mixed-placeholder minute masked (P6) |
| 4 | Minute held a TX or RX value removed by the plateau filter |
| 8 | Blackout total loss reconstructed in `trsl` (RX missing, TX may be observed) |
| 16 | Short-gap interpolation of `trsl` (TX and/or RX missing) |

At bit-8 and bit-16 minutes `trsl` ≠ `tx − rx`.
