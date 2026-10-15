"""Pure-array validation metrics; shared by tests and rainfall comparison."""
import numpy as np


def rain_metrics(computed, reference, atol=1e-4, rtol=1e-4):
    """Report both numeric errors AND differing missing-value masks.

    Validity mismatch is not hidden by dropping NaNs before computing RMSE.
    """
    computed=np.asarray(computed,dtype=float)
    reference=np.asarray(reference,dtype=float)
    if computed.shape!=reference.shape:raise ValueError('Arrays must have equal shapes')
    a=np.isfinite(computed);b=np.isfinite(reference);shared=a&b
    mismatch=int(np.count_nonzero(a^b))
    d=computed[shared]-reference[shared]
    ok=np.isclose(computed[shared],reference[shared],atol=atol,rtol=rtol)
    return {
        'n_total':int(computed.size),'finite_computed':int(a.sum()),
        'finite_reference':int(b.sum()),'n_shared_finite':int(shared.sum()),
        'validity_mismatches':mismatch,
        'out_of_tolerance':int((~ok).sum()),
        'rmse_mm_h':float(np.sqrt(np.mean(d*d))) if len(d) else None,
        'mean_bias_mm_h':float(np.mean(d)) if len(d) else None,
        'max_abs_error_mm_h':float(np.max(np.abs(d))) if len(d) else None,
        'numeric_pass':bool(len(d)>0 and mismatch==0 and np.all(ok)),
        'atol':atol,'rtol':rtol,
    }


def flag_metrics(computed, reference):
    """Compare flags, including missing sentinels (not just 0/1 values)."""
    a=np.asarray(computed);b=np.asarray(reference)
    if a.shape!=b.shape:raise ValueError('Flag arrays must have equal shapes')
    eq=(a==b) | (np.isnan(a.astype(float)) & np.isnan(b.astype(float)))
    return {'n_total':int(a.size),'n_disagree':int((~eq).sum()),'pass':bool(np.all(eq))}
