"""Simple diagnostic plot; no manuscript-style changes to original figures."""
from __future__ import annotations
from pathlib import Path
import numpy as np


def comparison_plot(computed, published, outfile: Path, label: str, rmse):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    computed=np.asarray(computed,dtype=float)
    published=np.asarray(published,dtype=float)
    use=np.isfinite(computed)&np.isfinite(published)
    if not np.any(use):return None
    fig,ax=plt.subplots(figsize=(5.4,5))
    hi=float(max(np.nanmax(computed[use]),np.nanmax(published[use]),1.0))
    ax.hexbin(published[use],computed[use],gridsize=70,mincnt=1,bins='log',cmap='viridis')
    ax.plot([0,hi],[0,hi],'k--',lw=1)
    ax.set(xlabel='Published CML rainfall (mm h$^{-1}$)',
           ylabel='Recomputed CML rainfall (mm h$^{-1}$)',
           xlim=(0,hi),ylim=(0,hi),title=f'{label}\nRMSE = {rmse:.4g} mm h$^{{-1}}$' if rmse is not None else label)
    fig.tight_layout();outfile.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(outfile,dpi=180);plt.close(fig)
    return outfile
