"""Print compact header info for each of the seven products."""
import sys
import json
import netCDF4
import numpy as np

for path in sys.argv[1:]:
    with netCDF4.Dataset(path) as ds:
        print("=" * 20, path.split("/")[-1])
        print(" dims:", {k: len(v) for k, v in ds.dimensions.items()})
        for name, v in ds.variables.items():
            extra = ""
            if "time" in name.lower() and v.ndim == 1:
                vals = v[:]
                try:
                    units = v.units
                    t0 = netCDF4.num2date(vals[0], units)
                    t1 = netCDF4.num2date(vals[-1], units)
                    if len(vals) > 1:
                        step = float(np.diff(np.asarray(vals[:1000], dtype="f8")).min())
                    else:
                        step = None
                    extra = f" range {t0}..{t1} step~{step} units={units}"
                except Exception as e:
                    extra = f" (time decode failed: {e})"
            print(f"  {name} {v.dtype} {v.dimensions} fill={getattr(v,'_FillValue',None)}{extra}")
        print(" global attrs:", json.dumps({k: str(ds.getncattr(k))[:80] for k in ds.ncattrs()}, indent=1)[:800])
