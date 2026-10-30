"""Radar rainfall along the CML paths from DWD RADKLIM YW 2017.002 (5 minutes).

Runs the steps of main() of the unchanged original/code/radar/dwd_radolan_processor.py with configurable paths:
  1. (with --download) download the RADKLIM YW 2017.002 netCDF archives for 2021 from DWD Open Data, including the
     2021 supplement files published by DWD, into extras.radklim_dir, and extract them;
  2. intersection weights of every CML path with the radar grid (poligrain, path length weights);
  3. path-averaged RR time series for every daily radar file found under extras.radklim_dir, with
     poligrain get_grid_time_series_at_intersections as in process_single_radar_file(); the weights stay in memory
     (sparse) instead of the processor's temporary weights file, which current sparse versions cannot write, and
     the series is named rr, the name merge_intersection_files() expects;
  4. merge of the daily files (merge_intersection_files).
The CML paths are those of the released radar file (1132 paths, endpoint coordinates).
Outputs in <output_dir>/radar: RADKLIM_5min/*_CML_RR.nc (daily) and cml_rr_RADKLIM_2021_CZ.nc (merged).
Usage: python scripts/radar_cml_paths.py [--config config/user_paths.json] [--download]
Needs poligrain and tqdm. The 2021 archives are about 15 GB.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.stdout.reconfigure(encoding='utf-8', errors='replace')  # the processor prints symbols
sys.path.insert(0, str(ROOT / 'src'))
from essd_cz.config import ORIGINAL, extra, load_config  # noqa: E402

# imported by name so that the worker processes can import it too
sys.path.insert(0, str(ORIGINAL / 'radar'))
import dwd_radolan_processor as m  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
    ap.add_argument('--download', action='store_true', help='download and extract the RADKLIM archives first')
    ap.add_argument('--radklim-dir', default=None, help='default: extras.radklim_dir of the configuration')
    a = ap.parse_args()
    cfg = load_config(a.config)
    radklim = str(Path(a.radklim_dir).resolve() if a.radklim_dir else extra(a.config, 'radklim_dir'))
    out = cfg['output_dir'] / 'radar'
    out.mkdir(parents=True, exist_ok=True)
    if a.download:
        m.download_radolan_data(years=m.YEARS, base_url=m.BASE_URL, root_dir=radklim)
        m.extract_all_archives(root_dir=radklim)
    sample = m.find_sample_radar_file(root_dir=radklim)
    if not sample:
        raise SystemExit(f'No radar .nc files under {radklim} (use --download)')
    weights, _ = m.calculate_intersection_weights(str(cfg['datasets']['radar']), sample)
    import glob, os
    import poligrain as plg
    import xarray as xr
    daily = out / 'RADKLIM_5min'
    daily.mkdir(exist_ok=True)
    files = sorted(glob.glob(f'{radklim}/**/*.nc', recursive=True))
    files = [f for f in files if 'RADKLIM' not in f and '_temp_' not in f]
    for k, f in enumerate(files, 1):
        target = daily / os.path.basename(f).replace('.nc', '_CML_RR.nc')
        if target.exists():
            continue
        with xr.open_dataset(f) as ds:
            ts = plg.spatial.get_grid_time_series_at_intersections(grid_data=ds['RR'], intersect_weights=weights)
            ts.rename('rr').to_netcdf(target)  # merge_intersection_files() expects the name rr
        print(f'{k}/{len(files)} {target.name}', flush=True)
    m.merge_intersection_files(input_dir=str(out / 'RADKLIM_5min'), output_file=str(out / 'cml_rr_RADKLIM_2021_CZ.nc'))


if __name__ == '__main__':
    main()
