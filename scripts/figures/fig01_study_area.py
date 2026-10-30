"""Figure 1: observation networks (CML paths, radar coverage, rain gauges, weather stations).

Data:
  * CML paths: the 2574 CMLs of the released RAW file (endpoint coordinates);
  * rain gauges: the 424 stations of the precipitation file;
  * weather stations: the union of the temperature (251) and humidity (247) station positions;
  * radar coverage: RADKLIM YW 2017.002 cells with finite data in at least one of the daily files
    1 April to 13 July 2021 (extras.radklim_dir of the configuration; DWD Climate Data Center).
Map: Natural Earth borders and rivers via cartopy, Europe inset, scale bar.
Output: <output_dir>/figures/f01_study_area.png/.pdf.
Usage: python scripts/figures/fig01_study_area.py [--config config/user_paths.json] [--coverage-cache FILE.npz]
"""
import argparse, sys
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np, xarray as xr
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
from essd_cz.config import extra, load_config  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--config', default=None, help='configuration file (default config/user_paths.json)')
ap.add_argument('--radar-nc', nargs='+', default=None,
                help='daily RADKLIM files (default: YW_2017.002_20210401..20210713 in extras.radklim_dir)')
ap.add_argument('--coverage-cache', default=None, help='npz with the coverage mask and lon/lat; written if missing, read if present')
ap.add_argument('--out', default=None)
a = ap.parse_args()
cfg = load_config(a.config)
a.cml_nc, a.gauge = str(cfg['datasets']['raw']), str(cfg['datasets']['gauge'])
a.temp, a.hum = str(cfg['datasets']['temperature']), str(cfg['datasets']['humidity'])
if not a.out:
    (cfg['output_dir'] / 'figures').mkdir(parents=True, exist_ok=True)
    a.out = str(cfg['output_dir'] / 'figures' / 'f01_study_area.png')
if a.radar_nc is None and not (a.coverage_cache and Path(a.coverage_cache).exists()):
    days = {f'YW_2017.002_{d:%Y%m%d}.nc' for d in __import__('pandas').date_range('2021-04-01', '2021-07-13')}
    a.radar_nc = sorted((str(f) for f in extra(a.config, 'radklim_dir').rglob('YW_2017.002_2021*.nc') if f.name in days),
                        key=lambda f: Path(f).name)
    if len(a.radar_nc) != len(days):
        sys.exit(f'{len(a.radar_nc)} of {len(days)} daily RADKLIM files found in extras.radklim_dir')
a.radar_nc = a.radar_nc or []

C_CML, C_RADAR, C_GAUGE, C_STATION = '#222222', '#d6e6f2', '#d95f02', '#1f6fb0'
EXT = [11.9, 19.05, 48.4, 51.2]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})

qc = xr.open_dataset(a.cml_nc)[['site_a_longitude', 'site_a_latitude', 'site_b_longitude', 'site_b_latitude']]
qc = qc.sel(cml_id=sorted(qc['cml_id'].values.astype(str)))
segs = np.stack([np.c_[qc.site_a_longitude.values, qc.site_a_latitude.values],
                 np.c_[qc.site_b_longitude.values, qc.site_b_latitude.values]], axis=1)
g = xr.open_dataset(a.gauge)
glon = np.asarray(g['lon'].values, float); glat = np.asarray(g['lat'].values, float)
st = []
for f in (a.temp, a.hum):
    d = xr.open_dataset(f)
    st.append(np.c_[np.asarray(d['lon'].values, float), np.asarray(d['lat'].values, float)])
st = np.unique(np.round(np.vstack(st), 6), axis=0)

# radar coverage: RADKLIM polar stereographic grid -> lon/lat of cells with finite data on the day
valid = None
cache = Path(a.coverage_cache) if a.coverage_cache else None
if cache is not None and cache.exists():
    z = np.load(cache); valid, rlon, rlat = z['valid'], z['lon'], z['lat']; a.radar_nc = []
    print('coverage read from', cache, 'days', int(z['ndays']))
for f in a.radar_nc:
    rd = xr.open_dataset(f)
    var = next(v for v in ['YW', 'RR', 'rainfall_rate', 'R', 'rainrate', 'precip'] if v in rd)
    if valid is None:
        valid = np.zeros(rd[var].shape[1:], bool)
    for t0 in range(0, rd.sizes['time'], 24):   # read in time chunks (memory)
        valid |= np.isfinite(np.asarray(rd[var].isel(time=slice(t0, t0 + 24)).values)).any(axis=0)
if a.radar_nc and {'lat', 'lon'} <= set(rd.variables):
    rlon, rlat = np.asarray(rd['lon'].values), np.asarray(rd['lat'].values)
elif a.radar_nc:
    X, Y = np.meshgrid(rd['x'].values, rd['y'].values)
    stp = ccrs.Stereographic(central_latitude=90, central_longitude=10, true_scale_latitude=60,
                             globe=ccrs.Globe(ellipse=None, semimajor_axis=6370040, semiminor_axis=6370040))
    p = ccrs.PlateCarree().transform_points(stp, X, Y); rlon, rlat = p[..., 0], p[..., 1]

if cache is not None and a.radar_nc:
    np.savez_compressed(cache, valid=valid, lon=rlon, lat=rlat, ndays=len(a.radar_nc))
proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(13.6, 6.75))
ax = fig.add_axes([0.05, 0.06, 0.92, 0.76], projection=proj)
ax.set_extent(EXT, crs=proj)
ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='#f3f3ef', edgecolor='none', zorder=0)
ax.add_feature(cfeature.OCEAN.with_scale('10m'), facecolor='#eaf2f8', edgecolor='none', zorder=0)
ax.contourf(rlon, rlat, valid.astype(float), levels=[0.5, 1.5], colors=[C_RADAR], alpha=0.55, transform=proj, zorder=1)
ax.add_feature(cfeature.RIVERS.with_scale('10m'), edgecolor='#a9cbe6', linewidth=0.6, zorder=1.5)
ax.add_feature(cfeature.BORDERS.with_scale('10m'), edgecolor='#555555', linewidth=0.7, zorder=2)
cz = [r for r in shpreader.Reader(shpreader.natural_earth('10m', 'cultural', 'admin_0_countries')).records()
      if r.attributes.get('ADM0_A3') == 'CZE']
ax.add_geometries([r.geometry for r in cz], crs=proj, facecolor='none', edgecolor='#222222', linewidth=1.4, zorder=2.5)
ax.add_collection(LineCollection(segs, colors=C_CML, linewidths=0.75, alpha=0.9, transform=proj, zorder=3))
ax.scatter(glon, glat, marker='^', s=46, c=C_GAUGE, edgecolors='white', linewidths=0.6, transform=proj, zorder=4)
ax.scatter(st[:, 0], st[:, 1], marker='o', s=58, c=C_STATION, edgecolors='white', linewidths=0.8, transform=proj, zorder=5)
ax.gridlines(crs=proj, draw_labels=False, linestyle='--', linewidth=0.5, color='0.6', alpha=0.8,
             xlocs=np.arange(12, 20, 1), ylocs=np.arange(48.5, 51.6, 0.5))
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
ax.set_xticks([12, 14, 16, 18], crs=proj); ax.set_yticks([49, 50, 51], crs=proj)
ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol='°')); ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol='°'))
ax.tick_params(labelsize=12, direction='out', length=5)
for sp in ax.spines.values():
    sp.set_linewidth(1.2); sp.set_edgecolor('#333333')

# scale bar (km at 49 N)
lat0, x0, y0 = 49.0, 17.05, 48.68
dlon = 100.0 / (111.32 * np.cos(np.deg2rad(lat0)))
ax.plot([x0, x0 + dlon], [y0, y0], color='k', lw=2.2, transform=proj, zorder=6)
for k, lab in ((0, '0'), (0.5, '50'), (1, '100 km')):
    ax.plot([x0 + k * dlon] * 2, [y0 - 0.025, y0 + 0.025], color='k', lw=1.6, transform=proj, zorder=6)
    ax.text(x0 + k * dlon, y0 + 0.07, lab, ha='center', va='bottom', fontsize=12, transform=proj, zorder=6)

# Europe inset
ia = fig.add_axes([0.80, 0.585, 0.16, 0.28], projection=proj)
ia.set_extent([-11, 35, 34, 71], crs=proj)
ia.add_feature(cfeature.LAND.with_scale('50m'), facecolor='#d9d9d6', edgecolor='none')
ia.add_feature(cfeature.OCEAN.with_scale('50m'), facecolor='white')
ia.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor='#666666', linewidth=0.4)
ia.add_feature(cfeature.COASTLINE.with_scale('50m'), edgecolor='#444444', linewidth=0.4)
ia.plot([EXT[0], EXT[1], EXT[1], EXT[0], EXT[0]], [EXT[2], EXT[2], EXT[3], EXT[3], EXT[2]], color='k', lw=1.4, transform=proj)
for sp in ia.spines.values():
    sp.set_edgecolor('#333333')

handles = [Line2D([], [], color=C_CML, lw=2.2, label='CML'),
           Patch(facecolor=C_RADAR, edgecolor='#9ec3df', label='Radar coverage'),
           Line2D([], [], marker='^', ls='none', color=C_GAUGE, markersize=9, label='Rain gauge'),
           Line2D([], [], marker='o', ls='none', color=C_STATION, markersize=9, label='Weather station')]
fig.legend(handles=handles, loc='upper center', ncol=4, frameon=False, fontsize=14, bbox_to_anchor=(0.5, 0.995))
fig.savefig(a.out, dpi=300, bbox_inches='tight', pad_inches=0.05)
fig.savefig(str(Path(a.out).with_suffix('.pdf')), bbox_inches='tight', pad_inches=0.05)
print('CMLs', len(segs), 'gauges', len(glon), 'weather stations (union T/RH)', len(st), 'radar valid cells', int(valid.sum()))
