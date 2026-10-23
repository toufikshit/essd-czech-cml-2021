"""Reader-friendly terminal commands; use `python run.py --help`."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from .config import load_paths,make_output_dir
from .catalog import inspect_all
from .compare import compare
from .workflows import availability,rebuild,validation,export_native_1min


def build_parser():
    p=argparse.ArgumentParser(description='Walk through and reproduce the ESSD Czech 2021 dataset workflows',
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--config',type=Path,required=True,help='JSON with seven explicit NetCDF paths')
    p.add_argument('--output-dir',type=Path,default=Path('results'),help='All generated data/figures go here')
    subs=p.add_subparsers(dest='task',required=True)
    subs.add_parser('check',help='Check seven NetCDF roles and Apr--Oct time coverage (metadata only)')
    a=subs.add_parser('availability',help='Reproduce five-sensor availability PNG/PDF/CSV')
    a.add_argument('--dpi',type=int,default=300)
    c=subs.add_parser('compare-rainfall',help='Full-history one-link recomputation vs published rain')
    c.add_argument('--link-index',type=int,action='append',help='QC CML index (repeatable; default 3)')
    c.add_argument('--start',default='2021-04-01',help='Start of the comparison period')
    c.add_argument('--days',type=int,default=14)
    c.add_argument('--threads',type=int,default=2)
    c.add_argument('--atol',type=float,default=1e-4)
    c.add_argument('--rtol',type=float,default=1e-4)
    c.add_argument('--no-plot',action='store_true')
    d5=subs.add_parser('derive-5min',help='Materialize the analysis 5-minute rainfall from the released native 1-minute product')
    d5.add_argument('--link-batch',type=int,default=128,help='Links per read batch (memory control)')
    rr=subs.add_parser('rebuild-all-rainfall',help='Very expensive 2219-link legacy 5-min rerun; needs explicit consent')
    rr.add_argument('--workers',type=int,default=4)
    rr.add_argument('--i-understand-full-source-cost',action='store_true',help='Allow 2020--2021 source to be processed; no original input modified')
    ex=subs.add_parser('export-1min-rainfall',help='Author step: regenerate the published native 1-minute rainfall from FULL-HISTORY sources with mandatory 5-min comparison')
    ex.add_argument('--workers',type=int,default=8)
    v=subs.add_parser('validate-gauge',help='Paper gauge plots (original detector required for exact replay)')
    v.add_argument('--detector',type=Path,help='Exact historical helper, if recovered')
    v.add_argument('--allow-unverified-v3',action='store_true',help='Use explicitly unverified experimental v3 helper')
    v.add_argument('--include-event',action='store_true')
    v.add_argument('--only-check-sources',action='store_true')
    return p


def main(argv=None):
    a=build_parser().parse_args(argv)
    paths=load_paths(a.config)
    if a.task=='check':
        results=inspect_all(paths)
        for item in results:
            print(f"{item['role']:12} {item.get('selected_slots','?'):>7} native slots "
                  f"entities={item.get('entities','?')}, "
                  f"physically cropped={item.get('is_physically_cropped','?')}: "
                  f"{'PASS' if item['ok'] else 'FAIL'}")
            for issue in item['problems']:print('    ERROR:',issue)
        print('NOTE: time-axis coverage is not a numerical-data identity test.')
        return 0 if all(row['ok'] for row in results) else 1
    target=make_output_dir(a.output_dir,paths)
    if a.task=='availability':availability(paths,target/'availability',dpi=a.dpi)
    elif a.task=='compare-rainfall':
        report=compare(paths,target/'rainfall_comparison',link_indices=a.link_index or [3],
                       start=a.start,days=a.days,atol=a.atol,rtol=a.rtol,
                       numba_threads=a.threads,make_plot=not a.no_plot)
        print('COMPARISON RESULT:', 'PASS (within numeric tolerance)' if report['all_pass'] else 'DIFFERENCES DETECTED')
        return 0 if report['all_pass'] else 1
    elif a.task=='derive-5min':
        from .rain5 import derive_5min_file,derived_5min_name,is_native_1min
        if not is_native_1min(paths['rain']):
            raise RuntimeError('Configured rain role is already a 5-minute product; nothing to derive.')
        dest=target/'derived_5min'/derived_5min_name(paths['rain'].name)
        derive_5min_file(paths['rain'],dest,link_batch=a.link_batch)
        print('DERIVED:',dest)
    elif a.task=='export-1min-rainfall':
        out=export_native_1min(paths,target/'native_1min_export',workers=a.workers)
        print('EXPORTED:',out)
    elif a.task=='rebuild-all-rainfall':
        from netCDF4 import Dataset
        with Dataset(paths['qc']) as ds:nt=len(ds.dimensions['time'])
        if nt!=308160 and not a.i_understand_full_source_cost:
            raise RuntimeError(f'QC source has {nt:,} minutes, not 308,160. Rerun would be full source history. '
                               'Pass --i-understand-full-source-cost ONLY if that is intended.')
        rebuild(paths,target/'full_rainfall',workers=a.workers)
    elif a.task=='validate-gauge':
        validation(paths,target,a.detector,allow_unverified_v3=a.allow_unverified_v3,
                   skip_event=not a.include_event,only_check_sources=a.only_check_sources)
    return 0


if __name__=='__main__':
    try:sys.exit(main())
    except (RuntimeError,ValueError,FileNotFoundError,KeyError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        sys.exit(2)
