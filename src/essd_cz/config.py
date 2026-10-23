"""Resolve user-provided input paths. Never guess which rain product to use."""
from __future__ import annotations
from pathlib import Path
import json

ROLES = ('raw', 'qc', 'rain', 'gauge', 'radar', 'temperature', 'humidity')
# Public release files (April to October 2021), as published on Zenodo.
# The single published CML rainfall product is the NATIVE 1-MINUTE NetCDF; the
# analysis 5-minute series are derived from it (src/essd_cz/rain5.py).
EXPECTED_NAMES = {
    'raw': 'CZ_CML_RAW_1MIN_20210401_20211031.nc',
    'qc': 'CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc',
    'rain': 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc',
    'gauge': 'CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc',
    'radar': 'CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc',
    'temperature': 'CHMI_10MIN_TEMPERATURE_20210401_20211031.nc',
    'humidity': 'CHMI_10MIN_HUMIDITY_20210401_20211031.nc',
}
# The author's full history source files (2020 to 2021), used by the hpc
# profile. The full-history rain role stays the legacy 5-minute reference
# product; it is an internal regression benchmark, not a Zenodo file.
FULL_HISTORY_NAMES = {
    'raw': 'CZ_CML_RAW_1MIN_20200101_20211231.nc',
    'qc': 'CZ_CML_PREPROCESSED_1MIN_20200101_20211231.nc',
    'rain': 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_5MIN_20200101_20211231.nc',
    'gauge': 'CZ_RAIN_GAUGE_10MIN_20200101_20211231.nc',
    'radar': 'CZ_RADAR_CMLPATH_5MIN_20200101_20211231.nc',
    'temperature': 'CHMI_10MIN_TEMPERATURE_2021.nc',
    'humidity': 'CHMI_10MIN_HUMIDITY_2021.nc',
}
ACCEPTED_NAMES = {role: (EXPECTED_NAMES[role], FULL_HISTORY_NAMES[role]) for role in ROLES}


def load_paths(config_file: str | Path) -> dict[str, Path]:
    """Return absolute paths for the seven original files; never copy/modify them.

    Relative paths refer to the directory containing config_file, so a user
    may relocate the code repository and downloaded data together.
    """
    config_file = Path(config_file).expanduser().resolve()
    settings = json.loads(config_file.read_text(encoding='utf-8'))
    provided = settings.get('datasets', {})
    missing = [r for r in ROLES if r not in provided]
    if missing:
        raise ValueError(f'Missing dataset path(s) in {config_file}: {missing}')
    found = {}
    for role in ROLES:
        p = Path(provided[role]).expanduser()
        path = (p if p.is_absolute() else config_file.parent / p).resolve()
        if path.name not in ACCEPTED_NAMES[role]:
            raise ValueError(f'{role}: wrong filename {path.name!r}; expected one of {ACCEPTED_NAMES[role]}')
        if not path.is_file():
            raise FileNotFoundError(f'{role}: {path}')
        found[role] = path
    return found


def make_output_dir(path: str | Path, inputs: dict[str, Path]) -> Path:
    """Reject output directories that collide with inputs."""
    output = Path(path).expanduser().resolve()
    if output in {p for p in inputs.values()}:
        raise ValueError('Output path must not be an input NetCDF')
    output.mkdir(parents=True, exist_ok=True)
    return output
