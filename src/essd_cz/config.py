"""Read the data paths of config/user_paths.json. Every script and notebook uses this file."""
from __future__ import annotations
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / 'config' / 'user_paths.json'
POPULATION_DIR = ROOT / 'config' / 'population'
ORIGINAL = ROOT / 'original' / 'code'

ROLES = ('raw', 'qc', 'rain', 'gauge', 'radar', 'temperature', 'humidity')
# The seven released files (1 April to 31 October 2021).
EXPECTED_NAMES = {
    'raw': 'CZ_CML_RAW_1MIN_20210401_20211031.nc',
    'qc': 'CZ_CML_PREPROCESSED_1MIN_20210401_20211031.nc',
    'rain': 'CZ_CML_RAINRATE_GRAF_ITU2005_PASTOREK_1MIN_20210401_20211031.nc',
    'gauge': 'CZ_RAIN_GAUGE_10MIN_20210401_20211031.nc',
    'radar': 'CZ_RADAR_CMLPATH_5MIN_20210401_20211031.nc',
    'temperature': 'CHMI_10MIN_TEMPERATURE_20210401_20211031.nc',
    'humidity': 'CHMI_10MIN_HUMIDITY_20210401_20211031.nc',
}
# Optional inputs: the three WVD supporting files of the deposit and a folder with
# daily RADKLIM YW files (DWD), used for the radar processing and Figs. 1 and 7.
EXTRAS = ('wvd_v6_nc', 'wvd_v6_diagnostics_csv', 'wvd_quantization_csv', 'radklim_dir')


def _resolve(value, base: Path) -> Path:
    p = Path(value).expanduser()
    return (p if p.is_absolute() else base / p).resolve()


def load_config(config_file: str | Path | None = None) -> dict:
    """Return {'datasets': {role: Path}, 'extras': {name: Path or None}, 'output_dir': Path}.

    Relative paths are relative to the folder of the configuration file.
    """
    config_file = Path(config_file or DEFAULT_CONFIG).expanduser().resolve()
    if not config_file.is_file():
        raise FileNotFoundError(f'{config_file} not found. Copy config/user_paths.example.json to '
                                'config/user_paths.json and set your data paths.')
    settings = json.loads(config_file.read_text(encoding='utf-8'))
    base = config_file.parent
    provided = settings.get('datasets', {})
    missing = [r for r in ROLES if r not in provided]
    if missing:
        raise ValueError(f'Missing dataset path(s) in {config_file}: {missing}')
    datasets = {}
    for role in ROLES:
        path = _resolve(provided[role], base)
        if path.name != EXPECTED_NAMES[role]:
            raise ValueError(f'{role}: wrong file name {path.name!r}; expected {EXPECTED_NAMES[role]!r}')
        datasets[role] = path
    extras = {k: (_resolve(v, base) if v else None) for k, v in settings.get('extras', {}).items()}
    for k in EXTRAS:
        extras.setdefault(k, None)
    output_dir = _resolve(settings.get('output_dir', '../results'), base)
    return {'datasets': datasets, 'extras': extras, 'output_dir': output_dir}


def load_paths(config_file: str | Path | None = None) -> dict[str, Path]:
    """Paths of the seven released files; each must exist."""
    datasets = load_config(config_file)['datasets']
    for role, path in datasets.items():
        if not path.is_file():
            raise FileNotFoundError(f'{role}: {path}')
    return datasets


def output_dir(config_file: str | Path | None = None, *parts: str) -> Path:
    """Folder for generated files: output_dir of the configuration plus parts (created)."""
    out = load_config(config_file)['output_dir'].joinpath(*parts)
    out.mkdir(parents=True, exist_ok=True)
    return out


def extra(config_file: str | Path | None, name: str) -> Path:
    """Path of an optional input; stops with a clear message if it is not set."""
    value = load_config(config_file)['extras'].get(name)
    if value is None:
        raise FileNotFoundError(f'Set extras.{name} in the configuration file to run this step.')
    if not value.exists():
        raise FileNotFoundError(f'extras.{name}: {value} not found')
    return value
