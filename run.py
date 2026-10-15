#!/usr/bin/env python3
"""Start here: `python run.py --help` then follow README.md."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from essd_cz.cli import main
if __name__=='__main__':
    try:sys.exit(main())
    except (RuntimeError,ValueError,FileNotFoundError,KeyError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        sys.exit(2)
