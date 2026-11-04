"""Label-only text changes for the publication figures.

The figures are drawn by the unchanged producers in original/code/figures. To align their axis labels, legends and
colour-bar labels with the terminology of the paper, the exact text strings below are replaced in the finished figure
just before it is saved. Nothing else (data, limits, colours, sizes, fonts) is touched.

Use in a wrapper:   relabel.install(relabel.FIG4)
Use for a producer that runs as a separate process:
    python -m essd_cz.relabel FIG7 <producer.py> <producer arguments>
"""
import runpy
import sys

import matplotlib.figure
import matplotlib.text

FIG2 = {'Entities with valid data at native timestep (%)': 'CMLs, radar paths or stations with valid data [%]'}
FIG4 = {'TRSL [dB]': 'Total signal loss [dB]', 'CML (Graf)': 'CML', 'Gauge': 'Rain gauge',
        'Temperature [°C]': 'Air temperature [°C]', 'Humidity [%]': 'Relative humidity [%]',
        'Sublink 1': 'Sublink 0', 'Sublink 2': 'Sublink 1'}  # sublink numbers as sublink_id in the files
FIG5 = {'CML Graf/RSTD detection': 'CML rain detection'}
FIG6 = {r'RR$_{gauge}$ [mm h$^{-1}$]': r'Gauge rain rate [mm h$^{-1}$]',
        r'RR$_{CML}$ [mm h$^{-1}$]': r'CML rain rate [mm h$^{-1}$]'}
FIG7 = {'Hourly rain sum [mm]': '1 h rainfall accumulation [mm]'}

_savefig = matplotlib.figure.Figure.savefig


def install(mapping):
    """Replace exact text strings of every figure saved from now on."""
    def savefig(self, *args, **kwargs):
        for t in self.findobj(matplotlib.text.Text):
            if t.get_text() in mapping:
                t.set_text(mapping[t.get_text()])
        return _savefig(self, *args, **kwargs)
    matplotlib.figure.Figure.savefig = savefig


if __name__ == '__main__':
    install(globals()[sys.argv[1]])
    sys.argv = sys.argv[2:]
    runpy.run_path(sys.argv[0], run_name='__main__')
