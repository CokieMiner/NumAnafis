"""Retained-case plots and observed input-bin summaries."""
import argparse
import math
import os
import tempfile
import random
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'num-anafis-matplotlib'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from .inputs import DB_PATH, FUNCTIONS
from .reporting import LABELS, run_metadata
from .storage import Store

CHARTS_DIR = DB_PATH.parent/'charts'
DEFAULT_BINS = 80
DEFAULT_MAX_SLICES = 12
DEFAULT_DPI = 160
ULP_FLOOR = 1e-3
ULP_CAP = 1e12
CLASS_COLORS = {'rounded': '#227c9d', 'within_1': '#17a398', 'within_5': '#e9c46a',
                'within_10': '#f4a261', 'severe': '#d62828', 'signed_zero': '#9d0208'}


def run_viz_cli(argv=None):
    parser = argparse.ArgumentParser(description='Plot retained precision cases from SQLite')
    parser.add_argument('--db', type=Path, default=DB_PATH)
    parser.add_argument('--run')
    parser.add_argument('--backend', choices=('32', '64'))
    parser.add_argument('--function')
    parser.add_argument('--output', type=Path, default=CHARTS_DIR)
    parser.add_argument('--report', action='store_true', help='print observed-bin summaries instead of plotting')
    parser.add_argument('--max-points', type=int, default=100000, help='per-function plot limit; counts remain in report')
    parser.add_argument('--bins', type=int, default=DEFAULT_BINS)
    parser.add_argument('--max-slices', type=int, default=DEFAULT_MAX_SLICES)
    parser.add_argument('--dpi', type=int, default=DEFAULT_DPI)
    parser.add_argument('--xscale', choices=('linear', 'log', 'symlog'), default='symlog')
    args = parser.parse_args(argv)
    if not args.db.is_file() or args.max_points < 1 or args.bins < 2 or args.max_slices < 1 or args.dpi < 1:
        parser.error('database must exist and plot sizes must be positive (at least two bins)')
    store = Store(args.db, readonly=True)
    grouped = defaultdict(lambda: defaultdict(list))
    try:
        run_id, _, _ = run_metadata(store, args.run)
        reservoirs = {}
        for entry in store.entries(run_id, args.backend, args.function):
            key = (entry['backend'], entry['function'])
            seen, rng = reservoirs.setdefault(key, [0, random.Random(str(key))])
            rows = grouped[key[0]][key[1]]
            entry = enrich(entry)
            if len(rows) < args.max_points:
                rows.append(entry)
            else:
                index = rng.randrange(seen+1)
                if index < args.max_points:
                    rows[index] = entry
            reservoirs[key][0] += 1
        for key, (seen, _) in reservoirs.items():
            if seen > args.max_points:
                print(f'f{key[0]}/{key[1]}: plotting {args.max_points} of {seen} retained cases (deterministic reservoir)')
    finally:
        store.close()
    if not grouped:
        parser.error('no retained cases match the selection')
    if args.report:
        print_domain_report(grouped, args.bins)
        return 0
    for width, functions in sorted(grouped.items()):
        out_dir = args.output/width
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, entries in sorted(functions.items()):
            output = out_dir/f'{name}.png'
            if plot_function(_normalize_entries(entries, name), width, name, output, args.max_slices, args.xscale, args.dpi):
                print(output)
            else:
                print(f'f{width}/{name}: no finite measured errors to plot; inspect the report')
    return 0


def plot_function(entries, width, name, output, max_slices=12, xscale='symlog', dpi=DEFAULT_DPI):
    discrete, continuous = _classify_args(name)
    measured = [entry for entry in entries if entry.get('_ulp') is not None
                and math.isfinite(entry['our_value']) and math.isfinite(entry['x'])
                and all(math.isfinite(entry[key]) for key in continuous)]
    if not measured:
        return False
    groups = _group_by_discrete(measured, discrete) if discrete else {(): measured}
    keys = sorted(groups, key=lambda key: max(entry['_ulp'] for entry in groups[key]), reverse=True)[:max_slices]
    cols = min(3, len(keys))
    rows = math.ceil(len(keys)/cols)
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 3.5*rows), squeeze=False)
    for axis, key in zip(axes.flat, keys):
        data = groups[key]
        transform, label, axis_scale = plot_coordinates([entry['x'] for entry in data], xscale)
        data = [dict(entry, x=transform(entry['x'])) for entry in data]
        axis.set_xscale(axis_scale)
        low, high = min(entry['x'] for entry in data), max(entry['x'] for entry in data)
        if low != high:
            axis.set_xlim(low, high)
        if len(continuous) > 1:
            y_key = continuous[1]
            y_factor = max(abs(entry[y_key]) for entry in data)
            y_factor = y_factor if y_factor > 1e250 else 1.
            image = axis.scatter([entry['x'] for entry in data], [entry[y_key]/y_factor for entry in data],
                                 c=[_safe_ulp(entry['_ulp']) for entry in data], cmap='viridis',
                                 norm=LogNorm(ULP_FLOOR, ULP_CAP), s=10, rasterized=True)
            fig.colorbar(image, ax=axis, label='ULP error (display clipped)')
            axis.set_ylabel(y_key if y_factor == 1 else f'{y_key} / {y_factor:.6g}')
        else:
            for category in dict.fromkeys(entry['_class'] for entry in data):
                selected = [entry for entry in data if entry['_class'] == category]
                axis.scatter([entry['x'] for entry in selected],
                             [_safe_ulp(entry['_ulp']) for entry in selected], s=8,
                             color=CLASS_COLORS.get(category, 'black'), label=LABELS.get(category, category), rasterized=True)
            axis.axhline(1, color='grey', linestyle='--', linewidth=.7)
            axis.axhline(10, color='grey', linestyle=':', linewidth=.7)
            axis.set_yscale('log')
            axis.set_ylabel('ULP error (log scale)')
            axis.legend(fontsize=7)
        axis.set_xlabel(label)
        axis.set_title(', '.join(f'{label}={value}' for label, value in zip(discrete, key)) or 'Measured inputs')
        axis.grid(alpha=.2)
    for axis in list(axes.flat)[len(keys):]:
        axis.set_visible(False)
    omitted = len(entries)-len(measured)
    fig.suptitle(f'{name} / f{width}\nRetained cases; {len(keys)}/{len(groups)} parameter slices')
    caption = (f'{len(entries)} retained cases; {omitted} exceptional/unscored cases omitted.\n'
               f'Zero errors displayed at {ULP_FLOOR:g}; ULP display clipped at {ULP_CAP:g}.')
    fig.text(.02, .01, caption, fontsize=7)
    fig.tight_layout(rect=(0,.06,1,.94))
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return True


def print_domain_report(grouped, bins=80):
    print('Observed input bins. Class counts apply only to retained tested points; no interval guarantee.')
    for width, functions in sorted(grouped.items()):
        for name, entries in sorted(functions.items()):
            print(f'\nf{width}/{name}:')
            for bucket in _bin_entries(entries, bins):
                values = [entry['_ulp'] for entry in bucket['entries'] if entry['_ulp'] is not None]
                maximum = max(values, default=None)
                counts = dict(Counter(entry['_class'] for entry in bucket['entries']))
                print(f'  [{bucket["start"]:.8g}, {bucket["end"]:.8g}]: max={maximum} ULP; counts={counts}')


def _normalize_entries(entries, name):
    normalized = []
    spec = FUNCTIONS[name]
    for raw in entries:
        entry = dict(raw)
        entry['x'] = entry['input']
        for index, (value, (_, _, integer)) in enumerate(zip(entry['extras'], spec[1:])):
            entry[f'p{index+1}'] = value
        normalized.append(entry)
    return normalized


def _classify_args(name):
    discrete, continuous = [], ['x']
    for index, (_, _, integer) in enumerate(FUNCTIONS[name][1:]):
        (discrete if integer else continuous).append(f'p{index+1}')
    return discrete, continuous


def _group_by_discrete(entries, keys):
    groups = defaultdict(list)
    for entry in entries:
        groups[tuple(entry[key] for key in keys)].append(entry)
    return dict(groups)


def _bin_entries(entries, bins):
    finite = [entry for entry in entries if math.isfinite(entry['input'])]
    if not finite:
        return []
    lo = min(entry['input'] for entry in finite)
    hi = max(entry['input'] for entry in finite)
    if lo == hi:
        return [{'start': lo, 'end': hi, 'entries': finite}]
    # Use observed extrema; do not relabel out-of-domain samples as boundary cases.
    weights = np.linspace(0., 1., max(2, bins)+1)
    edges = (1-weights)*lo+weights*hi
    edges = np.unique(edges)
    buckets = [[] for _ in range(len(edges)-1)]
    for entry in finite:
        index = min(len(buckets)-1, max(0, int(np.searchsorted(edges, entry['input'], side='right'))-1))
        buckets[index].append(entry)
    return [{'start': float(edges[i]), 'end': float(edges[i+1]), 'entries': rows}
            for i, rows in enumerate(buckets) if rows]


def _safe_ulp(value):
    if value is None or not math.isfinite(value):
        return ULP_CAP
    return min(ULP_CAP, max(ULP_FLOOR, value))


def enrich(entry):
    entry = dict(entry)
    entry['_class'] = entry['classification']
    entry['_ulp'] = entry.get('ulp_error')
    return entry


def plot_coordinates(values, scale):
    """Keep plotted coordinates finite across the full represented input range."""
    if scale == 'log' and any(value <= 0 for value in values):
        scale = 'symlog'
    largest = max(abs(value) for value in values)
    if largest > 1e250 and scale == 'symlog':
        return (lambda x: math.copysign(math.log1p(abs(x))/math.log(10), x),
                'sign(x) log10(1 + |x|)', 'linear')
    if largest > 1e250 and scale == 'linear':
        return (lambda x: x/largest, f'x / {largest:.6g}', 'linear')
    return (lambda x: x, 'x', scale)
