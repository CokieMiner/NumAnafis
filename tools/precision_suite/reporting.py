"""Full-coverage accuracy reports and streaming exports from SQLite."""
import argparse
from pathlib import Path

from .inputs import DB_PATH, FUNCTIONS
from .storage import Store, dumps, loads

LABELS = {
    'native_error': 'Native evaluation failed (unscored)', 'rounded': 'Matches rounded reference', 'within_1': 'Error ≤ 1 ULP',
    'within_5': 'Error ≤ 5 ULP', 'within_10': 'Error ≤ 10 ULP',
    'severe': 'Error > 10 ULP', 'mismatch': 'Exceptional-value mismatch',
    'expected_nonfinite': 'Expected NaN/infinity', 'domain': 'Real-domain rejection', 'signed_zero': 'Signed-zero mismatch',
    'pole': 'Pole (unscored)', 'reference_error': 'Reference failed (unscored)',
    'reference_unstable': 'Reference unstable (unscored)',
}


def run_report_cli(argv=None):
    parser = argparse.ArgumentParser(description='Summarize precision runs from SQLite')
    parser.add_argument('--db', type=Path, default=DB_PATH)
    parser.add_argument('--run', help='run ID; default: latest run')
    parser.add_argument('--backend', choices=('32', '64'))
    parser.add_argument('--function')
    parser.add_argument('--worst', type=int, default=5)
    parser.add_argument('--only-severe', action='store_true')
    args = parser.parse_args(argv)
    if args.worst < 0 or args.function is not None and args.function not in FUNCTIONS:
        parser.error('worst must be nonnegative and function must be supported')
    if not args.db.is_file():
        parser.error(f'database not found: {args.db}')
    store = Store(args.db, readonly=True)
    try:
        run_id, configuration, completed = run_metadata(store, args.run)
        print(f'Run {run_id}: {"complete" if completed else "incomplete"}; mode={configuration["mode"]}; seed={configuration["seed"]}')
        print(f'Reference: mpmath {configuration["mpmath_version"]}, adaptive precision; agreement is not certification.')
        print_error_definition()
        if configuration['mode'] == 'exhaustive':
            print('Coverage is the configured bit range and fixed extra arguments; plot cases are retained samples.')
        for width, name, cursor, summary in store.summaries(run_id):
            if args.backend and width != int(args.backend) or args.function and name != args.function:
                continue
            scored = sum(summary['counts'].get(key, 0) for key in ('rounded', 'within_1', 'within_5', 'within_10', 'severe', 'signed_zero', 'mismatch'))
            maximum = f'{summary["max_ulp_hp"]} ULP' if scored else 'unscored'
            print(f'\nf{width}/{name}: {summary["count"]} tested; max error={maximum}')
            for category, count in sorted(summary['counts'].items()):
                if args.only_severe and category not in ('native_error', 'severe', 'mismatch', 'signed_zero', 'reference_error', 'reference_unstable'):
                    continue
                print(f'  {LABELS.get(category, category):36} {count:10d} {100*count/max(1,summary["count"]):7.3f}%')
            for record in summary.get('failures', [])[:max(0, args.worst)]:
                print(f'  {record["classification"]}: input=0x{record["input_bits"]} output={record["output_bits"]} '
                      f'x={record["input"]!r} extras={record["extras"]!r}; {record.get("reference_message") or record.get("reference_hp")}')
            for record in summary['worst'][:max(0, args.worst)]:
                if args.only_severe and record['classification'] not in ('severe', 'mismatch'):
                    continue
                print(f'  x={record["input"]!r} extras={record["extras"]!r} bits=0x{record["input_bits"]} '
                      f'ours={record["our_value"]!r} reference={record.get("reference_hp",record.get("reference"))} '
                      f'error={record.get("ulp_error_hp", record["ulp_error"])} ULP; distance={record.get("ulp_distance")}; '
                      f'absolute={record.get("abs_error_hp", record.get("abs_error"))}; relative={record.get("rel_error_hp", record.get("rel_error"))}')
    finally:
        store.close()
    return 0


def run_export_cli(argv=None):
    parser = argparse.ArgumentParser(description='Export retained precision cases from SQLite')
    parser.add_argument('--db', type=Path, default=DB_PATH)
    parser.add_argument('--run')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    if not args.db.is_file() or args.db.resolve() == args.output.resolve():
        parser.error('database must exist; export must use a different path')
    store = Store(args.db, readonly=True)
    try:
        run_id, configuration, completed = run_metadata(store, args.run)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('w') as stream:
            stream.write('{"run_id":'+dumps(run_id)+',"configuration":'+dumps(configuration)+',"completed":'+dumps(completed))
            stream.write(',"summaries":'+dumps(store.summaries(run_id))+',"retained_cases":[')
            separator = ''
            for entry in store.entries(run_id):
                stream.write(separator+dumps(entry))
                separator = ','
            stream.write(']}\n')
    finally:
        store.close()
    return 0


def run_metadata(store, run_id=None):
    """Validate a run before accessing its retained rows or coverage summaries."""
    run_id = run_id or store.latest()
    row = store.db.execute('SELECT configuration,completed FROM runs WHERE id=?', (run_id,)).fetchone()
    if row is None:
        raise ValueError(f'unknown run {run_id}')
    return run_id, loads(row[0]), bool(row[1])


def print_error_definition():
    """The denominator follows the reference binade, with fixed subnormal spacing."""
    print('ULP error = |output - high-precision reference| / target spacing at the reference.')
    print('Distance counts representable steps to the rounded reference. Absolute/relative errors are separate.')
    print('Signed zero is checked separately; no absolute EPSILON threshold is applied.')
