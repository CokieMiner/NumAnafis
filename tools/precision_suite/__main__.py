"""Command-line orchestration for sampled and exhaustive primitive verification."""
import argparse
import hashlib
import itertools
import math
import platform
import sqlite3
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp

from .inputs import DB_PATH, FUNCTIONS, SAMPLES, TARGETED_SAMPLES, Case, finalize_args, from_bits, sample_cases
from .numerical import ReferenceEvaluator, build_runner, native_batch
from .storage import Store


def main(argv=None):
    arguments = sys.argv[1:] if argv is None else argv
    parser = argparse.ArgumentParser(description='Primitive f32/f64 accuracy verification, SQLite reports, and plots')
    parser.add_argument('command', choices=('run', 'report', 'plot', 'export'))
    if not arguments or arguments[0] in ('-h', '--help'):
        parser.print_help()
        return 0
    command = parser.parse_args(arguments[:1]).command
    try:
        if command == 'run':
            return run_verification_cli(arguments[1:])
        if command == 'report':
            from .reporting import run_report_cli
            return run_report_cli(arguments[1:])
        if command == 'plot':
            from .plotting import run_viz_cli
            return run_viz_cli(arguments[1:])
        from .reporting import run_export_cli
        return run_export_cli(arguments[1:])
    except (ValueError, RuntimeError, OSError, sqlite3.Error, subprocess.SubprocessError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 2


def run_verification_cli(argv=None):
    args, names, widths, extras = parse_run_options(argv)
    if args.list:
        for name, spec in sorted(FUNCTIONS.items()):
            print(f'{name}: '+', '.join('integer' if integer else 'real' for _, _, integer in spec))
        return 0
    runner = build_runner(args.no_std)
    configuration = {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()
                     if key not in ('db', 'resume', 'fail_ulp', 'list', 'max_batches')}
    configuration.update(runner_sha256=hashlib.sha256(Path(runner).read_bytes()).hexdigest(),
                         reference_sha256=source_hash(), mpmath_version=mp.__version__,
                         python_version=platform.python_version(), platform=platform.platform(),
                         machine=platform.machine(), libc=list(platform.libc_ver()))
    store = Store(args.db)
    try:
        run_id = store.start(uuid.uuid4().hex[:12], configuration, datetime.now(timezone.utc).isoformat(), args.resume)
    except Exception:
        store.close()
        raise
    reference = ReferenceEvaluator(args.reference_timeout)
    print(f'Run {run_id}; results: {args.db}', flush=True)
    failed = False
    committed = 0
    try:
        for width in widths:
            for name in names:
                cursor, summary = store.progress(run_id, width, name)
                if args.mode == 'exhaustive':
                    span = args.end_bits-args.start_bits
                    start = args.start_bits + span*args.shard//args.shards
                    end = args.start_bits + span*(args.shard+1)//args.shards
                    iterator = (Case(finalize_args(name, [from_bits(word, 32), *extras], 32), 'exhaustive', word)
                                for word in range(start+cursor, end))
                    total = end-start
                else:
                    iterator = itertools.islice(sample_cases(name, width, args.samples, args.targeted, args.seed, args.mode), cursor, None)
                    total = None
                started = time.monotonic()
                initial = cursor
                for chunk in batches(iterator, args.batch_size):
                    values = native_batch(runner, name, chunk, width, args.timeout)
                    records = [reference.evaluate(name, case.args, our, width, run_id, case.source,
                                              args.reference_bits, args.max_reference_bits, case.input_bits, output_bits, native_error)
                               for case, (our, output_bits, native_error) in zip(chunk, values)]
                    cursor += len(records)
                    summary = store.commit_batch(run_id, width, name, cursor, records, args.mode == 'exhaustive')
                    committed += 1
                    if args.max_batches is not None and committed >= args.max_batches:
                        print(f'Checkpoint saved; resume with --resume {run_id} and the same options', flush=True)
                        return 0
                    rate = (cursor-initial)/max(.001, time.monotonic()-started)
                    eta = f'; ETA {(total-cursor)/rate/3600:.2f}h' if total is not None else ''
                    print(f'  f{width}/{name}: {cursor} cases; {rate:.0f}/s; max {summary["max_ulp"]:.3g} ULP{eta}', flush=True)
                if args.fail_ulp is not None and summary['max_ulp'] > args.fail_ulp:
                    failed = True
                if any(summary['counts'].get(key) for key in ('native_error', 'mismatch', 'signed_zero', 'reference_error', 'reference_unstable')):
                    failed = True
        store.complete(run_id)
        return 1 if failed else 0
    except KeyboardInterrupt:
        print(f'Interrupted; resume with --resume {run_id} and the same options', file=sys.stderr)
        return 130
    finally:
        reference.close()
        store.close()


def parse_run_options(argv):
    parser = argparse.ArgumentParser(description='Compare F32Ext/F64Ext with high-precision real references')
    parser.add_argument('--backend', choices=('32', '64', 'both'), default='both')
    parser.add_argument('--function', help='comma-separated function names; default: all')
    parser.add_argument('--list', action='store_true', help='list supported functions without building')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--samples', type=positive, default=SAMPLES)
    parser.add_argument('--targeted', type=int, default=TARGETED_SAMPLES)
    parser.add_argument('--mode', choices=('random', 'grid', 'bits', 'exhaustive'), default='random')
    parser.add_argument('--start-bits', type=bit_index, default=0)
    parser.add_argument('--end-bits', type=bit_index, default=1 << 32, help='exclusive raw-bit endpoint')
    parser.add_argument('--shards', type=positive, default=1)
    parser.add_argument('--shard', type=int, default=0)
    parser.add_argument('--extras', type=str, help='comma-separated fixed parameters for exhaustive functions')
    parser.add_argument('--batch-size', type=positive, default=256)
    parser.add_argument('--reference-bits', type=positive, default=128)
    parser.add_argument('--max-reference-bits', type=positive, default=1024)
    parser.add_argument('--timeout', type=float, default=5, help='native batch and isolated-input timeout in seconds')
    parser.add_argument('--reference-timeout', type=float, default=5, help='per-case reference timeout in seconds')
    parser.add_argument('--db', type=Path, default=DB_PATH)
    parser.add_argument('--max-batches', type=positive, help='stop after this many committed batches; resume later')
    parser.add_argument('--resume', help='resume the given run ID with identical configuration')
    parser.add_argument('--no-std', action='store_true', help='evaluate libm paths rather than std paths')
    parser.add_argument('--fail-ulp', type=float, help='exit nonzero if measured error exceeds this ULP threshold')
    args = parser.parse_args(argv)
    if args.list:
        return args, [], [], []
    names = args.function.split(',') if args.function else sorted(FUNCTIONS)
    if len(names) != len(set(names)) or any(name not in FUNCTIONS for name in names):
        parser.error('function list contains duplicate or unsupported names; see --list')
    if (args.targeted < 0 or not math.isfinite(args.timeout) or args.timeout <= 0
            or args.reference_bits < 96 or args.max_reference_bits < 2*args.reference_bits):
        parser.error('invalid targeted count, timeout, or reference precision (minimum 96 bits, maximum at least twice initial)')
    if not 0 <= args.shard < args.shards or args.end_bits-args.start_bits < args.shards:
        parser.error('invalid shard or bit range')
    if not math.isfinite(args.reference_timeout) or args.reference_timeout <= 0:
        parser.error('reference-timeout must be finite and positive')
    if args.fail_ulp is not None and (not math.isfinite(args.fail_ulp) or args.fail_ulp < 0):
        parser.error('fail-ulp must be finite and nonnegative')
    widths = [32, 64] if args.backend == 'both' else [int(args.backend)]
    extras = [float(item) for item in args.extras.split(',')] if args.extras is not None else []
    if args.mode == 'exhaustive':
        if widths != [32]:
            parser.error('exhaustive evaluation requires --backend 32')
        for name in names:
            if len(FUNCTIONS[name])-1 != len(extras):
                parser.error(f'{name} requires {len(FUNCTIONS[name])-1} fixed --extras')
            for value, (_, _, integer) in zip(extras, FUNCTIONS[name][1:]):
                if integer and (not math.isfinite(value) or value != int(value) or not -(1 << 31) <= value < 1 << 31):
                    parser.error('integer extras must fit i32 exactly')
            if name in ('hermite', 'polygamma', 'zeta_deriv') and extras[0] < 0:
                parser.error(f'{name} requires a nonnegative order')
            if name in ('assoc_legendre', 'spherical_harmonic') and (extras[0] < 0 or abs(extras[1]) > extras[0]):
                parser.error(f'{name} requires degree >= 0 and |order| <= degree')
    elif (args.extras is not None or args.shards != 1 or args.shard
          or args.start_bits or args.end_bits != 1 << 32):
        parser.error('fixed extras, bit ranges, and sharding require exhaustive mode')
    return args, names, widths, extras


def source_hash():
    digest = hashlib.sha256()
    for path in sorted(Path(__file__).parent.glob('*.py')):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def batches(iterator, size):
    while chunk := list(itertools.islice(iterator, size)):
        yield chunk


def positive(value):
    number = int(value, 0)
    if number <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def bit_index(value):
    number = int(value, 0)
    if not 0 <= number <= 1 << 32:
        raise argparse.ArgumentTypeError('bit index must be between 0 and 2^32')
    return number


if __name__ == '__main__':
    sys.exit(main())
