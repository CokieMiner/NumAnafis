"""Persistence, exhaustive coverage, command failures, reports, and plots."""
import contextlib
import io
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.precision_suite.__main__ import main
from tools.precision_suite.reporting import run_export_cli, run_report_cli
from tools.precision_suite.storage import Store, dumps, loads

ROOT = Path(__file__).resolve().parents[3]


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.path = Path(self.directory.name)
        self.db = self.path/'results.sqlite3'
        self.store = Store(self.db)
        self.configuration = {'mode': 'exhaustive', 'seed': 0, 'mpmath_version': 'test', 'extras': 'nan'}
        self.store.start('test', self.configuration, 'today')
        self.closed = False

    def tearDown(self):
        if not self.closed:
            self.store.close()
        self.directory.cleanup()

    def test_nonfinite_encoding_preserves_arbitrary_strings(self):
        data = {'label': 'nan', 'values': [math.nan, math.inf, -math.inf, -0.]}
        encoded = dumps(data)
        json.loads(encoded, parse_constant=lambda value: self.fail(value))
        decoded = loads(encoded)
        self.assertEqual(decoded['label'], 'nan')
        self.assertTrue(math.isnan(decoded['values'][0]))
        self.assertEqual(decoded['values'][1:3], [math.inf, -math.inf])
        self.assertEqual(math.copysign(1., decoded['values'][3]), -1.)

    def test_exact_resume_configuration_and_completion(self):
        self.assertEqual(self.store.start('unused', self.configuration, 'later', 'test'), 'test')
        with self.assertRaises(ValueError):
            self.store.start('unused', {'changed': True}, 'later', 'test')
        with self.assertRaises(ValueError):
            self.store.start('unused', self.configuration, 'later', 'unknown')
        self.store.complete('test')
        with self.assertRaises(ValueError):
            self.store.start('unused', self.configuration, 'later', 'test')

    def test_checkpoint_counts_and_retention_are_independent(self):
        records = [self.record(i, i) for i in range(1000)]
        summary = self.store.commit_batch('test', 32, 'sin', 1000, records, True)
        self.assertEqual(summary['count'], 1000)
        self.assertEqual(sum(summary['counts'].values()), 1000)
        self.assertEqual(summary['max_ulp'], 999)
        self.assertEqual(len(list(self.store.entries('test'))), 1)
        self.assertEqual(list(self.store.entries('test'))[0]['ulp_error'], 999)
        self.assertEqual(self.store.progress('test', 32, 'sin')[0], 1000)
        self.store.commit_batch('test', 32, 'sin', 1001, [self.record(1000, 1)], True)
        self.assertEqual(list(self.store.entries('test'))[0]['ulp_error'], 999)

    def test_worst_ranking_preserves_errors_above_float_range(self):
        first, second = self.record(1, math.inf), self.record(2, math.inf)
        first['ulp_error_hp'], second['ulp_error_hp'] = '1e400', '1e500'
        summary = self.store.commit_batch('test', 32, 'sin', 2, [first, second], True)
        self.assertEqual(summary['worst'][0]['input_bits'], second['input_bits'])
        self.assertEqual(summary['max_ulp_hp'], '1e500')
        self.assertEqual(list(self.store.entries('test'))[0]['input_bits'], second['input_bits'])

    def test_ranking_handles_extreme_negative_decimal_exponents(self):
        record = self.record(1, 0.)
        record['ulp_error_hp'] = '1.3e-73786976349690225781'
        self.store.commit_batch('test', 32, 'sin', 1, [record], False)
        self.assertEqual(self.store.progress('test', 32, 'sin')[1]['max_ulp_hp'], record['ulp_error_hp'])

    def test_failed_commit_rolls_back_all_rows(self):
        malformed = self.record(2, 1)
        malformed['unserializable'] = object()
        with self.assertRaises(TypeError):
            self.store.commit_batch('test', 32, 'sin', 2, [self.record(1, 1), malformed], False)
        self.assertEqual(self.store.progress('test', 32, 'sin')[0], 0)
        self.assertEqual(list(self.store.entries('test')), [])
        with self.assertRaises(ValueError):
            self.store.commit_batch('test', 32, 'sin', 10, [self.record(1, 1)], False)

    def test_sampled_runs_keep_every_row_and_failures(self):
        failure = self.record(3, None)
        failure.update(classification='reference_error', reference_message='test failure')
        summary = self.store.commit_batch('test', 32, 'sin', 3,
                                           [self.record(1, 1), self.record(2, 2), failure], False)
        self.assertEqual(len(list(self.store.entries('test'))), 3)
        self.assertEqual(summary['counts']['reference_error'], 1)
        self.assertEqual(summary['failures'], [failure])

    def test_report_and_strict_export(self):
        self.store.commit_batch('test', 32, 'sin', 1, [self.record(1, 11)], True)
        with contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertEqual(run_report_cli(['--db', str(self.db)]), 0)
        self.assertIn('not certification', output.getvalue())
        exported = self.path/'export.json'
        self.assertEqual(run_export_cli(['--db', str(self.db), '--output', str(exported)]), 0)
        data = json.loads(exported.read_text())
        self.assertEqual(data['summaries'][0][3]['count'], 1)
        self.assertEqual(len(data['retained_cases']), 1)
        self.assertEqual(data['configuration']['extras'], 'nan')

    def test_plotting_and_observed_bounds(self):
        from tools.precision_suite.plotting import _bin_entries, run_viz_cli
        self.store.commit_batch('test', 32, 'sin', 2,
                                [self.record(1, 0, -100.), self.record(2, 1, 100.)], False)
        bounds = _bin_entries(list(self.store.entries('test')), 2)
        self.assertEqual(bounds[0]['start'], -100.)
        self.assertEqual(bounds[-1]['end'], 100.)
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(run_viz_cli(['--db', str(self.db), '--output', str(self.path/'charts')]), 0)
        self.assertGreater((self.path/'charts/32/sin.png').stat().st_size, 1000)

    def test_plots_and_bins_accept_extreme_f64_inputs(self):
        from tools.precision_suite.plotting import _bin_entries, _normalize_entries, enrich, plot_function
        maximum = sys.float_info.max
        records = [self.record(1, 1, -maximum), self.record(2, 2, maximum)]
        bounds = _bin_entries(records, 2)
        self.assertTrue(all(math.isfinite(item['start']) and math.isfinite(item['end']) for item in bounds))
        entries = _normalize_entries([enrich(record) for record in records], 'sin')
        self.assertTrue(plot_function(entries, 64, 'sin', self.path/'wide.png'))

    def test_cli_rejects_invalid_inputs_before_build(self):
        for arguments in (['--function', 'unknown'], ['--backend', 'both', '--mode', 'exhaustive'],
                          ['--fail-ulp', 'nan'], ['--timeout', 'nan'], ['--reference-timeout', '-1']):
            with patch('tools.precision_suite.__main__.build_runner') as build:
                with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                    main(['run', *arguments])
                self.assertEqual(error.exception.code, 2)
                build.assert_not_called()

    def test_exhaustive_resume_has_no_duplicates(self):
        self.store.close()
        self.closed = True
        scan_db = self.path/'scan.sqlite3'
        options = ['--function', 'sin', '--backend', '32', '--mode', 'exhaustive',
                   '--start-bits', '0', '--end-bits', '64', '--batch-size', '16', '--db', str(scan_db)]
        first = self.cli('run', *options, '--max-batches', '1')
        self.assertEqual(first.returncode, 0, first.stderr)
        store = Store(scan_db, readonly=True)
        run_id = store.latest()
        self.assertEqual(store.progress(run_id, 32, 'sin')[0], 16)
        store.close()
        second = self.cli('run', *options, '--resume', run_id)
        self.assertEqual(second.returncode, 0, second.stderr)
        store = Store(scan_db, readonly=True)
        self.assertEqual(store.progress(run_id, 32, 'sin')[1]['count'], 64)
        self.assertEqual(store.db.execute('SELECT completed FROM runs WHERE id=?', (run_id,)).fetchone()[0], 1)
        store.close()

    def test_contiguous_shards_cover_the_range(self):
        self.store.close()
        self.closed = True
        tested = 0
        for shard in range(3):
            path = self.path/f'shard{shard}.sqlite3'
            result = self.cli('run', '--function', 'sin', '--backend', '32', '--mode', 'exhaustive',
                              '--start-bits', '0', '--end-bits', '10', '--shards', '3', '--shard', str(shard),
                              '--db', str(path))
            self.assertEqual(result.returncode, 0, result.stderr)
            store = Store(path, readonly=True)
            tested += store.summaries(store.latest())[0][3]['count']
            store.close()
        self.assertEqual(tested, 10)

    def record(self, bits, ulp, x=1.):
        return {'run_id': 'test', 'backend': '32', 'function': 'sin', 'input': x, 'extras': [],
                'input_bits': f'{bits:08x}', 'output_bits': '00000000', 'our_value': 0.,
                'reference': 0., 'reference_hp': '0', 'abs_error': 0., 'rel_error': None,
                'ulp_distance': 0, 'ulp_error': ulp, 'classification': 'severe' if ulp and ulp > 10 else 'rounded'}

    def cli(self, *arguments):
        return subprocess.run([sys.executable, '-m', 'tools.precision_suite', *arguments], cwd=ROOT,
                              capture_output=True, text=True, timeout=30)
