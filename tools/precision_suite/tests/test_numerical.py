"""Rounding, reference contracts, and native-protocol regressions."""
import math
import random
import re
import sys
import tempfile
import unittest
from pathlib import Path

import mpmath as mp

from tools.precision_suite.inputs import Case, FUNCTIONS, float_bits, from_bits, round_input, sample_cases
from tools.precision_suite.numerical import (
    FORMATS, ReferenceEvaluator, build_runner, error_metrics, evaluate_reference,
    mathematical_reference, native_batch, ordered_bits, round_reference_bits,
)


class NumericalTests(unittest.TestCase):
    def test_rounding_matches_represented_values(self):
        rng = random.Random(71)
        for width in (32, 64):
            for _ in range(1000):
                bits = rng.getrandbits(width)
                value = from_bits(bits, width)
                if math.isfinite(value) and value != 0:
                    with mp.workprec(200):
                        self.assertEqual(round_reference_bits(mp.mpf(value), width), bits)

    def test_nearest_even_ties_and_binade_carry(self):
        with mp.workprec(200):
            for width, (p, emin, emax, _) in FORMATS.items():
                one_bits = float_bits(1., width)
                ulp = mp.ldexp(mp.mpf(1), 1-p)
                self.assertEqual(round_reference_bits(1+ulp/2, width), one_bits)
                self.assertEqual(round_reference_bits(1+3*ulp/2, width), one_bits+2)
                self.assertEqual(round_reference_bits(2-ulp/2, width), float_bits(2., width))
                tiny = mp.ldexp(mp.mpf(1), emin-p+1)
                self.assertEqual(round_reference_bits(tiny/2, width), 0)
                self.assertEqual(round_reference_bits(3*tiny/2, width), 2)
                maximum = mp.ldexp(mp.mpf(2)-ulp, emax)
                self.assertEqual(round_reference_bits(maximum, width), float_bits(math.inf, width)-1)
                self.assertEqual(round_reference_bits(maximum+mp.ldexp(ulp, emax-1), width), float_bits(math.inf, width))
                self.assertEqual(round_reference_bits(-tiny/2, width), 1 << (width-1))

    def test_f32_reference_avoids_double_rounding(self):
        with mp.workprec(200):
            reference = 1+mp.ldexp(mp.mpf(1), -24)+mp.ldexp(mp.mpf(1), -80)
            self.assertEqual(float_bits(float(reference), 32), float_bits(1., 32))
            self.assertEqual(round_reference_bits(reference, 32), float_bits(1., 32)+1)

    def test_error_uses_unrounded_reference(self):
        with mp.workprec(200):
            reference = 1+mp.ldexp(mp.mpf(1), -25)
            metric = error_metrics(1., reference, 32, 200)
            self.assertEqual(metric['ulp_error'], .25)
            self.assertEqual(metric['ulp_distance'], 0)
            self.assertTrue(metric['matches_rounded_reference'])
            self.assertIsInstance(metric['ulp_distance'], int)

    def test_binade_spacing_changes_at_power_of_two(self):
        with mp.workprec(200):
            for width, (p, _, _, _) in FORMATS.items():
                lower = mp.ldexp(mp.mpf(1), -p)
                below = from_bits(float_bits(1., width)-1, width)
                at_boundary = error_metrics(below, mp.mpf(1), width, 200)
                below_boundary = error_metrics(1., mp.mpf(1)-lower, width, 200)
                self.assertEqual(at_boundary['ulp_error'], .5)
                self.assertEqual(below_boundary['ulp_error'], 1.)
                self.assertEqual(at_boundary['ulp_distance'], 1)
                self.assertEqual(below_boundary['ulp_distance'], 1)

    def test_distance_is_monotone_across_zero(self):
        for width in (32, 64):
            values = [-math.inf, -2., -1., -from_bits(1, width), -0., 0., from_bits(1, width), 1., 2., math.inf]
            positions = [ordered_bits(float_bits(value, width), width) for value in values]
            self.assertEqual(positions, sorted(positions))
            self.assertEqual(positions[4], positions[5])
            self.assertEqual(positions[6]-positions[3], 2)

    def test_errors_too_small_for_python_float_remain_recorded(self):
        with mp.workprec(200):
            metric = error_metrics(0., mp.ldexp(mp.mpf(1), -2000), 64, 200)
            self.assertEqual(metric['abs_error'], 0.)
            self.assertGreater(mp.mpf(metric['abs_error_hp']), 0)
            self.assertGreater(mp.mpf(metric['ulp_error_hp']), 0)

    def test_tiny_absolute_error_is_not_excused(self):
        with mp.workprec(200):
            metric = error_metrics(0., mp.mpf('1e-20'), 32, 200)
        self.assertLess(metric['abs_error'], 2**-23)
        self.assertGreater(metric['ulp_error'], 1_000_000)
        self.assertEqual(metric['classification'], 'severe')

    def test_zero_spacing_and_signed_zero(self):
        with mp.workprec(200):
            self.assertEqual(error_metrics(from_bits(1, 32), mp.mpf(0), 32, 200)['ulp_error'], 1.)
            self.assertEqual(error_metrics(0., -0., 32, 200)['classification'], 'signed_zero')
            self.assertEqual(round_reference_bits(-0., 64), 1 << 63)
        self.assertEqual(error_metrics(math.nan, mp.inf, 32, 128)['classification'], 'mismatch')
        self.assertEqual(error_metrics(math.inf, mp.inf, 32, 128)['classification'], 'expected_nonfinite')

    def test_reference_domains_and_ieee_limits(self):
        for name, args in [('sqrt', [-1.]), ('elliptic_k', [2.]), ('lambert_wm1', [1.])]:
            reference, _, status, _ = evaluate_reference(name, args, 64)
            self.assertTrue(mp.isnan(reference))
            self.assertEqual(status, 'domain')
        self.assertEqual(evaluate_reference('ln', [0.], 64)[0], -math.inf)
        self.assertEqual(evaluate_reference('atanh', [1.], 64)[0], math.inf)
        self.assertEqual(evaluate_reference('gamma', [-2.], 64)[2], 'pole')
        self.assertEqual(evaluate_reference('powf', [-1., math.inf], 64)[0], 1.)
        self.assertEqual(float_bits(evaluate_reference('sin', [-0.], 32)[0], 32), 1 << 31)

    def test_function_conventions(self):
        with mp.workprec(200):
            self.assertEqual(mathematical_reference('elliptic_k', [.5]), mp.ellipk(mp.mpf('.5')))
            self.assertEqual(mathematical_reference('lgamma', [-.5]), mp.log(abs(mp.gamma(-.5))))
            self.assertEqual(mathematical_reference('assoc_legendre', [1., 3, 1]), 0)
            self.assertEqual(mathematical_reference('assoc_legendre', [.5, 1, -1]), mp.sqrt(mp.mpf('.75'))/2)
            theta, phi = mp.mpf(.75), mp.mpf(.5)
            expected = -mp.spherharm(2, 1, theta, phi).real
            self.assertEqual(mathematical_reference('spherical_harmonic', [.75, 2, -1, .5]), expected)
            self.assertAlmostEqual(float(mathematical_reference('zeta_deriv', [2., 1])), float(mp.diff(mp.zeta, 2)))

    def test_precision_agreement_is_required(self):
        self.assertEqual(evaluate_reference('sin', [1.], 64, 128, 128)[2], 'reference_unstable')
        self.assertEqual(evaluate_reference('sin', [1.], 64, 128, 256)[2], 'stable')

    def test_sampling_is_reproducible_and_target_rounded(self):
        first = list(sample_cases('sin', 32, 8, 2, 71, 'bits'))
        second = list(sample_cases('sin', 32, 8, 2, 71, 'bits'))
        self.assertEqual([case.input_bits for case in first], [case.input_bits for case in second])
        for case in first:
            if not math.isnan(case.args[0]):
                self.assertEqual(case.args[0], round_input(case.args[0], 32))
        roots = list(sample_cases('lgamma', 64, 1, 0, 71))
        for root in (1., 2.):
            for value in (math.nextafter(root, -math.inf), root, math.nextafter(root, math.inf)):
                self.assertIn(float_bits(value, 64), [case.input_bits for case in roots])


class NativeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner = build_runner()

    def test_catalog_covers_both_traits(self):
        root = Path(__file__).resolve().parents[3]
        for width in (32, 64):
            source = (root/f'src/ext/f{width}ext.rs').read_text().split(f'impl F{width}Ext')[0]
            self.assertEqual(set(re.findall(r'fn (\w+)\(', source)), set(FUNCTIONS))

    def test_native_width_and_integer_parameters(self):
        for width in (32, 64):
            case = Case([.5, 2], 'test', float_bits(.5, width))
            value, bits, error = native_batch(self.runner, 'hermite', [case], width)[0]
            self.assertEqual(value, -1.)
            self.assertEqual(bits, float_bits(-1., width))
            self.assertIsNone(error)

    def test_original_signaling_nan_input_is_preserved(self):
        case = Case([from_bits(0x7f800001, 32)], 'test', 0x7f800001)
        value, _, error = native_batch(self.runner, 'sin', [case], 32)[0]
        self.assertTrue(math.isnan(value))
        self.assertIsNone(error)
        self.assertEqual(case.input_bits, 0x7f800001)

    def test_native_timeout_preserves_successful_prefix_and_suffix(self):
        with tempfile.TemporaryDirectory() as directory:
            runner = Path(directory)/'runner'
            runner.write_text(f'#!{sys.executable}\n'+
                'import sys, struct, time\n'
                'while row := sys.stdin.buffer.read(32):\n'
                '    word = struct.unpack("<4Q", row)[0]\n'
                '    if word == 0x7f800000: time.sleep(2)\n'
                '    sys.stdout.buffer.write(struct.pack("<Q", word))\n'
                '    sys.stdout.buffer.flush()\n')
            runner.chmod(0o755)
            cases = [Case([x], 'test', float_bits(x, 32)) for x in (1., math.inf, 2.)]
            outputs = native_batch(runner, 'gamma', cases, 32, timeout=.1)
            self.assertEqual(outputs[0][0], 1.)
            self.assertIn('exceeded', outputs[1][2])
            self.assertEqual(outputs[2][0], 2.)

    def test_reference_timeout_is_recorded(self):
        from unittest.mock import patch
        evaluator = ReferenceEvaluator(5)
        try:
            arguments = ('sin', [1.], math.sin(1.), 64, 'test', 'test', 128, 256,
                         float_bits(1., 64), float_bits(math.sin(1.), 64), None)
            evaluator.evaluate(*arguments)
            with patch.object(evaluator.connection, 'poll', return_value=False):
                record = evaluator.evaluate(*arguments)
            self.assertEqual(record['classification'], 'reference_error')
            self.assertIn('exceeded', record['reference_message'])
        finally:
            evaluator.close()

    def test_reference_worker_can_restart(self):
        evaluator = ReferenceEvaluator(5)
        try:
            record = evaluator.evaluate('sin', [1.], math.sin(1.), 64, 'test', 'test', 128, 256,
                                        float_bits(1., 64), float_bits(math.sin(1.), 64), None)
            self.assertEqual(record['reference_status'], 'stable')
            evaluator.close()
            record = evaluator.evaluate('sqrt', [4.], 2., 64, 'test', 'test', 128, 256,
                                        float_bits(4., 64), float_bits(2., 64), None)
            self.assertEqual(record['classification'], 'rounded')
        finally:
            evaluator.close()
