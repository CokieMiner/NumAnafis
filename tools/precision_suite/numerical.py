"""Native evaluation, adaptive references, and IEEE error measurement."""
import hashlib
import json
import math
import multiprocessing
import os
import tempfile
from collections import deque
import shutil
import struct
import subprocess
from pathlib import Path

import mpmath as mp

from .inputs import CRATE_DIR, FUNCTIONS, INTEGER_POLE_FUNCS, float_bits, from_bits

FORMATS = {32: (24, -126, 127, 127), 64: (53, -1022, 1023, 1023)}


def build_runner(no_std=False):
    command = ['cargo', 'build', '--release', '--example', 'precision_eval', '--message-format=json']
    if no_std:
        command.append('--no-default-features')
    completed = subprocess.run(command, cwd=CRATE_DIR, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError('Cargo build failed:\n'+completed.stderr)
    for line in completed.stdout.splitlines():
        event = json.loads(line)
        if event.get('reason') == 'compiler-artifact' and event.get('target', {}).get('name') == 'precision_eval':
            if event.get('executable'):
                digest = hashlib.sha256(Path(event['executable']).read_bytes()).hexdigest()
                executable = CRATE_DIR/'target/precision-runners'/digest/'precision_eval'
                executable.parent.mkdir(parents=True, exist_ok=True)
                if not executable.exists():
                    with tempfile.NamedTemporaryFile(dir=executable.parent, delete=False) as temporary:
                        temp_path = Path(temporary.name)
                    try:
                        shutil.copy2(event['executable'], temp_path)
                        os.replace(temp_path, executable)
                    finally:
                        temp_path.unlink(missing_ok=True)
                return executable
    raise RuntimeError('Cargo did not return the precision evaluator executable')


def evaluate_case(name, args, our, width, run_id, source, ref_bits=128, max_ref_bits=1024,
                  input_bits=None, output_bits=None, native_error=None):
    if native_error:
        reference, precision, status, message = None, ref_bits, 'native_error', native_error
    else:
        reference, precision, status, message = evaluate_reference(name, args, width, ref_bits, max_ref_bits)
    word = input_bits if input_bits is not None else float_bits(args[0], width)
    record = {'run_id': run_id, 'backend': str(width), 'precision_bits': 24 if width == 32 else 53,
              'function': name, 'input': args[0], 'extras': args[1:], 'our_value': our,
              'input_bits': f'{word:0{width//4}x}',
              'output_bits': f'{output_bits:0{width//4}x}' if output_bits is not None else None, 'source': source,
              'targeted': source in ('edge', 'targeted'), 'reference_bits': precision,
              'reference_status': status, 'reference_message': message}
    if reference is None:
        record.update(reference=None, ulp_error=None, classification=status)
    else:
        record.update(error_metrics(our, reference, width, precision))
        with mp.workprec(precision):
            record['reference_hp'] = mp.nstr(reference, 60)
    return record


def native_batch(runner, name, cases, width, timeout=5):
    """Preserve successful outputs; isolate timeouts without recursive batch splitting."""
    pending = deque([cases])
    outputs = []
    while pending:
        chunk = pending.popleft()
        payload = bytearray()
        for case in chunk:
            words = [int(value) & ((1 << 64)-1) if integer else float_bits(value, width)
                     for value, (_, _, integer) in zip(case.args, FUNCTIONS[name])]
            words[0] = case.input_bits
            payload.extend(struct.pack('<4Q', *(words+[0]*(4-len(words)))))
        try:
            result = subprocess.run([str(runner), str(width), name], input=payload,
                                    capture_output=True, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            prefix = exc.stdout or b''
            if len(prefix) % 8:
                raise RuntimeError('native evaluator returned a partial output word') from exc
            completed = list(struct.iter_unpack('<Q', prefix))
            outputs.extend((from_bits(word, width), word, None) for word, in completed)
            index = len(completed)
            if index == len(chunk):
                raise RuntimeError('native evaluator timed out after returning the batch') from exc
            if len(chunk) == 1:
                outputs.append((None, None, f'native evaluation exceeded {timeout:g}s'))
            else:
                if index+1 < len(chunk):
                    pending.appendleft(chunk[index+1:])
                # Retry alone: an aggregate timeout does not prove an input is stuck.
                pending.appendleft(chunk[index:index+1])
            continue
        if result.returncode:
            raise RuntimeError('native evaluator failed: '+result.stderr.decode(errors='replace').strip())
        if len(result.stdout) != 8*len(chunk):
            raise RuntimeError('native evaluator returned an incomplete batch')
        outputs.extend((from_bits(word, width), word, None)
                       for word, in struct.iter_unpack('<Q', result.stdout))
    return outputs


class ReferenceEvaluator:
    """Isolate reference work so a difficult input cannot block an entire scan."""
    def __init__(self, timeout=5):
        self.timeout = timeout
        self.context = multiprocessing.get_context('spawn')
        self.connection = None
        self.process = None

    def evaluate(self, *arguments):
        if self.process is None:
            parent, child = self.context.Pipe()
            self.process = self.context.Process(target=reference_worker, args=(child,), daemon=True)
            self.process.start()
            child.close()
            self.connection = parent
        self.connection.send(arguments)
        if self.connection.poll(self.timeout):
            try:
                return self.connection.recv()
            except EOFError as exc:
                self.close()
                raise RuntimeError('reference worker exited unexpectedly') from exc
        self.close()
        # Construct the complete input record while bypassing further reference work.
        timed_out = list(arguments)
        timed_out[-1] = f'reference evaluation exceeded {self.timeout:g}s'
        record = evaluate_case(*timed_out)
        record.update(classification='reference_error', reference_status='reference_error')
        return record

    def close(self):
        if self.process is not None:
            self.process.terminate()
            self.process.join()
            self.connection.close()
            self.process.close()
            self.process = None
            self.connection = None


def evaluate_reference(name, args, width, bits=128, max_bits=1024):
    """Confirm target-rounded agreement at two precisions; this is not certification."""
    if any(math.isnan(a) for a in args):
        if name == 'powf' and (args[1] == 0 or args[0] == 1):
            return mp.mpf(1), bits, 'stable', None
        return mp.nan, bits, 'nonfinite', None
    if not real_domain(name, args):
        return mp.nan, bits, 'domain', None
    exceptional = exceptional_reference(name, args)
    if exceptional is not None:
        return exceptional, bits, 'nonfinite' if not math.isfinite(exceptional) else 'stable', None
    if at_pole(name, args):
        return None, bits, 'pole', 'pole sign/exception policy is not scored'
    precision = bits
    previous = None
    while precision <= max_bits:
        try:
            with mp.workprec(precision):
                value = mathematical_reference(name, args)
                if isinstance(value, mp.mpc):
                    if value.imag:
                        return None, precision, 'reference_error', 'unexpected complex reference'
                    value = value.real
                rounded = round_reference_bits(value, width)
                if previous == rounded:
                    return value, precision, 'stable', None
                previous = rounded
        except (ValueError, ZeroDivisionError, OverflowError, ArithmeticError) as exc:
            return None, precision, 'reference_error', f'{type(exc).__name__}: {exc}'
        precision *= 2
    return None, max_bits, 'reference_unstable', 'target rounding did not stabilize'


def mathematical_reference(name, args):
    x = mp.mpf(args[0])
    extras = [mp.mpf(a) for a in args[1:]]
    if name == 'lgamma':
        return mp.log(abs(mp.gamma(x)))
    if name in ('trigamma', 'tetragamma'):
        return mp.polygamma(1 if name == 'trigamma' else 2, x)
    if name in ('lambert_w0', 'lambert_wm1'):
        return mp.lambertw(x, 0 if name == 'lambert_w0' else -1)
    if name == 'polygamma':
        return mp.polygamma(int(args[1]), x)
    if name == 'zeta_deriv':
        return mp.zeta(x, derivative=int(args[1]))
    if name.startswith('bessel_'):
        return getattr(mp, 'bessel'+name[-1])(int(args[1]), x)
    if name == 'hermite':
        return mp.hermite(int(args[1]), x)
    if name == 'assoc_legendre':
        degree, order = int(args[1]), int(args[2])
        positive = (-1)**abs(order) * (1-x*x)**(mp.mpf(abs(order))/2) * mp.diff(lambda t: mp.legendre(degree, t), x, abs(order))
        return positive if order >= 0 else (-1)**order * mp.factorial(degree+order)/mp.factorial(degree-order) * positive
    if name == 'spherical_harmonic':
        degree, order = int(args[1]), int(args[2])
        # The crate defines the cosine projection, without sqrt(2) scaling.
        value = mp.spherharm(degree, abs(order), x, extras[2]).real
        return value * (-1 if order < 0 and abs(order) % 2 else 1)
    if name in ('elliptic_k', 'elliptic_e'):
        return (mp.ellipk if name == 'elliptic_k' else mp.ellipe)(x)
    if name == 'cbrt':
        return mp.sign(x)*mp.root(abs(x), 3)
    if name == 'round':
        value = mp.sign(x)*mp.floor(abs(x)+mp.mpf('.5'))
        return -0. if value == 0 and args[0] < 0 else value
    if name in ('trunc', 'ceil'):
        value = mp.ceil(x) if name == 'ceil' or x < 0 else mp.floor(x)
        return -0. if value == 0 and args[0] < 0 else value
    if name == 'fract':
        return x-mathematical_reference('trunc', args)
    if name == 'atan2' and args[0] == 0 and math.copysign(1., args[1]) < 0:
        return mp.pi*math.copysign(1., args[0])
    aliases = {'exp_m1': 'expm1', 'ln_1p': 'log1p', 'ln': 'log', 'powf': 'power'}
    return getattr(mp, aliases.get(name, name))(x, *extras)


def exceptional_reference(name, args):
    """Known IEEE limits and signed zeros; None requests high-precision evaluation."""
    x = args[0]
    zero_preserving = {'sin', 'tan', 'asin', 'atan', 'sinh', 'tanh', 'asinh', 'atanh',
                       'erf', 'lambert_w0', 'exp_m1', 'ln_1p', 'sqrt', 'cbrt', 'floor', 'ceil', 'round', 'trunc'}
    if x == 0 and name in zero_preserving:
        return x
    if name == 'atan2' and x == 0 and math.copysign(1., args[1]) > 0:
        return x
    if (name == 'ln' and x == 0) or (name == 'ln_1p' and x == -1):
        return -math.inf
    if name == 'atanh' and abs(x) == 1:
        return math.copysign(math.inf, x)
    if name == 'lambert_wm1' and x == 0:
        return -math.inf
    if name == 'powf':
        y = args[1]
        if y == 0 or x == 1:
            return 1.
        if math.isinf(y):
            return 1. if abs(x) == 1 else (math.inf if (abs(x) > 1) == (y > 0) else 0.)
        if x == 0 and y < 0:
            odd = math.isfinite(y) and y.is_integer() and abs(y) % 2 == 1
            return math.copysign(math.inf, x) if odd else math.inf
    if math.isinf(x):
        if name in ('sin', 'cos', 'tan', 'fract', 'spherical_harmonic'):
            return math.nan
        if name in ('floor', 'ceil', 'trunc', 'round', 'cbrt', 'sinh', 'asinh'):
            return x
        if name in ('exp', 'exp_m1'):
            return math.inf if x > 0 else (0. if name == 'exp' else -1.)
        if name == 'cosh':
            return math.inf
        if name in ('erf', 'tanh'):
            return math.copysign(1., x)
        if name == 'erfc':
            return 0. if x > 0 else 2.
        if x > 0:
            limits = {'lgamma': math.inf, 'ln': math.inf, 'ln_1p': math.inf, 'sqrt': math.inf,
                      'gamma': math.inf, 'digamma': math.inf, 'trigamma': 0.,
                      'tetragamma': -0., 'zeta': 1., 'lambert_w0': math.inf,
                      'bessel_y': 0., 'bessel_k': 0.}
            return limits.get(name)
    return None


def real_domain(name, args):
    x = args[0]
    if name in ('asin', 'acos', 'assoc_legendre') and abs(x) > 1:
        return False
    if name == 'acosh' and x < 1:
        return False
    if name == 'atanh' and abs(x) > 1:
        return False
    if name in ('sqrt', 'ln') and x < 0:
        return False
    if name == 'ln_1p' and x < -1:
        return False
    if name in ('elliptic_k', 'elliptic_e') and x > 1:
        return False
    if name in ('bessel_y', 'bessel_k') and x <= 0:
        return False
    if name in ('lambert_w0', 'lambert_wm1'):
        with mp.workprec(128):
            if mp.mpf(x) < -1/mp.e:
                return False
        if name == 'lambert_wm1' and x > 0:
            return False
    if name == 'powf' and x < 0 and math.isfinite(args[1]) and args[1] != math.trunc(args[1]):
        return False
    return True


def at_pole(name, args):
    x = args[0]
    if not math.isfinite(x):
        return False
    if name in INTEGER_POLE_FUNCS:
        return x <= 0 and x == math.trunc(x)
    return (name in ('zeta', 'zeta_deriv', 'elliptic_k') and x == 1.)


def error_metrics(our, reference, width, ref_bits):
    """Measure error / 2**(max(emin,floor(log2(abs(reference))))-p+1).

    Zero and subnormals use the minimum subnormal quantum. At powers of two,
    the reference binade determines spacing; the smaller spacing below that
    boundary is not substituted. Representable-step distance is independent.
    """
    expected_bits = round_reference_bits(reference, width)
    expected = from_bits(expected_bits, width)
    actual_bits = float_bits(our, width)
    if math.isnan(expected) or math.isinf(expected):
        matches = (math.isnan(our) if math.isnan(expected) else our == expected)
        return {'reference': expected, 'rounded_reference_bits': f'{expected_bits:0{width//4}x}',
                'abs_error': None, 'rel_error': None, 'ulp_error': None,
                'ulp_distance': None, 'matches_rounded_reference': matches,
                'classification': 'expected_nonfinite' if matches else 'mismatch'}
    if not math.isfinite(our):
        return {'reference': expected, 'abs_error': math.inf, 'rel_error': math.inf,
                'ulp_error': math.inf, 'ulp_distance': None,
                'matches_rounded_reference': False, 'classification': 'mismatch'}
    with mp.workprec(ref_bits):
        reference = mp.mpf(reference)
        error = abs(mp.mpf(our)-reference)
        p, emin, _, _ = FORMATS[width]
        exponent = reference._mpf_[2]+reference._mpf_[3]-1 if reference else emin
        quantum = mp.ldexp(mp.mpf(1), max(emin, exponent)-p+1)
        ulp = error/quantum
        distance = abs(ordered_bits(actual_bits, width)-ordered_bits(expected_bits, width))
        correct = actual_bits == expected_bits
        zero_sign_mismatch = distance == 0 and not correct
        classification = ('signed_zero' if zero_sign_mismatch else 'rounded' if correct else 'within_1' if ulp <= 1 else
                          'within_5' if ulp <= 5 else 'within_10' if ulp <= 10 else 'severe')
        return {'reference': expected, 'rounded_reference_bits': f'{expected_bits:0{width//4}x}',
                'abs_error': float(error), 'rel_error': float(error/abs(reference)) if reference else None,
                'ulp_error': float(ulp), 'ulp_distance': distance,
                'abs_error_hp': mp.nstr(error, 40),
                'rel_error_hp': mp.nstr(error/abs(reference), 40) if reference else None,
                'ulp_error_hp': mp.nstr(ulp, 40), 'ulp_spacing_hp': mp.nstr(quantum, 40),
                'matches_rounded_reference': correct, 'classification': classification}


def round_reference_bits(value, width):
    if isinstance(value, float):
        if value == 0:
            return float_bits(value, width)
        value = mp.mpf(value)
    p, emin, emax, bias = FORMATS[width]
    sign = int(value < 0)
    sign_word = sign << (width-1)
    infinity = ((1 << (width-p))-1) << (p-1)
    if mp.isnan(value):
        return infinity | (1 << (p-2))
    if mp.isinf(value):
        return sign_word | infinity
    if value == 0:
        return sign_word
    # mpmath stores finite values exactly as signed integer * 2**exponent.
    _, mantissa, exponent, bit_count = value._mpf_
    mantissa, exponent, bit_count = int(mantissa), int(exponent), int(bit_count)
    leading = exponent+bit_count-1
    if leading > emax:
        return sign_word | infinity
    min_exp = emin-p+1
    if leading < min_exp-1:
        return sign_word
    shift = max(bit_count-p, min_exp-exponent)
    if shift > 0:
        quotient = mantissa >> shift
        remainder = mantissa-(quotient << shift)
        half = 1 << (shift-1)
        rounded = quotient + int(remainder > half or (remainder == half and quotient & 1))
    else:
        rounded = mantissa << -shift
    if rounded == 0:
        return sign_word
    exponent += shift
    leading = exponent+rounded.bit_length()-1
    if leading > emax:
        return sign_word | infinity
    if leading < emin:
        return sign_word | (rounded << (exponent-min_exp))
    padding = p-rounded.bit_length()
    significand = rounded << padding if padding >= 0 else rounded >> -padding
    fraction = significand - (1 << (p-1))
    return sign_word | ((leading+bias) << (p-1)) | fraction


def ordered_bits(bits, width):
    sign = 1 << (width-1)
    return sign-(bits & (sign-1)) if bits & sign else sign+bits


def reference_worker(connection):
    try:
        while True:
            try:
                arguments = connection.recv()
            except EOFError:
                break
            connection.send(evaluate_case(*arguments))
    finally:
        connection.close()
