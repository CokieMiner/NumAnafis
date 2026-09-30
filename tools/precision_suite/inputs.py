"""Function domains and deterministic represented-input generation."""
import math
import random
import struct
from dataclasses import dataclass
from pathlib import Path

CRATE_DIR = Path(__file__).resolve().parents[2]

DB_PATH = CRATE_DIR / "tools/results/precision.sqlite3"

SAMPLES = 1000

TARGETED_SAMPLES = 100

FUNCTIONS = {name: [(-10., 10., False)] for name in
             "sin cos tan atan sinh cosh tanh asinh cbrt floor ceil round trunc fract erf erfc".split()}

FUNCTIONS.update({
    "asin": [(-1., 1., False)], "acos": [(-1., 1., False)],
    "acosh": [(1., 100., False)], "atanh": [(-.999, .999, False)],
    "exp": [(-100., 100., False)], "exp_m1": [(-100., 100., False)],
    "ln": [(1e-30, 100., False)], "ln_1p": [(-.999, 100., False)],
    "sqrt": [(0., 100., False)],
    "atan2": [(-10., 10., False), (-10., 10., False)],
    "powf": [(.001, 100., False), (-5., 5., False)],
    "gamma": [(-30.5, 200., False)], "lgamma": [(-30.5, 200., False)],
    "digamma": [(-30.5, 30.5, False)], "trigamma": [(-30.5, 30.5, False)],
    "tetragamma": [(-30.5, 30.5, False)], "zeta": [(-20.5, 50., False)],
    "lambert_w0": [(-1/math.e, 1000., False)],
    "lambert_wm1": [(-1/math.e, -1e-8, False)],
    "elliptic_k": [(-10., .999, False)], "elliptic_e": [(-10., .999, False)],
    "bessel_j": [(-100., 100., False), (-10, 10, True)],
    "bessel_y": [(.001, 100., False), (-10, 10, True)],
    "bessel_i": [(-100., 100., False), (-10, 10, True)],
    "bessel_k": [(.001, 100., False), (-10, 10, True)],
    "polygamma": [(-15.5, 15.5, False), (0, 5, True)],
    "beta": [(-15.5, 15.5, False), (-15.5, 15.5, False)],
    "zeta_deriv": [(-10.5, 5., False), (0, 5, True)],
    "hermite": [(-10., 10., False), (0, 10, True)],
    "assoc_legendre": [(-.999, .999, False), (0, 10, True), (-10, 10, True)],
    "spherical_harmonic": [(0., math.pi, False), (0, 10, True), (-10, 10, True), (0., 2*math.pi, False)],
})

INTEGER_POLE_FUNCS = {"gamma", "lgamma", "digamma", "trigamma", "tetragamma", "polygamma"}


@dataclass(frozen=True)
class Case:
    args: list
    source: str
    input_bits: int


def sample_cases(name, width, samples, targeted, seed, mode='random'):
    rng = random.Random(f'{seed}:{width}:{name}')
    spec = FUNCTIONS[name]
    def extras():
        values = [rng.randint(lo, hi) if integer else rng.uniform(lo, hi)
                  for lo, hi, integer in spec[1:]]
        if name in ('assoc_legendre', 'spherical_harmonic'):
            values[1] = rng.randint(-values[0], values[0])
        return values
    for i in range(samples):
        lo, hi, _ = spec[0]
        if mode == 'bits':
            raw_bits = rng.getrandbits(width)
            x = from_bits(raw_bits, width)
        elif mode == 'grid':
            x = lo + (hi-lo) * i/max(1, samples-1)
        else:
            x = rng.uniform(lo, hi)
        arguments = finalize_args(name, [x, *extras()], width)
        yield Case(arguments, mode, raw_bits if mode == 'bits' else float_bits(arguments[0], width))
    anchors = [0., -0., 1., -1., spec[0][0], spec[0][1]]
    if name == 'lgamma':
        anchors += [2.]
    if name in INTEGER_POLE_FUNCS:
        anchors += [-float(k) for k in range(31)]
    if name in ('lambert_w0', 'lambert_wm1'):
        anchors += [-1/math.e]
    anchors += [math.inf, -math.inf, math.nan]
    for anchor in anchors:
        for x in neighbors(anchor, width):
            arguments = finalize_args(name, [x, *extras()], width)
            yield Case(arguments, 'edge', float_bits(arguments[0], width))
    for _ in range(targeted):
        if name in INTEGER_POLE_FUNCS:
            x = -rng.randint(0, 30) + rng.choice([-1, 1])*10**rng.uniform(-15, -2)
        elif name in ('elliptic_k', 'elliptic_e', 'zeta', 'zeta_deriv'):
            x = 1. + rng.choice([-1, 1])*10**rng.uniform(-15, -2)
        elif name in ('lambert_w0', 'lambert_wm1'):
            x = -1/math.e + 10**rng.uniform(-15, -2)
        else:
            x = rng.choice([-1, 1])*10**rng.uniform(-35, -2)
        arguments = finalize_args(name, [x, *extras()], width)
        yield Case(arguments, 'targeted', float_bits(arguments[0], width))


def finalize_args(name, args, width):
    spec = FUNCTIONS[name]
    if len(args) != len(spec):
        raise ValueError(f'{name} needs {len(spec)} arguments')
    return [int(value) if is_int else round_input(value, width)
            for value, (_, _, is_int) in zip(args, spec)]


def neighbors(value, width):
    rounded = round_input(value, width)
    bits = float_bits(rounded, width)
    if math.isnan(rounded):
        return [rounded]
    if rounded == 0:
        return [-from_bits(1, width), -0., 0., from_bits(1, width)]
    return [from_bits(candidate, width) for candidate in (bits-1, bits, bits+1)
            if 0 <= candidate < 1 << width]


def round_input(value, width):
    return from_bits(float_bits(value, width), width)


def float_bits(value, width):
    try:
        return int.from_bytes(struct.pack('<f' if width == 32 else '<d', value), 'little')
    except OverflowError:
        return float_bits(math.copysign(math.inf, value), width)


def from_bits(bits, width):
    return struct.unpack('<f' if width == 32 else '<d', bits.to_bytes(width//8, 'little'))[0]
