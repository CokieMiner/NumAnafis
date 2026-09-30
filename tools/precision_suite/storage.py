"""Transactional SQLite checkpoints, retained cases, and complete-run summaries."""
import json
import math
import sqlite3
from collections import Counter
from decimal import Decimal
from pathlib import Path
SCHEMA_VERSION = 1


class Store:
    def __init__(self, path, readonly=False):
        self.path = Path(path)
        if not readonly:
            self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.path.resolve().as_uri()+'?mode=ro', uri=True) if readonly else sqlite3.connect(self.path, timeout=30)
        version = self.db.execute('PRAGMA user_version').fetchone()[0]
        if version not in (0, SCHEMA_VERSION):
            self.db.close()
            raise ValueError(f'unsupported database schema {version}')
        if readonly:
            if version != SCHEMA_VERSION:
                self.db.close()
                raise ValueError('file is not a precision-suite database')
            return
        self.db.execute('PRAGMA foreign_keys=ON')
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.executescript('''
            CREATE TABLE IF NOT EXISTS runs (
                id TEXT PRIMARY KEY, configuration TEXT NOT NULL, created TEXT NOT NULL,
                completed INTEGER NOT NULL DEFAULT 0);
            CREATE TABLE IF NOT EXISTS progress (
                run_id TEXT NOT NULL REFERENCES runs(id), width INTEGER NOT NULL,
                function TEXT NOT NULL, cursor INTEGER NOT NULL, summary TEXT NOT NULL,
                PRIMARY KEY(run_id,width,function));
            CREATE TABLE IF NOT EXISTS cases (
                run_id TEXT NOT NULL, width INTEGER NOT NULL, function TEXT NOT NULL,
                slot INTEGER NOT NULL, position INTEGER NOT NULL, record TEXT NOT NULL,
                PRIMARY KEY(run_id,width,function,slot));
            CREATE INDEX IF NOT EXISTS case_function ON cases(function,width,run_id);
        ''')
        self.db.execute(f'PRAGMA user_version={SCHEMA_VERSION}')
        self.db.commit()

    def start(self, run_id, configuration, created, resume=None):
        if resume:
            row = self.db.execute('SELECT id,configuration,completed FROM runs WHERE id=?', (resume,)).fetchone()
            if row is None:
                raise ValueError(f'unknown run {resume}')
            if loads(row[1]) != configuration:
                raise ValueError('resume configuration or evaluator/reference code differs from the saved run')
            if row[2]:
                raise ValueError('run is already complete')
            return row[0]
        with self.db:
            self.db.execute('INSERT INTO runs(id,configuration,created) VALUES(?,?,?)',
                            (run_id, dumps(configuration), created))
        return run_id

    def progress(self, run_id, width, name):
        row = self.db.execute('SELECT cursor,summary FROM progress WHERE run_id=? AND width=? AND function=?',
                              (run_id, width, name)).fetchone()
        return (row[0], loads(row[1])) if row else (0, {'count': 0, 'counts': {}, 'max_ulp': 0., 'max_ulp_hp': '0', 'worst': [], 'failures': []})

    def commit_batch(self, run_id, width, name, cursor, records, exhaustive):
        # Serialize the read/update so two resumptions cannot count the same batch.
        with self.db:
            self.db.execute('BEGIN IMMEDIATE')
            previous, summary = self.progress(run_id, width, name)
            if cursor != previous+len(records):
                raise ValueError('checkpoint cursor does not follow the previous batch')
            counts = Counter(summary['counts'])
            scored = [record for record in records if record.get('ulp_error') is not None]
            counts.update(record['classification'] for record in records)
            summary['counts'] = dict(counts)
            summary['count'] += len(records)
            summary['max_ulp'] = max([summary['max_ulp'], *[record['ulp_error'] for record in scored]])
            prior = {'ulp_error': summary['max_ulp'], 'ulp_error_hp': summary['max_ulp_hp']}
            maximum = max([prior, *scored], key=score)
            summary['max_ulp_hp'] = maximum.get('ulp_error_hp', str(maximum['ulp_error']))
            summary['worst'] = sorted(summary['worst']+scored, key=score, reverse=True)[:20]
            summary['failures'] = (summary['failures']+[record for record in records if
                record['classification'] in ('native_error', 'mismatch', 'signed_zero', 'reference_error', 'reference_unstable')])[:20]
            retained = {}
            for offset, record in enumerate(records):
                position = previous+offset
                # At most one candidate per bin per batch; complete counts are independent.
                slot = int(record['input_bits'], 16) >> (width-12) if exhaustive else position
                if slot not in retained or score(record) > score(retained[slot][1]):
                    retained[slot] = (position, record)
            for slot, (position, record) in retained.items():
                if exhaustive:
                    existing = self.db.execute('SELECT record FROM cases WHERE run_id=? AND width=? AND function=? AND slot=?',
                                                (run_id, width, name, slot)).fetchone()
                    if existing and score(record) <= score(loads(existing[0])):
                        continue
                self.db.execute('INSERT OR REPLACE INTO cases VALUES(?,?,?,?,?,?)',
                                (run_id, width, name, slot, position, dumps(record)))
            self.db.execute('INSERT OR REPLACE INTO progress VALUES(?,?,?,?,?)',
                            (run_id, width, name, cursor, dumps(summary)))
        return summary

    def complete(self, run_id):
        with self.db:
            self.db.execute('UPDATE runs SET completed=1 WHERE id=?', (run_id,))

    def latest(self):
        row = self.db.execute('SELECT id FROM runs ORDER BY rowid DESC LIMIT 1').fetchone()
        if row is None:
            raise ValueError('database has no runs')
        return row[0]

    def summaries(self, run_id):
        return [(width, name, cursor, loads(summary)) for width, name, cursor, summary in
                self.db.execute('SELECT width,function,cursor,summary FROM progress WHERE run_id=? ORDER BY width,function', (run_id,))]

    def close(self):
        self.db.close()

    def entries(self, run_id=None, width=None, function=None):
        query = 'SELECT record FROM cases WHERE run_id=?'
        arguments = [run_id or self.latest()]
        if width:
            query += ' AND width=?'
            arguments.append(int(width))
        if function:
            query += ' AND function=?'
            arguments.append(function)
        for record, in self.db.execute(query+' ORDER BY width,function,position', arguments):
            yield loads(record)


def score(record):
    """Order positive errors even when the decimal exponent exceeds machine limits."""
    value = record.get('ulp_error')
    if value is None:
        return (-1, 0, Decimal(0))
    text = record.get('ulp_error_hp', str(value)).lower()
    if text in ('inf', 'infinity'):
        return (2, 0, Decimal(0))
    mantissa, _, exponent = text.partition('e')
    decimal = Decimal(mantissa)
    if not decimal:
        return (0, 0, Decimal(0))
    sign, digits, shift = decimal.as_tuple()
    power = int(exponent or '0')+shift+len(digits)-1
    normalized = Decimal((sign, digits, 1-len(digits)))
    return (1, power, normalized)


def dumps(value):
    return json.dumps(encode(value), allow_nan=False, separators=(',', ':'))


def loads(value):
    return decode(json.loads(value))


def encode(value):
    if isinstance(value, float) and not math.isfinite(value):
        return {'$float': 'nan' if math.isnan(value) else 'inf' if value > 0 else '-inf'}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def decode(value):
    if isinstance(value, dict) and set(value) == {'$float'}:
        return float(value['$float'])
    if isinstance(value, dict):
        return {key: decode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [decode(item) for item in value]
    return value
