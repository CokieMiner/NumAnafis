"""Audit positive and negative controls, including lexer evasions."""

import unittest
from unittest.mock import patch

from tools.audit.import_rules import FileAnalyzer
from tools.audit.common import ROOT

from tools.audit.import_checks import import_structure_findings
from tools.audit.items import cfg_requires_test, top_level_items
from tools.audit.rust_source import clean_rust_code
from tools.audit.structure_checks import structural_findings


class StructureTests(unittest.TestCase):
    def findings(self, text, registry=True):
        return {f.kind for f in structural_findings(text, "src/demo/mod.rs" if registry else "src/demo/value.rs", registry=registry, test_file=False)}

    def test_valid_registry(self):
        text = '//! Registry.\nuse super::Thing;\nmod value;\npub use value::Value;\n#[cfg(test)]\nmod tests;'
        self.assertEqual(self.findings(text), set())

    def test_registry_code_even_on_same_line(self):
        self.assertIn("implementation_in_module_registry", self.findings("mod value; const HIDDEN: u8 = 1;"))

    def test_multiline_function_declaration_is_code(self):
        self.assertIn("implementation_in_module_registry", self.findings("pub\nasync\nfn entry() {}"))

    def test_tests_must_be_last_and_gated(self):
        kinds = self.findings("mod tests; mod value;")
        self.assertIn("test_module_without_test_gate", kinds)
        self.assertIn("module_registry_order", kinds)

    def test_cfg_any_is_not_a_test_gate(self):
        self.assertFalse(cfg_requires_test('cfg(any(test, feature = "x"))'))
        self.assertTrue(cfg_requires_test('cfg(all(test, feature = "x"))'))
        self.assertFalse(cfg_requires_test("cfg(not(test))"))
        self.assertTrue(cfg_requires_test("cfg(not(not(test)))"))

    def test_inline_tests_are_found_outside_registry(self):
        kinds = self.findings("#[cfg(test)] mod checks { #[test] fn check() {} }", registry=False)
        self.assertIn("inline_test_module", kinds)
        self.assertIn("test_in_production_file", kinds)

    def test_expect_needs_reason_in_cfg_attr(self):
        for text in ('#[expect(dead_code)] fn f() {}', '#[cfg_attr(test, expect(dead_code))] fn f() {}',
                     '#[expect(dead_code, reason = "")] fn f() {}', '#[expect(dead_code /* reason = "fake" */)] fn f() {}'):
            with self.subTest(text=text):
                self.assertIn("expect_without_reason", self.findings(text, registry=False))

    def test_expect_real_reason(self):
        self.assertEqual(self.findings('#[expect(dead_code, reason = "target-specific entry")] fn f() {}', registry=False), set())

    def test_comment_and_string_decoys(self):
        text = '// mod fake { todo!(); }\n/* nested /* fn f() {} */ comment */\nmod real;\n#[doc = r#"fn fake() {}"#]\npub use real::Thing;'
        self.assertEqual(self.findings(text), set())

    def test_escaped_newline_preserves_positions(self):
        source = 'const TEXT: &str = "first\\\nsecond";\nfn f() {}'
        cleaned = clean_rust_code(source)
        self.assertEqual(len(source), len(cleaned))
        self.assertEqual([i for i,c in enumerate(source) if c == '\n'], [i for i,c in enumerate(cleaned) if c == '\n'])

    def test_unknown_macro_is_not_a_registry_escape(self):
        self.assertIn("implementation_in_module_registry", self.findings("generate_code!();"))

    def test_invalid_source_is_not_silently_skipped(self):
        self.assertIn("unparsed_rust_structure", self.findings("mod broken {"))


class ImportTests(unittest.TestCase):
    def kinds(self, text, parent=None):
        return {f.kind for f in import_structure_findings(text, "src/demo/value.rs", parent)}

    def test_multiline_use_tree_preserves_scope_and_group(self):
        text = 'use core::cmp::Ordering;\n\nuse crate::{\n    Thing,\n};\n\nuse super::Other;'
        path = ROOT/'src/demo/value.rs'
        with patch('pathlib.Path.read_text', return_value=text):
            analyzer = FileAnalyzer(path, set())
            analyzer.analyze()
        self.assertEqual([item[3] for item in analyzer.private_imports], [0, 3, 4])
        self.assertNotIn('import_group_spacing', {item[0] for item in analyzer.finding_keys})

    def test_grouped_imports_do_not_bypass_order(self):
        self.assertIn("grouped_import_path_order", self.kinds("use core::ops::{Add, Sub};\nuse core::fmt::{Debug, Display};"))

    def test_rustfmt_single_before_grouped_import(self):
        self.assertNotIn("grouped_import_path_order", self.kinds("use super::Z;\nuse super::{A, B};"))

    def test_duplicate_imports_and_tree_entries(self):
        self.assertIn("duplicate_import_in_tree", self.kinds("use core::{fmt::Debug, fmt::Debug};"))
        self.assertIn("duplicate_import", self.kinds("use core::fmt::Debug;\nuse core::fmt::Debug;"))

    def test_distinct_cfg_alternatives_are_not_duplicates(self):
        self.assertNotIn("duplicate_import", self.kinds('#[cfg(feature="a")] use core::fmt::Debug;\n#[cfg(not(feature="a"))] use core::fmt::Debug;'))
        self.assertNotIn("duplicate_import", self.kinds('#[cfg(target_arch="arm")] use core::fmt::Debug;\n#[cfg(target_arch="x86")] use core::fmt::Debug;'))

    def test_private_sibling_prelude_is_found(self):
        self.assertIn("private_sibling_prelude", self.kinds("use super::Thing;", "mod child;\nuse child::Thing;"))
        self.assertNotIn("private_sibling_prelude", self.kinds("use super::Thing;", "mod child;\npub use child::Thing;"))

    def test_each_import_needs_own_line(self):
        self.assertIn("multiple_imports_on_line", self.kinds("use core::fmt::Debug; use core::fmt::Display;"))

    def test_use_tree_does_not_end_at_its_brace(self):
        items = top_level_items("use core::{fmt::{Debug, Display}, ops::Add};\nfn value() {}")
        self.assertEqual([item.kind for item in items], ["use", "fn"])


if __name__ == "__main__":
    unittest.main()
