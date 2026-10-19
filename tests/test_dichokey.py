import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import dichokey  # noqa: E402

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

GOOD_KEY = """\
key: Tiny
1a leaf = needle -> 2
1b leaf = broad -> 4
2a bundles = yes -> 3
2b bundles = no -> Taxus baccata
3a cone = woody -> Pinus sylvestris
3b cone = berry -> Juniperus communis
4a lobes = deep -> Quercus robur
4b lobes = shallow -> Acer campestre
"""


def run_cli(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = dichokey.main(argv)
    return code, out.getvalue(), err.getvalue()


def messages(diags):
    return [(d.severity, d.line, d.message) for d in diags]


class ParseTests(unittest.TestCase):
    def test_parses_title_and_leads(self):
        key = dichokey.parse_key(GOOD_KEY)
        self.assertEqual(key.title, "Tiny")
        self.assertEqual(sorted(key.couplets), [1, 2, 3, 4])
        lead = key.couplets[1][0]
        self.assertEqual((lead.char, lead.state, lead.target, lead.line), ("leaf", "needle", 2, 2))
        self.assertEqual(key.couplets[2][1].target, "Taxus baccata")
        self.assertEqual(key.parse_diagnostics, [])

    def test_bad_line_is_reported_with_its_line_number(self):
        key = dichokey.parse_key("1a leaf = needle -> 2\n\n# note\n3a cone -> Pinus\n")
        self.assertEqual(messages(key.parse_diagnostics)[0][:2], ("error", 4))

    def test_duplicate_lead_is_rejected_and_points_at_first(self):
        key = dichokey.parse_key("1a leaf = needle -> 2\n1a leaf = broad -> 3\n")
        self.assertEqual(len(key.couplets[1]), 1)
        self.assertIn("first at line 1", key.parse_diagnostics[0].message)
        self.assertEqual(key.parse_diagnostics[0].line, 2)


class LintTests(unittest.TestCase):
    def test_good_key_is_clean(self):
        self.assertEqual(dichokey.lint(dichokey.parse_key(GOOD_KEY)), [])

    def test_missing_couplet_one(self):
        diags = dichokey.lint(dichokey.parse_key("2a x = y -> A\n2b x = z -> B\n"))
        self.assertTrue(any("no couplet 1" in d.message for d in diags))

    def test_incomplete_couplet_and_undefined_target(self):
        diags = dichokey.lint(dichokey.parse_key("1a leaf = needle -> 7\n"))
        msgs = [d.message for d in diags]
        self.assertIn("couplet 1 has no lead 'b'", msgs)
        self.assertIn("lead 1a points to undefined couplet 7", msgs)

    def test_mixed_characters_and_identical_states(self):
        diags = dichokey.lint(dichokey.parse_key("1a leaf = x -> A\n1b cone = X -> B\n"))
        msgs = [d.message for d in diags]
        self.assertIn("couplet 1 leads test different characters", msgs)
        self.assertIn("couplet 1 has two leads with the same state", msgs)

    def test_cycle_is_an_error(self):
        diags = dichokey.lint(dichokey.parse_key("1a a = x -> 2\n1b a = y -> A\n2a b = x -> 1\n2b b = y -> B\n"))
        errors = [d for d in diags if d.severity == "error"]
        self.assertTrue(any("makes a cycle" in d.message and d.line == 3 for d in errors))

    def test_unreachable_couplet_and_duplicate_taxon_are_warnings(self):
        text = "1a a = x -> A\n1b a = y -> A\n5a b = x -> B\n5b b = y -> C\n"
        diags = dichokey.lint(dichokey.parse_key(text))
        warnings = [d.message for d in diags if d.severity == "warning"]
        self.assertIn("couplet 5 is never reached from couplet 1", warnings)
        self.assertTrue(any("already identified at line 1" in m for m in warnings))

    def test_broken_example_has_expected_errors(self):
        key = dichokey.parse_key((EXAMPLES / "broken.key").read_text(encoding="utf-8"))
        errors = [d for d in dichokey.lint(key) if d.severity == "error"]
        self.assertTrue(any(d.line == 3 and "same state" in d.message for d in errors))
        self.assertTrue(any("undefined couplet 9" in d.message and d.line == 4 for d in errors))
        self.assertTrue(any("cannot parse" in d.message and d.line == 6 for d in errors))
        self.assertTrue(any("makes a cycle" in d.message and d.line == 7 for d in errors))


class IdentifyTests(unittest.TestCase):
    def setUp(self):
        self.key = dichokey.parse_key(GOOD_KEY)

    def test_identifies_taxon_with_full_path(self):
        result = dichokey.identify(self.key, {"leaf": "needle", "bundles": "yes", "cone": "woody"})
        self.assertEqual(result.status, "identified")
        self.assertEqual(result.taxon, "Pinus sylvestris")
        self.assertEqual([s.couplet for s in result.path], [1, 2, 3])

    def test_matching_ignores_case_and_extra_spaces(self):
        result = dichokey.identify(self.key, {"Leaf": "  BROAD ", "lobes": "Deep"})
        self.assertEqual(result.taxon, "Quercus robur")

    def test_missing_character_is_undetermined_with_candidates(self):
        result = dichokey.identify(self.key, {"leaf": "needle"})
        self.assertEqual(result.status, "undetermined")
        self.assertEqual((result.couplet, result.needs), (2, "bundles"))
        self.assertEqual(result.candidates, ("Juniperus communis", "Pinus sylvestris", "Taxus baccata"))

    def test_value_fitting_no_lead_is_no_match(self):
        result = dichokey.identify(self.key, {"leaf": "needle", "bundles": "maybe"})
        self.assertEqual(result.status, "no-match")
        self.assertEqual(result.expected, ("yes", "no"))

    def test_parse_observations_detects_conflicts_with_line_numbers(self):
        obs, diags = dichokey.parse_observations("leaf = needle\n# c\nleaf = broad\nbad line\n")
        self.assertEqual(obs, {"leaf": "needle"})
        self.assertEqual([d.line for d in diags], [3, 4])


class ToolTests(unittest.TestCase):
    def test_separate_reports_the_parting_couplet(self):
        key = dichokey.parse_key(GOOD_KEY)
        info = dichokey.separate(key, "Pinus sylvestris", "Quercus robur")
        self.assertEqual((info["couplet"], info["char"]), (1, "leaf"))
        self.assertEqual((info["state_a"], info["state_b"]), ("needle", "broad"))

    def test_separate_unknown_taxon_raises(self):
        with self.assertRaises(ValueError):
            dichokey.separate(dichokey.parse_key(GOOD_KEY), "Pinus sylvestris", "Nothing here")

    def test_tree_marks_repeated_couplets(self):
        text = "1a a = x -> 2\n1b a = y -> 2\n2a b = p -> A\n2b b = q -> B\n"
        lines = dichokey.tree_lines(dichokey.parse_key(text))
        self.assertTrue(any("(see couplet 2)" in line for line in lines))

    def test_stats_counts_depth_and_characters(self):
        data = dichokey.stats(dichokey.parse_key(GOOD_KEY))
        self.assertEqual((data["couplets"], data["leads"], data["taxa"]), (4, 8, 5))
        self.assertEqual(data["max_depth"], 3)
        self.assertEqual(data["characters"], {"leaf": 1, "bundles": 1, "cone": 1, "lobes": 1})

    def test_dot_export_has_edges_and_quoted_labels(self):
        dot = dichokey.dot_text(dichokey.parse_key(GOOD_KEY))
        self.assertTrue(dot.startswith("digraph key {"))
        self.assertIn("c1 -> c2", dot)
        self.assertIn('label="a) needle"', dot)

    def test_batch_identifies_rows_and_reports_bad_rows(self):
        key = dichokey.parse_key(GOOD_KEY)
        text = "specimen,leaf,bundles,cone\ns1,needle,yes,woody\ns2,needle,no\ns3,needle,yes,\n,needle,no,\n"
        results, diags = dichokey.run_batch(key, text)
        self.assertEqual([name for name, _ in results], ["s1", "s3"])
        self.assertEqual(results[0][1].taxon, "Pinus sylvestris")
        self.assertEqual(results[1][1].status, "undetermined")
        self.assertEqual(results[1][1].needs, "cone")
        self.assertEqual([d.line for d in diags], [3, 5])

    def test_batch_without_specimen_column_fails(self):
        _, diags = dichokey.run_batch(dichokey.parse_key(GOOD_KEY), "name,leaf\nx,needle\n")
        self.assertIn("no 'specimen' column", diags[0].message)


class CliTests(unittest.TestCase):
    def test_lint_example_gate_exit_codes(self):
        code, out, _ = run_cli(["lint", str(EXAMPLES / "trees.key")])
        self.assertEqual(code, 0)
        self.assertIn("OK (4 couplets, 5 taxa)", out)
        code, out, _ = run_cli(["lint", str(EXAMPLES / "broken.key")])
        self.assertEqual(code, 1)
        self.assertIn("same state", out)
        self.assertIn("cannot parse", out)

    def test_identify_json_output(self):
        key_path = str(EXAMPLES / "trees.key")
        code, out, _ = run_cli(["identify", key_path, "--obs", "leaf=broad", "--obs", "lobes=shallow", "--json"])
        self.assertEqual(code, 0)
        self.assertIn('"taxon": "Acer campestre"', out)

    def test_batch_example_tallies(self):
        code, out, _ = run_cli(["batch", str(EXAMPLES / "trees.key"), str(EXAMPLES / "specimens.csv")])
        self.assertEqual(code, 0)
        self.assertIn("5 specimens: 3 identified, 1 undetermined, 1 no match, 0 row error(s)", out)
        self.assertIn("spec-05\tundetermined (needs cone)", out)

    def test_batch_strict_fails_on_unresolved_rows(self):
        code, _, _ = run_cli(["batch", str(EXAMPLES / "trees.key"), str(EXAMPLES / "specimens.csv"), "--strict"])
        self.assertEqual(code, 1)

    def test_bad_key_blocks_identify(self):
        code, _, err = run_cli(["identify", str(EXAMPLES / "broken.key"), "--obs", "leaf=needle"])
        self.assertEqual(code, 1)
        self.assertIn("has errors", err)

    def test_missing_file_is_reported_not_raised(self):
        code, _, err = run_cli(["stats", "/nonexistent/none.key"])
        self.assertEqual(code, 1)
        self.assertIn("dichokey:", err)

    def test_dot_writes_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "key.dot"
            code, _, _ = run_cli(["dot", str(EXAMPLES / "trees.key"), "-o", str(target)])
            self.assertEqual(code, 0)
            self.assertIn("digraph key", target.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
