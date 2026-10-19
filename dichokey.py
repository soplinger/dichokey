"""Dichokey: offline dichotomous identification-key toolkit.

A key is a plain-text file of couplets. Each couplet has an "a" lead and a
"b" lead that test the same character. A lead names either the next couplet
or the taxon it identifies:

    key: Common trees
    1a leaf = needle -> 2
    1b leaf = broad -> 4

Run `python3 dichokey.py --help` for the commands.
"""

import argparse
import csv
import io
import json
import re
import sys
from collections import Counter, namedtuple
from dataclasses import dataclass, field
from pathlib import Path

LEAD_RE = re.compile(r"^(\d+)([ab])\s+([^\s=]+)\s*=\s*(.+?)\s*->\s*(\S.*)$")
TITLE_RE = re.compile(r"^key:\s*(.+)$", re.IGNORECASE)
OBS_RE = re.compile(r"^([^\s=]+)\s*=\s*(.+)$")

Step = namedtuple("Step", "couplet char state target")


def norm(text):
    """Case- and whitespace-insensitive form used for every comparison."""
    return " ".join(text.split()).lower()


@dataclass
class Diagnostic:
    severity: str  # "error" or "warning"
    message: str
    line: int | None = None

    def format(self, source):
        where = f"{source}:{self.line}" if self.line else source
        return f"{where}: {self.severity}: {self.message}"


@dataclass
class Lead:
    number: int
    letter: str
    char: str
    state: str
    target: object  # int couplet number, or a taxon name (str)
    line: int


@dataclass
class Key:
    title: str = ""
    couplets: dict = field(default_factory=dict)  # number -> [Lead, ...]
    parse_diagnostics: list = field(default_factory=list)


@dataclass
class Result:
    status: str  # "identified", "undetermined" or "no-match"
    path: tuple = ()
    taxon: str | None = None
    couplet: int | None = None
    needs: str | None = None
    observed: str | None = None
    expected: tuple = ()
    candidates: tuple = ()

    def to_dict(self):
        return {
            "status": self.status,
            "taxon": self.taxon,
            "couplet": self.couplet,
            "needs": self.needs,
            "observed": self.observed,
            "expected": list(self.expected),
            "candidates": list(self.candidates),
            "path": [step._asdict() for step in self.path],
        }


def parse_key(text):
    """Parse key text. Malformed lines are recorded, not raised."""
    key = Key()
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        title = TITLE_RE.match(line)
        if title:
            key.title = title.group(1).strip()
            continue
        match = LEAD_RE.match(line)
        if not match:
            key.parse_diagnostics.append(Diagnostic(
                "error", f"cannot parse {line!r}; expected e.g. '1a leaf = needle -> 2'", lineno))
            continue
        number, letter, char, state, target = match.groups()
        number = int(number)
        target = int(target) if target.isdigit() else " ".join(target.split())
        leads = key.couplets.setdefault(number, [])
        previous = [lead for lead in leads if lead.letter == letter]
        if previous:
            key.parse_diagnostics.append(Diagnostic(
                "error", f"duplicate lead {number}{letter} (first at line {previous[0].line})", lineno))
            continue
        leads.append(Lead(number, letter, norm(char), " ".join(state.split()), target, lineno))
    return key


def lint(key):
    """Return every problem with the key as Diagnostic objects."""
    diags = list(key.parse_diagnostics)
    if not key.couplets:
        return diags + [Diagnostic("error", "key has no couplets")]
    if 1 not in key.couplets:
        diags.append(Diagnostic("error", "no couplet 1; the key must start at couplet 1"))
    taxa_seen = {}
    for number in sorted(key.couplets):
        leads = key.couplets[number]
        letters = {lead.letter for lead in leads}
        for letter in sorted({"a", "b"} - letters):
            diags.append(Diagnostic("error", f"couplet {number} has no lead '{letter}'", leads[0].line))
        if len({lead.char for lead in leads}) > 1:
            diags.append(Diagnostic("error", f"couplet {number} leads test different characters", leads[1].line))
        if len({norm(lead.state) for lead in leads}) < len(leads):
            diags.append(Diagnostic("error", f"couplet {number} has two leads with the same state", leads[1].line))
        for lead in leads:
            if isinstance(lead.target, int):
                if lead.target not in key.couplets:
                    diags.append(Diagnostic(
                        "error", f"lead {number}{lead.letter} points to undefined couplet {lead.target}", lead.line))
                continue
            name = norm(lead.target)
            if name in taxa_seen:
                diags.append(Diagnostic(
                    "warning", f"taxon {lead.target!r} already identified at line {taxa_seen[name]}", lead.line))
            else:
                taxa_seen[name] = lead.line
    if not taxa_seen:
        diags.append(Diagnostic("error", "key has no taxa"))
    diags.extend(_cycle_diagnostics(key))
    diags.extend(_unreachable_diagnostics(key))
    return diags


def _cycle_diagnostics(key):
    diags = []
    state = {}

    def visit(number):
        state[number] = "active"
        for lead in key.couplets[number]:
            target = lead.target
            if not isinstance(target, int) or target not in key.couplets:
                continue
            if state.get(target) == "active":
                diags.append(Diagnostic(
                    "error", f"lead {number}{lead.letter} returns to couplet {target}, which makes a cycle",
                    lead.line))
            elif target not in state:
                visit(target)
        state[number] = "done"

    for number in sorted(key.couplets):
        if number not in state:
            visit(number)
    return diags


def reachable_couplets(key):
    if 1 not in key.couplets:
        return set()
    seen = {1}
    stack = [1]
    while stack:
        number = stack.pop()
        for lead in key.couplets[number]:
            target = lead.target
            if isinstance(target, int) and target in key.couplets and target not in seen:
                seen.add(target)
                stack.append(target)
    return seen


def _unreachable_diagnostics(key):
    if 1 not in key.couplets:
        return []
    reachable = reachable_couplets(key)
    return [
        Diagnostic("warning", f"couplet {number} is never reached from couplet 1", key.couplets[number][0].line)
        for number in sorted(key.couplets) if number not in reachable
    ]


def leaves_from(key, number, prefix=(), active=frozenset()):
    """Yield (taxon, steps) for every taxon reachable from couplet `number`."""
    if number in active or number not in key.couplets:
        return
    for lead in key.couplets[number]:
        path = prefix + (Step(number, lead.char, lead.state, lead.target),)
        if isinstance(lead.target, int):
            yield from leaves_from(key, lead.target, path, active | {number})
        else:
            yield lead.target, path


def count_taxa(key):
    return len({norm(taxon) for taxon, _ in leaves_from(key, 1)})


def parse_observations(text, source="observations"):
    """Parse 'character = state' lines into a dict with normalised keys and values."""
    observations = {}
    diags = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = OBS_RE.match(line)
        if not match:
            diags.append(Diagnostic("error", f"cannot parse {line!r}; expected 'character = state'", lineno))
            continue
        char, state = norm(match.group(1)), norm(match.group(2))
        if char in observations and observations[char] != state:
            diags.append(Diagnostic("error", f"conflicting values for {char!r} in {source}", lineno))
            continue
        observations[char] = state
    return observations, diags


def identify(key, observations):
    """Walk the key with observed character states and report where it stops."""
    observed = {norm(char): norm(value) for char, value in observations.items()}
    path = []
    number = 1
    visited = set()
    while True:
        if number not in key.couplets:
            raise ValueError(f"couplet {number} is not defined")
        if number in visited:
            raise ValueError(f"couplet {number} loops back on itself")
        visited.add(number)
        leads = key.couplets[number]
        char = leads[0].char
        if char not in observed:
            candidates = sorted({taxon for taxon, _ in leaves_from(key, number)})
            return Result("undetermined", path=tuple(path), couplet=number, needs=char,
                          candidates=tuple(candidates))
        value = observed[char]
        matches = [lead for lead in leads if norm(lead.state) == value]
        if not matches:
            return Result("no-match", path=tuple(path), couplet=number, needs=char, observed=value,
                          expected=tuple(lead.state for lead in leads))
        lead = matches[0]
        path.append(Step(number, char, lead.state, lead.target))
        if isinstance(lead.target, int):
            number = lead.target
        else:
            return Result("identified", path=tuple(path), taxon=lead.target)


def find_taxon(key, name):
    wanted = norm(name)
    for taxon, path in leaves_from(key, 1):
        if norm(taxon) == wanted:
            return taxon, path
    raise ValueError(f"taxon {name!r} is not in the key")


def separate(key, name_a, name_b):
    """Return the couplet where two taxa's routes first part, and the states used there."""
    taxon_a, path_a = find_taxon(key, name_a)
    taxon_b, path_b = find_taxon(key, name_b)
    if norm(taxon_a) == norm(taxon_b):
        raise ValueError("both names are the same taxon")
    for step_a, step_b in zip(path_a, path_b):
        if step_a != step_b:
            if step_a.couplet != step_b.couplet:
                raise ValueError("routes diverge at different couplets; the key is inconsistent")
            return {
                "couplet": step_a.couplet,
                "char": step_a.char,
                "taxon_a": taxon_a,
                "state_a": step_a.state,
                "taxon_b": taxon_b,
                "state_b": step_b.state,
            }
    raise ValueError("one route is a prefix of the other; the key is inconsistent")


def tree_lines(key):
    if 1 not in key.couplets:
        return ["(no couplet 1)"]
    lines = []
    seen = set()

    def walk(number, depth):
        pad = "  " * depth
        if number in seen:
            lines.append(f"{pad}(see couplet {number})")
            return
        seen.add(number)
        leads = key.couplets.get(number)
        if leads is None:
            lines.append(f"{pad}(undefined couplet {number})")
            return
        lines.append(f"{pad}{number}. {leads[0].char}")
        for lead in leads:
            lines.append(f"{pad}  {lead.letter}) {lead.state} -> {lead.target}")
            if isinstance(lead.target, int):
                walk(lead.target, depth + 2)

    walk(1, 0)
    return lines


def stats(key):
    paths = list(leaves_from(key, 1)) if 1 in key.couplets else []
    depths = [len(path) for _, path in paths]
    return {
        "couplets": len(key.couplets),
        "leads": sum(len(leads) for leads in key.couplets.values()),
        "taxa": count_taxa(key),
        "max_depth": max(depths, default=0),
        "mean_depth": round(sum(depths) / len(depths), 2) if depths else 0.0,
        "characters": dict(Counter(leads[0].char for leads in key.couplets.values()).most_common()),
    }


def _dot(text):
    return text.replace("\\", "\\\\").replace('"', '\\"')


def dot_text(key):
    out = ["digraph key {", "  rankdir=TB;", '  node [fontname="Helvetica"];']
    for number in sorted(key.couplets):
        out.append(f'  c{number} [shape=box, label="{number}. {_dot(key.couplets[number][0].char)}"];')
    taxon_ids = {}
    for number in sorted(key.couplets):
        for lead in key.couplets[number]:
            if isinstance(lead.target, int):
                target = f"c{lead.target}"
            else:
                name = norm(lead.target)
                if name not in taxon_ids:
                    taxon_ids[name] = f"t{len(taxon_ids) + 1}"
                    out.append(f'  {taxon_ids[name]} [shape=ellipse, style=italic, label="{_dot(lead.target)}"];')
                target = taxon_ids[name]
            out.append(f'  c{number} -> {target} [label="{lead.letter}) {_dot(lead.state)}"];')
    out.append("}")
    return "\n".join(out) + "\n"


def run_batch(key, text):
    """Identify every specimen row. Returns ([(specimen, Result)], [Diagnostic])."""
    reader = csv.reader(io.StringIO(text))
    try:
        header = next(reader)
    except StopIteration:
        return [], [Diagnostic("error", "specimen file is empty")]
    names = [norm(column) for column in header]
    if "specimen" not in names:
        return [], [Diagnostic("error", "header has no 'specimen' column", 1)]
    specimen_col = names.index("specimen")
    results, diags = [], []
    for row in reader:
        line = reader.line_num
        if not any(cell.strip() for cell in row):
            continue
        if len(row) != len(header):
            diags.append(Diagnostic("error", f"row has {len(row)} fields; header has {len(header)}", line))
            continue
        specimen = row[specimen_col].strip()
        if not specimen:
            diags.append(Diagnostic("error", "row has no specimen name", line))
            continue
        observations = {names[i]: row[i] for i in range(len(row))
                        if i != specimen_col and row[i].strip()}
        results.append((specimen, identify(key, observations)))
    return results, diags


def describe(result):
    lines = []
    if result.status == "identified":
        lines.append(f"Identified: {result.taxon}")
    elif result.status == "undetermined":
        lines.append(f"Undetermined at couplet {result.couplet}: need a value for '{result.needs}'.")
        if result.candidates:
            lines.append("Still possible: " + ", ".join(result.candidates))
    else:
        lines.append(f"No match at couplet {result.couplet}: '{result.needs}' = '{result.observed}' "
                     f"fits none of {', '.join(result.expected)}.")
    if result.path:
        lines.append("Path:")
        for step in result.path:
            lines.append(f"  couplet {step.couplet}: {step.char} = {step.state} -> {step.target}")
    return "\n".join(lines)


def read_text(path):
    return Path(path).read_text(encoding="utf-8")


def load_checked(path):
    """Parse and lint a key; print blocking errors and return None if it is unusable."""
    key = parse_key(read_text(path))
    errors = [d for d in lint(key) if d.severity == "error"]
    if errors:
        for diag in errors:
            print(diag.format(path), file=sys.stderr)
        print(f"{path}: has errors; run 'dichokey.py lint {path}' for details", file=sys.stderr)
        return None
    return key


def cmd_lint(args):
    key = parse_key(read_text(args.file))
    diags = lint(key)
    for diag in diags:
        print(diag.format(args.file))
    errors = sum(diag.severity == "error" for diag in diags)
    warnings = len(diags) - errors
    if diags:
        print(f"{args.file}: {errors} error(s), {warnings} warning(s)")
    else:
        print(f"{args.file}: OK ({len(key.couplets)} couplets, {count_taxa(key)} taxa)")
    return 1 if errors or (args.strict and warnings) else 0


def cmd_identify(args):
    key = load_checked(args.file)
    if key is None:
        return 1
    observations = {}
    if args.obs_file:
        parsed, diags = parse_observations(read_text(args.obs_file), args.obs_file)
        if diags:
            for diag in diags:
                print(diag.format(args.obs_file), file=sys.stderr)
            return 1
        observations.update(parsed)
    for spec in args.obs:
        parsed, diags = parse_observations(spec, "--obs")
        if diags:
            print(f"dichokey: bad --obs {spec!r}; expected CHARACTER=STATE", file=sys.stderr)
            return 1
        observations.update(parsed)  # command-line values override the file
    result = identify(key, observations)
    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
    else:
        print(describe(result))
    return 0 if result.status == "identified" else 1


def cmd_batch(args):
    key = load_checked(args.file)
    if key is None:
        return 1
    results, diags = run_batch(key, read_text(args.specimens))
    for diag in diags:
        print(diag.format(args.specimens))
    counts = Counter()
    taxa = Counter()
    for specimen, result in results:
        counts[result.status] += 1
        if result.status == "identified":
            label = result.taxon
            taxa[result.taxon] += 1
        elif result.status == "undetermined":
            label = f"undetermined (needs {result.needs})"
        else:
            label = f"no match at couplet {result.couplet}"
        print(f"{specimen}\t{label}")
    print(f"{len(results)} specimens: {counts['identified']} identified, "
          f"{counts['undetermined']} undetermined, {counts['no-match']} no match, "
          f"{len(diags)} row error(s)")
    for name, count in taxa.most_common():
        print(f"  {name}: {count}")
    gate_failed = diags or (args.strict and (counts["undetermined"] or counts["no-match"]))
    return 1 if gate_failed else 0


def cmd_tree(args):
    key = load_checked(args.file)
    if key is None:
        return 1
    print("\n".join(tree_lines(key)))
    return 0


def cmd_taxa(args):
    key = load_checked(args.file)
    if key is None:
        return 1
    for taxon, path in leaves_from(key, 1):
        signature = "; ".join(f"{step.char} = {step.state}" for step in path)
        print(f"{taxon}: {signature}")
    return 0


def cmd_separate(args):
    key = load_checked(args.file)
    if key is None:
        return 1
    info = separate(key, args.taxon_a, args.taxon_b)
    print(f"{info['taxon_a']} and {info['taxon_b']} part at couplet {info['couplet']} "
          f"by {info['char']}: {info['state_a']} vs {info['state_b']}")
    return 0


def cmd_stats(args):
    key = load_checked(args.file)
    if key is None:
        return 1
    data = stats(key)
    if args.json:
        print(json.dumps(data, indent=2))
    else:
        print(f"couplets: {data['couplets']}")
        print(f"leads: {data['leads']}")
        print(f"taxa: {data['taxa']}")
        print(f"max depth: {data['max_depth']}")
        print(f"mean depth: {data['mean_depth']}")
        print("couplets per character:")
        for char, count in data["characters"].items():
            print(f"  {char}: {count}")
    return 0


def cmd_dot(args):
    key = load_checked(args.file)
    if key is None:
        return 1
    text = dot_text(key)
    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(
        prog="dichokey.py",
        description="Offline dichotomous identification-key toolkit.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("lint", help="check a key file; exit 1 on errors")
    p.add_argument("file")
    p.add_argument("--strict", action="store_true", help="also fail on warnings")
    p.set_defaults(func=cmd_lint)

    p = sub.add_parser("identify", help="identify a specimen from observed characters")
    p.add_argument("file")
    p.add_argument("--obs", action="append", default=[], metavar="CHAR=STATE",
                   help="an observed character state (repeatable)")
    p.add_argument("--obs-file", metavar="FILE", help="file of 'character = state' lines")
    p.add_argument("--json", action="store_true", help="print the result as JSON")
    p.set_defaults(func=cmd_identify)

    p = sub.add_parser("batch", help="identify every row of a specimens CSV")
    p.add_argument("file")
    p.add_argument("specimens", help="CSV with a 'specimen' column and one column per character")
    p.add_argument("--strict", action="store_true", help="also fail on undetermined or no-match rows")
    p.set_defaults(func=cmd_batch)

    p = sub.add_parser("tree", help="print the key as an indented tree")
    p.add_argument("file")
    p.set_defaults(func=cmd_tree)

    p = sub.add_parser("taxa", help="list every taxon with the characters that lead to it")
    p.add_argument("file")
    p.set_defaults(func=cmd_taxa)

    p = sub.add_parser("separate", help="show where two taxa's routes part")
    p.add_argument("file")
    p.add_argument("taxon_a")
    p.add_argument("taxon_b")
    p.set_defaults(func=cmd_separate)

    p = sub.add_parser("stats", help="depth, size and character usage of a key")
    p.add_argument("file")
    p.add_argument("--json", action="store_true", help="print the stats as JSON")
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("dot", help="export the key as Graphviz DOT")
    p.add_argument("file")
    p.add_argument("-o", "--output", metavar="FILE", help="write to FILE instead of stdout")
    p.set_defaults(func=cmd_dot)

    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (OSError, ValueError) as error:
        print(f"dichokey: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
