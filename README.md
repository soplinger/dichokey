# Dichokey

**Offline dichotomous identification keys, written in plain text and checked before you trust them.**

Dichokey takes a plain-text identification key (the "if leaf is needle, go to 2, otherwise go to 4" kind used in field guides), checks it for dead ends, cycles and unreachable couplets, identifies specimens from what you observed, and tells you where two look-alike taxa part company. It is one Python file with no dependencies and no network access.

## Features

- **Plain-text key format**: one couplet lead per line, `1a leaf = needle -> 2`, with a `key:` title and `#` comments.
- **Line-numbered lint**: malformed lines, duplicate leads, incomplete couplets, mixed characters, identical states, undefined targets, cycles and unreachable couplets are reported as `file:line: error|warning`.
- **CI gate**: `lint` exits 1 on errors, and `--strict` also fails on warnings.
- **Identify from observations**: `identify` walks the key with `--obs CHARACTER=STATE` values or an `--obs-file`, and prints the path it took.
- **Missing-character guidance**: when a needed character is not supplied, it names that character and lists every taxon still possible from that couplet.
- **No-match detection**: an observed state that fits neither lead is reported with the states that were expected.
- **Case- and space-insensitive matching**: `Leaf`, `  NEEDLE ` and `needle` are the same value.
- **JSON output**: `identify --json` gives status, taxon, path, missing character and candidates for scripts.
- **Batch identification**: `batch` runs a specimens CSV, reports row-level errors with line numbers, and tallies taxa.
- **Taxon signatures**: `taxa` lists every taxon with the character states that lead to it.
- **Separating two taxa**: `separate` shows the couplet and characters where two taxa's routes part.
- **Indented tree view**: `tree` prints the whole key, marking repeated couplets as `(see couplet N)`.
- **Key statistics**: `stats` reports couplets, leads, taxa, maximum and mean path depth, and how many couplets each character decides (`--json` available).
- **Graphviz export**: `dot` writes the key as a DOT graph to stdout or a file.

## Quick start

The repository includes an example park key and a specimens file.

```console
$ python3 dichokey.py lint examples/trees.key
examples/trees.key: OK (4 couplets, 5 taxa)

$ python3 dichokey.py identify examples/trees.key --obs leaf=needle --obs bundles=yes
Undetermined at couplet 3: need a value for 'cone'.
Still possible: Juniperus communis, Pinus sylvestris
Path:
  couplet 1: leaf = needle -> 2
  couplet 2: bundles = yes -> 3

$ python3 dichokey.py batch examples/trees.key examples/specimens.csv
spec-01	Pinus sylvestris
spec-02	Taxus baccata
spec-03	Quercus robur
spec-04	no match at couplet 4
spec-05	undetermined (needs cone)
5 specimens: 3 identified, 1 undetermined, 1 no match, 0 row error(s)
  Pinus sylvestris: 1
  Taxus baccata: 1
  Quercus robur: 1

$ python3 dichokey.py separate examples/trees.key "Pinus sylvestris" "Quercus robur"
Pinus sylvestris and Quercus robur part at couplet 1 by leaf: needle vs broad
```

A key with mistakes gets line-numbered diagnostics:

```console
$ python3 dichokey.py lint examples/broken.key
examples/broken.key:6: error: cannot parse '3a cone -> Pinus'; expected e.g. '1a leaf = needle -> 2'
examples/broken.key:3: error: couplet 1 has two leads with the same state
examples/broken.key:3: error: lead 1b points to undefined couplet 3
examples/broken.key:4: error: lead 2a points to undefined couplet 9
examples/broken.key:7: error: lead 4a returns to couplet 4, which makes a cycle
examples/broken.key:7: warning: couplet 4 is never reached from couplet 1
examples/broken.key: 5 error(s), 1 warning(s)
```

### Writing a key

```text
key: Common trees of the example park
# <number><a|b> <character> = <state> -> <next couplet or taxon>
1a leaf = needle -> 2
1b leaf = broad -> 4
2a bundles = yes -> 3
2b bundles = no -> Taxus baccata
3a cone = woody -> Pinus sylvestris
3b cone = berry -> Juniperus communis
4a lobes = deep -> Quercus robur
4b lobes = shallow -> Acer campestre
```

Every couplet starts at `1` and has an `a` and a `b` lead that test the same character. Each lead points either to another couplet number or to a taxon name.

### Specimen files

`batch` reads a CSV with a `specimen` column and one column per character. Empty cells mean "not observed".

```csv
specimen,leaf,bundles,cone,lobes
spec-01,needle,yes,woody,
spec-05,needle,yes,,
```

## Commands

| Command | Purpose | Exit codes |
| --- | --- | --- |
| `lint FILE [--strict]` | Check a key | 0 clean, 1 errors (or warnings with `--strict`) |
| `identify FILE [--obs C=S]... [--obs-file F] [--json]` | Identify one specimen | 0 identified, 1 otherwise |
| `batch FILE SPECIMENS.csv [--strict]` | Identify a whole CSV of specimens | 0 ok, 1 row errors (or unresolved rows with `--strict`) |
| `tree FILE` | Print the key as an indented tree | 0 |
| `taxa FILE` | List taxa with the states that lead to them | 0 |
| `separate FILE TAXON_A TAXON_B` | Show where two taxa's routes part | 0 ok, 1 unknown taxon |
| `stats FILE [--json]` | Depth, size and character usage | 0 |
| `dot FILE [-o OUT]` | Export Graphviz DOT | 0 |

Run `python3 dichokey.py COMMAND --help` for details. Command-line `--obs` values override the same characters in `--obs-file`. A key with lint errors is refused by every command except `lint`.

## How to run

Requires Python 3.10 or newer. There are no third-party dependencies.

```console
$ python3 dichokey.py --help
```

Run it from the project directory, or give the path to `dichokey.py`.

## How to test

The test suite uses the standard library `unittest` module, so it runs without installing anything:

```console
$ python3 -m unittest discover -s tests -v
```

If pytest is installed, this also works:

```console
$ python3 -m pytest
```
