#!/usr/bin/env python3
"""Flag the keys of a SWIFT parameter file that SWIFT would not read.

The valid keys are collected statically, without running SWIFT, from:

1. the string literals passed to ``parser_get_param_*``,
   ``parser_get_opt_param_*`` and ``parser_does_param_exist`` in ``src/``
   (adjacent literals, ``#define`` string macros, ternaries, format strings
   such as ``"%s:max_age"`` and keys built into a buffer with
   ``sprintf``/``strcpy`` in the same file are resolved), and
2. every key in ``examples/parameter_example.yml``.

A key that is in neither set is reported as UNKNOWN. A key found in a table
of retired keys in ``src/`` (an array whose name contains ``retired``, with
``{"old_key", "replacement"}`` entries) is reported as RETIRED, together
with its replacement. Keys that depend on a runtime value and cannot be
resolved are listed at the end, so the gap in the check is visible.

The key set is the union over every module and every ``#ifdef`` branch of
``src/``. A key read only by a module that the build does not use is
therefore not reported.

Examples
--------
Check one file::

    tools/check_param_keys.py params.yml

Check every parameter file under ``examples/`` and fail on any finding::

    tools/check_param_keys.py --all-examples --strict
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Pattern, Sequence, Set, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent

PARSER_CALL = re.compile(
    r"\bparser_(?:get_param|get_opt_param|does_param_exist)[a-z_0-9]*\s*\("
)
STRING_LITERAL = re.compile(r'"(?:[^"\\\n]|\\.)*"')
IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")
FORMAT_SPEC = re.compile(
    r"%(?:%|[-+ #0]*[0-9*]*(?:\.[0-9*]+)?(?:hh|h|ll|l|z|j|t)?[a-zA-Z])"
)
RETIRED_TABLE = re.compile(r"\b\w*retired\w*\s*\[\s*\]\s*=\s*\{")
MIN_PATTERN_LITERAL_CHARS = 3


@dataclass
class KeySet:
    """The keys SWIFT may read.

    Attributes
    ----------
    literals : set of str
        Keys given as plain strings.
    patterns : list of (regex, str)
        Keys given as format strings, as a compiled regex and the source text.
    """

    literals: Set[str] = field(default_factory=set)
    patterns: List[Tuple[Pattern[str], str]] = field(default_factory=list)

    def __contains__(self, key: str) -> bool:
        """Return whether ``key`` is a literal key or matches a pattern."""
        if key in self.literals:
            return True
        return any(rx.fullmatch(key) for rx, _ in self.patterns)

    def add(self, text: str) -> bool:
        """Add a key or a format-string key.

        Parameters
        ----------
        text : str
            The key, possibly with printf conversions.

        Returns
        -------
        bool
            False if ``text`` is too unspecific to use as a pattern.
        """
        if "%" not in text:
            self.literals.add(text)
            return True
        rx = format_to_regex(text)
        if rx is None:
            return False
        if all(src != text for _, src in self.patterns):
            self.patterns.append((rx, text))
        return True


@dataclass
class Unresolved:
    """A parser call whose key could not be resolved statically.

    Attributes
    ----------
    where : str
        ``file:line`` of the call.
    expr : str
        The key expression as written.
    """

    where: str
    expr: str


@dataclass
class Finding:
    """One problem in a parameter file.

    Attributes
    ----------
    path : Path
        The parameter file.
    line : int
        1-based line of the key.
    key : str
        The full ``Section:name`` key.
    kind : str
        ``"RETIRED"`` or ``"UNKNOWN"``.
    advice : str
        The replacement for a retired key, or the closest valid key for an
        unknown one (may be empty).
    """

    path: Path
    line: int
    key: str
    kind: str
    advice: str


def format_to_regex(fmt: str) -> Optional[Pattern[str]]:
    """Turn a printf format string into a regex for the keys it can build.

    Parameters
    ----------
    fmt : str
        A key with printf conversions, such as ``"%s:max_age"``.

    Returns
    -------
    re.Pattern or None
        The regex, or None if the format has fewer than three fixed
        characters and so would match almost any key.
    """
    out: List[str] = []
    fixed = 0
    pos = 0
    for m in FORMAT_SPEC.finditer(fmt):
        lit = fmt[pos : m.start()]
        out.append(re.escape(lit))
        fixed += len(lit)
        spec = m.group(0)
        if spec == "%%":
            out.append("%")
            fixed += 1
        elif spec[-1] in "diuxXo":
            out.append(r"-?\d+")
        else:
            out.append(r"[^:\s]+")
        pos = m.end()
    tail = fmt[pos:]
    out.append(re.escape(tail))
    fixed += len(tail)
    if fixed < MIN_PATTERN_LITERAL_CHARS:
        return None
    return re.compile("".join(out))


def strip_comments(text: str) -> str:
    """Blank out C comments, keeping strings, character literals and lines.

    Parameters
    ----------
    text : str
        C source.

    Returns
    -------
    str
        The source with every comment replaced by spaces, so offsets and line
        numbers are unchanged.
    """
    out: List[str] = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        two = text[i : i + 2]
        if two == "/*":
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j]))
            i = j
        elif two == "//":
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
                if j < n and text[j] == "\n":
                    break
            out.append(text[i : j + 1])
            i = j + 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def balanced_args(text: str, open_pos: int) -> Optional[Tuple[List[str], int]]:
    """Split the argument list that starts at the parenthesis at ``open_pos``.

    Parameters
    ----------
    text : str
        Source with comments removed.
    open_pos : int
        Index of the opening parenthesis.

    Returns
    -------
    tuple or None
        The top-level arguments as written and the index of the closing
        parenthesis, or None if the parentheses do not balance.
    """
    depth = 0
    args: List[str] = []
    start = open_pos + 1
    i = open_pos
    n = len(text)
    while i < n:
        c = text[i]
        if c in "\"'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            i = j
        elif c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
            if depth == 0:
                args.append(text[start:i].strip())
                return args, i
        elif c == "," and depth == 1:
            args.append(text[start:i].strip())
            start = i + 1
        i += 1
    return None


def unescape_c(lit: str) -> str:
    """Return the text of a C string literal without its quotes.

    Parameters
    ----------
    lit : str
        A literal such as ``'"a\\\\n"'``.

    Returns
    -------
    str
        The contents, with ``\\\\`` and ``\\"`` resolved.
    """
    return lit[1:-1].replace('\\"', '"').replace("\\\\", "\\")


def collect_macros(text: str) -> Dict[str, str]:
    """Collect object-like ``#define`` macros whose value is a string.

    Parameters
    ----------
    text : str
        C source with comments removed.

    Returns
    -------
    dict
        Macro name to the concatenated string it expands to.
    """
    macros: Dict[str, str] = {}
    joined = text.replace("\\\n", " ")
    for m in re.finditer(
        r"^[ \t]*#[ \t]*define[ \t]+([A-Za-z_]\w*)[ \t]+(.+)$", joined, re.M
    ):
        value = literal_concat(m.group(2), macros)
        if value is not None:
            macros[m.group(1)] = value
    return macros


def literal_concat(expr: str, macros: Dict[str, str]) -> Optional[str]:
    """Evaluate an expression made only of string literals and string macros.

    Parameters
    ----------
    expr : str
        The expression as written, e.g. ``'"Section:" KEY'``.
    macros : dict
        Known string macros.

    Returns
    -------
    str or None
        The concatenated string, or None if ``expr`` holds anything else.
    """
    pos = 0
    parts: List[str] = []
    expr = expr.strip()
    if not expr:
        return None
    while pos < len(expr):
        if expr[pos].isspace():
            pos += 1
            continue
        m = STRING_LITERAL.match(expr, pos)
        if m:
            parts.append(unescape_c(m.group(0)))
            pos = m.end()
            continue
        m = IDENTIFIER.match(expr, pos)
        if m and m.group(0) in macros:
            parts.append(macros[m.group(0)])
            pos = m.end()
            continue
        return None
    return "".join(parts)


def resolve_variable(
    name: str, text: str, macros: Dict[str, str], need_colon: bool = True
) -> List[str]:
    """Find the key strings a buffer or pointer variable can hold.

    Looks, in the same file, for ``sprintf``/``snprintf``/``strcpy`` into the
    variable, for plain assignments and for initialiser lists.

    Parameters
    ----------
    name : str
        The variable name.
    text : str
        Source of the file with comments removed.
    macros : dict
        String macros of the file.
    need_colon : bool
        Keep only strings that hold a colon, i.e. full ``Section:name`` keys.

    Returns
    -------
    list of str
        Candidate keys (may hold printf conversions).
    """
    found: List[str] = []
    esc = re.escape(name)
    for m in re.finditer(r"\b(sn?printf|check_snprintf|strcpy|strncpy)\s*\(", text):
        got = balanced_args(text, m.end() - 1)
        if got is None:
            continue
        args, _ = got
        if not args or args[0].lstrip("&*") != name:
            continue
        idx = 2 if m.group(1).endswith("snprintf") else 1
        if len(args) > idx:
            lit = literal_concat(args[idx], macros)
            if lit is not None:
                found.append(lit)
    for m in re.finditer(rf"\b{esc}\b\s*(?:\[[^\]]*\])*\s*=\s*([^;]+);", text):
        found += [unescape_c(s) for s in STRING_LITERAL.findall(m.group(1))]
    return [s for s in found if ":" in s or not need_colon]


def resolve_key_expr(
    expr: str, text: str, macros: Dict[str, str]
) -> Optional[List[str]]:
    """Resolve the key argument of a parser call.

    Parameters
    ----------
    expr : str
        The second argument as written.
    text : str
        Source of the file with comments removed.
    macros : dict
        String macros of the file.

    Returns
    -------
    list of str or None
        Every key (or format string) the expression can take, or None if it
        cannot be resolved.
    """
    lit = literal_concat(expr, macros)
    if lit is not None:
        return [lit]
    yml = re.fullmatch(r"YML_NAME\s*\((.*)\)", expr, re.S)
    if yml:
        inner = literal_concat(yml.group(1), macros)
        names = (
            [inner]
            if inner is not None
            else resolve_variable(yml.group(1).strip(), text, macros, need_colon=False)
        )
        if names:
            return [f"Lightcone%d:{n}" for n in names] + [
                f"LightconeCommon:{n}" for n in names
            ]
    if "?" in expr and STRING_LITERAL.search(expr):
        parts = [unescape_c(s) for s in STRING_LITERAL.findall(expr)]
        parts = [p for p in parts if ":" in p]
        if parts:
            return parts
    base = re.match(r"^[&*]?\s*([A-Za-z_]\w*)\s*(?:\[[^\]]*\])*$", expr)
    if base:
        got = resolve_variable(base.group(1), text, macros)
        if got:
            return got
    return None


def scan_source(src: Path) -> Tuple[KeySet, List[Unresolved], Dict[str, str]]:
    """Collect every key SWIFT reads and every retired key from ``src``.

    Parameters
    ----------
    src : Path
        The ``src`` directory (searched recursively for ``.c`` and ``.h``).

    Returns
    -------
    tuple
        The valid keys, the parser calls that could not be resolved, and the
        retired keys mapped to their replacement text.
    """
    keys = KeySet()
    unresolved: List[Unresolved] = []
    retired: Dict[str, str] = {}
    for path in sorted(src.rglob("*")):
        if path.suffix not in (".c", ".h") or path.name in ("parser.c", "parser.h"):
            continue
        raw = path.read_text(errors="replace")
        text = strip_comments(raw)
        macros = collect_macros(text)
        rel = path.relative_to(src.parent)
        for m in RETIRED_TABLE.finditer(text):
            retired.update(parse_retired_table(text, m.end() - 1, macros))
        for m in PARSER_CALL.finditer(text):
            got = balanced_args(text, m.end() - 1)
            if got is None or len(got[0]) < 2:
                continue
            expr = got[0][1]
            line = text.count("\n", 0, m.start()) + 1
            resolved = resolve_key_expr(expr, text, macros)
            ok = resolved is not None and all(keys.add(k) for k in resolved)
            if not ok:
                unresolved.append(Unresolved(f"{rel}:{line}", " ".join(expr.split())))
    return keys, unresolved, retired


def parse_retired_table(
    text: str, brace_pos: int, macros: Dict[str, str]
) -> Dict[str, str]:
    """Read ``{"old", "replacement"}`` entries from a C array initialiser.

    Parameters
    ----------
    text : str
        Source with comments removed.
    brace_pos : int
        Index of the opening brace of the initialiser.
    macros : dict
        String macros of the file.

    Returns
    -------
    dict
        Retired key to replacement text.
    """
    got = balanced_args(text, brace_pos)
    if got is None:
        return {}
    _, end = got
    body = text[brace_pos + 1 : end]
    table: Dict[str, str] = {}
    for entry in re.finditer(r"\{([^{}]*)\}", body):
        fields = balanced_args("(" + entry.group(1) + ")", 0)
        if fields is None or len(fields[0]) < 2:
            continue
        old = literal_concat(fields[0][0], macros)
        new = literal_concat(fields[0][1], macros)
        if old is not None and new is not None:
            table[old] = new
    return table


def read_param_keys(path: Path) -> List[Tuple[str, int]]:
    """List the keys of a SWIFT parameter file as the SWIFT parser reads them.

    Parameters
    ----------
    path : Path
        The parameter file.

    Returns
    -------
    list of (str, int)
        Each ``Section:name`` key with its 1-based line. A key at column 0
        with a value is standalone and has no section.
    """
    keys: List[Tuple[str, int]] = []
    section = ""
    for number, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if raw.startswith("#"):
            continue
        body = raw.split("#", 1)[0].rstrip()
        if not body.strip() or ":" not in body:
            continue
        name, _, value = body.partition(":")
        if raw[0] in " \t":
            keys.append((f"{section}{name.strip()}", number))
        elif value.strip():
            keys.append((name.strip(), number))
            section = ""
        else:
            section = name.strip() + ":"
    return keys


def is_param_file(path: Path) -> bool:
    """Return whether ``path`` is a SWIFT parameter file.

    Parameters
    ----------
    path : Path
        A YAML file.

    Returns
    -------
    bool
        True if the file has an ``InternalUnitSystem`` section, which every
        runnable parameter file needs.
    """
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return False
    return re.search(r"^InternalUnitSystem:", text, re.M) is not None


def check_file(path: Path, valid: KeySet, retired: Dict[str, str]) -> List[Finding]:
    """Check one parameter file.

    Parameters
    ----------
    path : Path
        The parameter file.
    valid : KeySet
        The keys SWIFT reads.
    retired : dict
        Retired key to replacement text.

    Returns
    -------
    list of Finding
        One entry per retired or unknown key, in file order.
    """
    findings: List[Finding] = []
    pool = sorted(valid.literals)
    for key, line in read_param_keys(path):
        if key in retired:
            findings.append(Finding(path, line, key, "RETIRED", retired[key]))
        elif key not in valid:
            close = difflib.get_close_matches(key, pool, n=1, cutoff=0.85)
            findings.append(
                Finding(path, line, key, "UNKNOWN", close[0] if close else "")
            )
    return findings


def find_param_files(root: Path) -> List[Path]:
    """Find every SWIFT parameter file under ``root``.

    Parameters
    ----------
    root : Path
        Directory to search.

    Returns
    -------
    list of Path
        YAML files with an ``InternalUnitSystem`` section, except
        ``parameter_example.yml`` (the reference, not an input) and the
        ``used_parameters*.yml``/``unused_parameters*.yml`` files that SWIFT
        writes (outputs, not inputs).
    """
    files = [
        p
        for p in sorted(root.rglob("*"))
        if p.suffix in (".yml", ".yaml")
        and p.is_file()
        and p.name != "parameter_example.yml"
        and not p.name.startswith(("used_parameters", "unused_parameters"))
        and is_param_file(p)
    ]
    return files


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the check.

    Parameters
    ----------
    argv : sequence of str, optional
        Command-line arguments; ``sys.argv[1:]`` when omitted.

    Returns
    -------
    int
        0, or 1 with ``--strict`` when anything was found.
    """
    ap = argparse.ArgumentParser(
        description="List the keys of SWIFT parameter files that SWIFT would "
        "not read. Does not run SWIFT."
    )
    ap.add_argument("files", nargs="*", type=Path, help="parameter files to check")
    ap.add_argument(
        "--all-examples",
        action="store_true",
        help="check every parameter file under examples/",
    )
    ap.add_argument(
        "--root",
        type=Path,
        default=REPO_ROOT,
        help="SWIFT checkout (default: this one)",
    )
    ap.add_argument(
        "--strict", action="store_true", help="exit 1 when any key is flagged"
    )
    ap.add_argument(
        "--summary-only",
        action="store_true",
        help="print the counts and the unresolved patterns, not each finding",
    )
    args = ap.parse_args(argv)

    root: Path = args.root
    valid, unresolved, retired = scan_source(root / "src")
    example = root / "examples" / "parameter_example.yml"
    if example.is_file():
        for key, _ in read_param_keys(example):
            valid.literals.add(key)

    files: List[Path] = list(args.files)
    if args.all_examples:
        files += find_param_files(root / "examples")
    if not files:
        ap.error("give at least one parameter file or --all-examples")

    findings: List[Finding] = []
    for path in files:
        found = check_file(path, valid, retired)
        findings += found
        if not args.summary_only:
            for f in found:
                if f.kind == "RETIRED":
                    note = f"use: {f.advice}"
                else:
                    note = (
                        f"closest valid key: {f.advice}" if f.advice else "no close key"
                    )
                print(f"{f.path}:{f.line}: {f.kind} {f.key} ({note})")

    n_retired = sum(f.kind == "RETIRED" for f in findings)
    n_unknown = sum(f.kind == "UNKNOWN" for f in findings)
    files_hit = len({f.path for f in findings})

    by_key: Dict[str, Set[Path]] = defaultdict(set)
    for f in findings:
        by_key[f"{f.kind} {f.key}"].add(f.path)
    if len(files) > 1 and by_key:
        print("\nDistinct findings (kind key: number of files):")
        for k in sorted(by_key):
            print(f"  {k}: {len(by_key[k])}")

    print(
        f"\nKnown keys: {len(valid.literals)} literal, {len(valid.patterns)} "
        f"patterns. Retired-key table: {len(retired)} entries."
    )
    print(
        f"Unresolved parser calls: {len(unresolved)} (keys not checked against these):"
    )
    for u in unresolved:
        print(f"  {u.where}: {u.expr}")
    print(
        f"\nSummary: {len(files)} file(s) checked, {files_hit} with findings; "
        f"{n_retired} RETIRED, {n_unknown} UNKNOWN."
    )
    return 1 if (args.strict and findings) else 0


if __name__ == "__main__":
    sys.exit(main())
