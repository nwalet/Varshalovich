#!/usr/bin/env python3
r"""Apply CONVERSION_RULES.md rule 7 (automated equation numbering) to a chapter file.

Finds every display-math block (`\[...\]`, `equation*`, `align*`, `gather*`)
that contains at least one `\tag{n}`, drops the star / promotes `\[...\]` to
`equation`, and replaces each `\tag{n}` with `\label{chap<C>:eq:<N>}` where
`<N>` is a fresh, chapter-wide sequential number in document order.

The source resets its manual numbering `n` back to 1 at each top-level
`\section` (see rule 6), so the same `n` is reused several times per
chapter; a depth-aware scan (tracking `\begin{...}/\end{...}` nesting so
nested `array`/`matrix` row breaks aren't mistaken for top-level
`align`/`gather` row breaks) splits each block into logical rows, adds
`\notag` to rows that had no `\tag` in the source, and records an
(old-scope, old-local-number) -> new-global-number mapping so prose
cross-references ("Eq. (9)", "Eqs. (26)-(27)", ...) can be resolved
separately with tools/build_eqref_map.py or by hand.

Usage: python3 tools/convert_equation_numbering.py Chap1.tex
Writes Chap1.tex in place and Chap1.eqmap.json alongside it (mapping
"scope,local_n" -> global_n, plus the scope boundary line numbers, for the
prose cross-reference pass).
"""
import json
import re
import sys
from pathlib import Path

BEGIN_RE = re.compile(r'^\\begin\{(equation\*|align\*|gather\*)\}\s*$')
END_RE = re.compile(r'^\\end\{(equation\*|align\*|gather\*)\}\s*$')
BRACKET_OPEN_RE = re.compile(r'^\\\[\s*$')
BRACKET_CLOSE_RE = re.compile(r'^\\\]\s*$')
TAG_RE = re.compile(r'\\tag\{(\d+)\}')
BEGIN_ENV_RE = re.compile(r'\\begin\{[a-zA-Z*]+\}')
END_ENV_RE = re.compile(r'\\end\{[a-zA-Z*]+\}')

ENV_RENAME = {
    'equation*': 'equation',
    'align*': 'align',
    'gather*': 'gather',
}


def split_top_level_rows(body):
    """Split a block body into rows at top-level (depth-0) '\\\\', tracking
    \\begin{...}/\\end{...} nesting so nested array/matrix row breaks aren't
    mistaken for top-level align/gather row breaks. Returns list of row
    strings (without the trailing '\\\\')."""
    rows = []
    depth = 0
    i = 0
    row_start = 0
    n = len(body)
    while i < n:
        m = BEGIN_ENV_RE.match(body, i)
        if m:
            depth += 1
            i = m.end()
            continue
        m = END_ENV_RE.match(body, i)
        if m:
            depth -= 1
            i = m.end()
            continue
        if depth == 0 and body[i:i + 2] == '\\\\':
            rows.append(body[row_start:i])
            i += 2
            row_start = i
            continue
        i += 1
    rows.append(body[row_start:])
    return rows


def process_block(kind, body, chapter, global_n, scope_state, mapping):
    """Returns (new_body, global_n)."""
    tags_in_order = TAG_RE.findall(body)

    if not tags_in_order:
        return None, global_n  # unchanged, caller keeps original text

    def record(local_n):
        nonlocal global_n
        local_n = int(local_n)
        if local_n <= scope_state['last_local']:
            scope_state['scope'] += 1
        scope_state['last_local'] = local_n
        global_n += 1
        mapping[f"{scope_state['scope']},{local_n}"] = global_n
        return global_n

    if kind == 'bracket':
        assert len(tags_in_order) == 1, f"bracket block has {len(tags_in_order)} tags, expected 1"
        gn = record(tags_in_order[0])
        new_body = TAG_RE.sub(lambda m: f"\\label{{chap{chapter}:eq:{gn}}}", body, count=1)
        return new_body, global_n

    if kind == 'equation*':
        assert len(tags_in_order) == 1, f"equation* block has {len(tags_in_order)} tags, expected 1"
        gn = record(tags_in_order[0])
        new_body = TAG_RE.sub(lambda m: f"\\label{{chap{chapter}:eq:{gn}}}", body, count=1)
        return new_body, global_n

    # align*/gather*: split into logical rows, tag-or-notag each
    rows = split_top_level_rows(body)
    new_rows = []
    for row in rows:
        row_tags = TAG_RE.findall(row)
        assert len(row_tags) <= 1, f"row has {len(row_tags)} tags: {row!r}"
        if row_tags:
            gn = record(row_tags[0])
            new_row = TAG_RE.sub(lambda m: f"\\label{{chap{chapter}:eq:{gn}}}", row, count=1)
        else:
            if row.strip():
                new_row = row.rstrip('\n') + ' \\notag'
                if row.endswith('\n'):
                    new_row += '\n'
            else:
                new_row = row
        new_rows.append(new_row)
    new_body = '\\\\'.join(new_rows)
    return new_body, global_n


def convert(path, chapter):
    lines = Path(path).read_text().split('\n')
    # re-add trailing newline marker per line except we'll join with \n at the end
    out = []
    mapping = {}
    scope_state = {'scope': 0, 'last_local': 0}
    global_n = 0
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        m_begin = BEGIN_RE.match(line)
        if m_begin:
            env = m_begin.group(1)
            j = i + 1
            while not END_RE.match(lines[j]):
                j += 1
            body = '\n'.join(lines[i + 1:j]) + '\n'
            new_body, global_n = process_block(env, body, chapter, global_n, scope_state, mapping)
            if new_body is None:
                out.extend(lines[i:j + 1])
            else:
                new_env = ENV_RENAME[env]
                out.append(f'\\begin{{{new_env}}}')
                out.extend(new_body.split('\n')[:-1])
                out.append(f'\\end{{{new_env}}}')
            i = j + 1
            continue
        if BRACKET_OPEN_RE.match(line):
            j = i + 1
            while not BRACKET_CLOSE_RE.match(lines[j]):
                j += 1
            body = '\n'.join(lines[i + 1:j]) + '\n'
            new_body, global_n = process_block('bracket', body, chapter, global_n, scope_state, mapping)
            if new_body is None:
                out.extend(lines[i:j + 1])
            else:
                out.append('\\begin{equation}')
                out.extend(new_body.split('\n')[:-1])
                out.append('\\end{equation}')
            i = j + 1
            continue
        out.append(line)
        i += 1

    Path(path).write_text('\n'.join(out))
    mapfile = Path(path).with_suffix('.eqmap.json')
    mapfile.write_text(json.dumps(mapping, indent=2, sort_keys=True))
    print(f"Wrote {path} ({global_n} equations labeled, {scope_state['scope'] + 1} scopes)")
    print(f"Wrote {mapfile}")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    target = sys.argv[1]
    # chapter number from filename, e.g. Chap1.tex -> 1
    m = re.search(r'(\d+)', Path(target).stem)
    if not m:
        print("Could not infer chapter number from filename; expected e.g. Chap1.tex")
        sys.exit(1)
    convert(target, m.group(1))
