#!/usr/bin/env python3
"""Inject verified DOIs (and a few authoritative metadata corrections) into
references.bib. DOIs were obtained from the CrossRef REST API and each was
verified by matching journal + volume + first page against the entry.

Two groups:
  * CROSSREF  -- fetched and field-matched individually.
  * PATTERN   -- APS (10.1103/<Journal>.<vol>.<firstpage>) and PTP
                 (10.1143/PTP.<vol>.<firstpage>) legacy DOIs, constructed from
                 the deterministic scheme after it was verified against 7 APS
                 and 4 PTP actual CrossRef lookups (all exact). Flagged so they
                 can be double-checked when the web quota resets.
"""
import re, sys

DOIS = {
    # --- CrossRef-verified ---
    'ref015': '10.1103/RevModPhys.2.305',
    'ref027': '10.1007/978-3-662-11761-3',
    'ref036': '10.3367/ufnr.0106.197201a.0003',
    'ref048': '10.1007/BF02750104',
    'ref049': '10.1143/PTP.11.143',
    'ref053': '10.1103/RevModPhys.34.829',
    'ref056': '10.1103/RevModPhys.24.249',
    'ref057': '10.1002/sapm1952311287',
    'ref058': '10.1063/1.1665333',
    'ref059': '10.1098/rspa.1951.0110',
    'ref060': '10.1016/S0031-8914(57)95547-7',
    'ref065': '10.1098/rspa.1953.0109',
    'ref067': '10.1007/BF01336904',
    'ref069': '10.1063/1.1724228',
    'ref070': '10.1090/qam/60649',
    'ref071': '10.1016/0003-4916(68)90146-2',
    'ref073': '10.1098/rspa.1951.0026',
    'ref074': '10.1103/PhysRev.93.318',
    'ref075': '10.1063/1.1664592',
    'ref078': '10.1063/1.1672352',
    'ref081': '10.1103/PhysRev.110.815',
    'ref082': '10.1143/PTP.20.798',
    'ref083': '10.1007/BF03157469',
    'ref085': '10.1002/sapm1957361157',
    'ref086': '10.1016/0375-9474(68)90895-6',
    'ref087': '10.1016/0003-4916(66)90040-6',
    'ref090': '10.1103/PhysRev.61.186',
    'ref094': '10.1007/BF02859841',
    'ref095': '10.1007/BF02724914',
    'ref096': '10.1002/sapm1958371215',
    'ref098': '10.1063/1.1704115',
    'ref099': '10.1143/PTP.13.405',
    'ref103': '10.1007/BF02771400',
    'ref104': '10.1016/0375-9474(67)90163-7',
    'ref107': '10.1088/0370-1298/70/12/111',
    'ref108': '10.1139/p64-101',
    'ref109': '10.2307/2371276',
    'ref111': '10.1063/1.1665855',
    'ref114': '10.1016/S0092-640X(71)80020-7',
    'ref116': '10.1139/p64-036',
    'ref117': '10.1139/p52-024',
    'ref119': '10.1139/p57-038',
    'ref121': '10.1103/PhysRev.107.186',
    'ref122': '10.1143/PTPS.26.64',
    'ref126': '10.1143/PTPS.13.1',
    'ref129': '10.1016/0550-3213(70)90140-9',
    'ref131': '10.1063/1.1665626',
    'ref133': '10.1088/0305-4470/11/3/009',
    'ref134': '10.1007/BF00420700',
    'ref135': '10.1063/1.522543',
    'ref136': '10.1088/0305-4470/9/8/009',
    'ref140': '10.1063/1.522427',
    'ref144': '10.1088/0305-4470/7/16/004',
    # --- pattern-derived (APS / PTP), not individually fetched ---
    'ref047': '10.1103/RevModPhys.28.432',
    'ref072': '10.1103/PhysRev.111.194',
    'ref091': '10.1103/PhysRev.62.438',
    'ref092': '10.1103/PhysRev.63.367',
    'ref093': '10.1103/PhysRev.84.910',
    'ref120': '10.1143/PTP.14.589',
    'ref125': '10.1143/PTP.13.540',
    'ref130': '10.1143/PTP.8.431',
}

# volume added where the entry lacked it (PTP Supplement)
VOL_ADD = {'ref122': '26', 'ref126': '13'}

# authoritative metadata corrections (old -> new), verified against matched DOI
YEAR_FIX = {'ref053': ('1952', '1962'),
            'ref096': ('1957', '1958'),
            'ref109': ('1939', '1941')}
PAGE_FIX = {'ref109': ('12--17', '57--63')}

def main():
    path = 'references.bib'
    with open(path, encoding='utf-8') as f:
        text = f.read()
    # split into entries keeping the leading @
    parts = re.split(r'(?m)^(?=@)', text)
    out = []
    changed = 0
    for part in parts:
        m = re.match(r'@\w+\{(ref\d+),', part)
        if not m:
            out.append(part); continue
        key = m.group(1)
        if key not in DOIS:
            out.append(part); continue
        # measure alignment from the keywords line
        km = re.search(r'(?m)^(\s*)keywords(\s*)=\s', part)
        indent = km.group(1)
        width = len('keywords') + len(km.group(2))   # chars before '='
        def field(name, val):
            return f"{indent}{name.ljust(width)}= {{{val}}},\n"
        # year / page corrections
        if key in YEAR_FIX:
            old, new = YEAR_FIX[key]
            part, n = re.subn(r'(year\s*=\s*\{)'+old+r'(\})', r'\g<1>'+new+r'\2', part)
            assert n == 1, f'year fix failed for {key}'
        if key in PAGE_FIX:
            old, new = PAGE_FIX[key]
            part, n = re.subn(re.escape('{'+old+'}'), '{'+new+'}', part)
            assert n == 1, f'page fix failed for {key}'
        # build the block to insert before the keywords line
        block = ''
        if key in VOL_ADD and not re.search(r'(?m)^\s*volume\s*=', part):
            block += field('volume', VOL_ADD[key])
        block += field('doi', DOIS[key])
        # insert immediately before the keywords line
        part = re.sub(r'(?m)^(\s*keywords\s*=\s)', block + r'\1', part, count=1)
        changed += 1
        out.append(part)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(''.join(out))
    print(f'entries updated: {changed} / {len(DOIS)} in map')

if __name__ == '__main__':
    main()
