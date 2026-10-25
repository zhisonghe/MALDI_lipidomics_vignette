"""Sum-composition reference of common lipid classes, generated from headgroup formula rules.

Each class maps (total carbons C, double bonds D) of its acyl/alkyl chains (for sphingolipids: sphingoid base
plus N-acyl chain) to a neutral formula. Rules were checked against known species (e.g. PC 34:1 = C42H82NO8P,
SM 34:1;O2 = C39H79N2O6P, HexCer 42:2;O2 = C48H91NO8, PE O-38:5 (= PE P-38:4) = C43H78NO7P, TG 52:2 = C55H102O6).
"""
import re
import pandas as pd

MASS = {"C": 12.0, "H": 1.00782503223, "N": 14.00307400443, "O": 15.99491461957, "P": 30.97376199842}
# class: (name template, C offset, H offset, extra atoms, C range, D range)
CLASSES = {
    "PC":     ("PC {c}:{d}",       8, 16, "NO8P", (26, 46), (0, 12)),
    "PC O-":  ("PC O-{c}:{d}",     8, 18, "NO7P", (28, 44), (0, 10)),
    "PE":     ("PE {c}:{d}",       5, 10, "NO8P", (26, 46), (0, 12)),
    "PE O-":  ("PE O-{c}:{d}",     5, 12, "NO7P", (28, 44), (0, 10)),
    "PS":     ("PS {c}:{d}",       6, 10, "NO10P", (28, 46), (0, 12)),
    "PI":     ("PI {c}:{d}",       9, 15, "O13P", (28, 46), (0, 12)),
    "PG":     ("PG {c}:{d}",       6, 11, "O10P", (28, 46), (0, 12)),
    "PA":     ("PA {c}:{d}",       3, 5, "O8P", (28, 46), (0, 12)),
    "LPC":    ("LPC {c}:{d}",      8, 18, "NO7P", (12, 24), (0, 6)),
    "LPE":    ("LPE {c}:{d}",      5, 12, "NO7P", (12, 24), (0, 6)),
    "SM":     ("SM {c}:{d};O2",    5, 13, "N2O6P", (28, 46), (0, 4)),
    "Cer":    ("Cer {c}:{d};O2",   0, 1, "NO3", (30, 46), (0, 4)),
    "HexCer": ("HexCer {c}:{d};O2", 6, 11, "NO8", (30, 46), (0, 4)),
    "DG":     ("DG {c}:{d}",       3, 4, "O5", (26, 44), (0, 10)),
    "TG":     ("TG {c}:{d}",       3, 2, "O6", (40, 62), (0, 12)),
    "CE":     ("CE {c}:{d}",      27, 44, "O2", (14, 24), (0, 6)),
}


def mono(formula):
    return sum(MASS[e] * int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula))


def lipid_reference():
    rows = []
    for cls, (tmpl, cc, hh, extra, (c0, c1), (d0, d1)) in CLASSES.items():
        for c in range(c0, c1 + 1):
            for d in range(d0, d1 + 1):
                h = 2 * c - 2 * d + hh
                if d > c // 2 or h <= 0:
                    continue
                formula = f"C{c + cc}H{h}{extra}"
                rows.append((tmpl.format(c=c, d=d), cls, formula, mono(formula)))
    return pd.DataFrame(rows, columns=["name", "class", "formula", "exact_mass"])


if __name__ == "__main__":
    ref = lipid_reference()
    check = {"PC 34:1": 759.5778, "SM 34:1;O2": 702.5676, "HexCer 42:2;O2": 809.6745, "PE O-38:5": 751.5516,
             "TG 52:2": 858.7676}
    for n, m in check.items():
        got = ref.loc[ref.name == n, "exact_mass"].iat[0]
        print(f"{n}: {got:.4f} (expected {m}) {'ok' if abs(got - m) < 2e-4 else 'MISMATCH'}")
    print(len(ref), "sum compositions")
