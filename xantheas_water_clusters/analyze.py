#!/usr/bin/env python3
"""Body-by-body binding energy, split into HF, CCSD correlation and (T).

    python3 analyze.py output_w6_prism__ri-ccsdt__none.json [--ref output_w1ref__ri-ccsdt.json]

Reads the `levels` of an MBE output. Each column is expanded on its own, so a
row is that body order's contribution to that part of the binding energy, in
kcal/mol, with its share of the column's total in parentheses -- the layout of
Tables 1 and 2 of Nguyen and Xantheas (2026).

The one-body row is sum_I E_I - N E(1) and needs --ref, the w1ref run (two
relaxed monomers at level 1). Without it the row is left out and the totals
are interaction energies at the cluster geometry rather than binding energies.

Under `ssfc` the one-body term here is each monomer in its own basis; the
paper's corrected one is in the cluster's basis. The rows of order two and
above are the same quantity as the paper's.
"""
import argparse
import json

KCAL = 627.5094740631


def columns(obj):
    """HF, correlation and its parts of a `levels` entry or the top level."""
    parts = obj.get("correlation_parts", {})
    c = {"HF": obj.get("scf_energy", obj["total_energy"])}
    if "correlation_energy" in obj:
        c["corr"] = obj["correlation_energy"]
    if "cc_doubles" in parts:
        c["CCSD"] = parts.get("cc_singles", 0.0) + parts["cc_doubles"]
    if "cc_triples" in parts:
        c["(T)"] = parts["cc_triples"]
    c["total"] = obj["total_energy"]
    return c


def mbe(path):
    doc = json.load(open(path))
    return doc["mbe"] if "mbe" in doc else next(iter(doc.values()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("output")
    ap.add_argument("--ref", help="output of the matching w1ref deck")
    ap.add_argument("--n", type=int, help="number of waters (default: the highest level present)")
    args = ap.parse_args()

    m = mbe(args.output)
    levels = sorted(m["levels"], key=lambda l: l["frag_level"])
    n = args.n or levels[-1]["frag_level"]
    rows = {l["frag_level"]: columns(l) for l in levels}
    names = list(rows[1])

    if args.ref:
        one = columns(mbe(args.ref)["levels"][0])
        rows[1] = {k: rows[1][k] - n * one[k] / 2.0 for k in names}
    else:
        del rows[1]

    total = {k: sum(r[k] for r in rows.values()) for k in names}
    print(f"{args.output}   counterpoise: {m.get('counterpoise', 'none')}   kcal/mol (% of column)")
    print("body " + "".join(f"{k:>20}" for k in names))
    for lvl in sorted(rows):
        cells = "".join(f"{rows[lvl][k] * KCAL:11.3f} ({100 * rows[lvl][k] / total[k]:6.1f})" for k in names)
        print(f"{lvl}B   " + cells)
    print("sum  " + "".join(f"{total[k] * KCAL:11.3f}{'':9}" for k in names))
    if not args.ref:
        print("no --ref: the one-body (monomer relaxation) row is absent")


if __name__ == "__main__":
    main()
