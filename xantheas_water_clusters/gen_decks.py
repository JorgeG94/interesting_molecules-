#!/usr/bin/env python3
"""Decks for the many-body expansion of CCSD(T) binding energies of water clusters.

The clusters, geometries and basis sets are those of Nguyen and Xantheas,
J. Chem. Theory Comput. 2026, 22, 9888 (MP2/aug-cc-pVXZ geometries, each run in
the basis it was optimised with). Their Tables 1 and 2 split the binding energy
into Hartree-Fock and MP2 correlation and expand each over bodies; these decks
do the same with RI-CCSD(T), whose output splits the correlation further into
singles, doubles and triples at every body order, and with RI-MP2, which
reproduces the paper's own numbers as a check of the setup.

    python3 gen_decks.py              # density-fitted SCF, the default
    python3 gen_decks.py --exact-scf  # four-centre SCF, as the paper's HF is

Per cluster, basis and method there are two decks, both an MBE carried to the
full order N so the sum is the supermolecular energy:

    <cluster>__<method>__none.json   every term in its own basis
    <cluster>__<method>__ssfc.json   every term of order >= 2 in the cluster's basis

(H2O)20 is taken to level 4 (the paper stops at 2), with a separate
unfragmented deck for the total.

The one-body term needs the relaxed monomer. `w1__opt.json` optimises it at
RI-MP2 (needs a build with DL-FIND); save the result as w1_<aVXZ>.xyz beside this script, or
put the SI's monomer there, and rerun this script: it then writes
`w1ref__<method>.json`, two copies of that monomer run as a level-1 expansion
so its output has the same SCF/correlation split as the clusters'.
"""
import argparse
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEOMS = HERE
MQC_BASIS_DIR = HERE.parent.parent / "metalquicha" / "basis_sets"

BASES = {"aVDZ": "aug-cc-pvdz", "aVTZ": "aug-cc-pvtz", "aVQZ": "aug-cc-pvqz", "aV5Z": "aug-cc-pv5z"}
METHODS = {"ri-mp2": "ri-mp2", "ri-ccsdt": "ri-ccsd(t)"}
# Spherical functions per water, orbital and RIFIT, for the size table only.
NBF = {"aVDZ": (41, 118), "aVTZ": (92, 198), "aVQZ": (172, 328), "aV5Z": (287, 496)}
FULL_ORDER_MAX = 6
LARGE_CLUSTER_LEVEL = 4


def dist(a, b):
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def comb(n, k):
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))


def read_xyz(path):
    lines = path.read_text().splitlines()
    n = int(lines[0])
    atoms = []
    for line in lines[2:2 + n]:
        s, x, y, z = line.split()[:4]
        atoms.append((s, float(x), float(y), float(z)))
    return atoms


def water_fragments(atoms):
    """One fragment per oxygen with its two nearest hydrogens, 0-based."""
    oxygens = [i for i, a in enumerate(atoms) if a[0] == "O"]
    frags = {o: [o] for o in oxygens}
    for i, a in enumerate(atoms):
        if a[0] != "H":
            continue
        d, o = min((dist(a[1:], atoms[o][1:]), o) for o in oxygens)
        if d > 1.1:
            raise SystemExit(f"hydrogen {i} is {d:.2f} A from any oxygen")
        frags[o].append(i)
    if any(len(f) != 3 for f in frags.values()):
        raise SystemExit("a water does not have two hydrogens")
    return [frags[o] for o in oxygens]


def deck(name, xyz, basis, method, args, fragments=None, level=None, cp="none", backend=None):
    scf = {"maxiter": 200, "tolerance": 1e-10, "guess": "sad", "density_fitting": not args.exact_scf}
    d = {
        "schema": {"name": name, "version": "1.0"},
        "molecules": [{"xyz": xyz, "molecular_charge": 0, "molecular_multiplicity": 1}],
        "model": {"method": method, "basis": basis, "aux_basis": basis + "-rifit"},
        "backend": backend or args.backend,
        "keywords": {"scf": scf, "correlation": {"freeze_core": True}},
        "system": {"logger": {"level": "Verbose"}},
        "driver": "Energy",
    }
    if "cc" in method:
        d["keywords"]["cc"] = {"tolerance": 1e-9, "maxiter": 100}
    if fragments is not None:
        mol = d["molecules"][0]
        mol["fragments"] = fragments
        mol["fragment_charges"] = [0] * len(fragments)
        mol["fragment_multiplicities"] = [1] * len(fragments)
        frag = {"method": "mbe", "level": level, "embedding": "none"}
        if cp != "none":
            frag["counterpoise"] = cp
        d["keywords"]["fragmentation"] = frag
    return d


def write(path, d):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(d, indent=2) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--exact-scf", action="store_true", help="four-centre SCF instead of a fitted one")
    ap.add_argument("--backend", default="terco")
    ap.add_argument("--out", default="decks")
    ap.add_argument("--basis-dir", default=str(MQC_BASIS_DIR),
                    help="metalquicha's basis_sets/, to skip a basis it does not have")
    args = ap.parse_args()
    out = HERE / args.out
    basis_dir = Path(args.basis_dir)

    rows, skipped = [], set()
    for xyz in sorted(GEOMS.glob("w*_aV?Z.xyz")):
        m = re.fullmatch(r"(w(\d+)(?:_[a-z]+)?)_(aV[DTQ5]Z)", xyz.stem)
        cluster, n, tag = m.group(1), int(m.group(2)), m.group(3)
        basis = BASES[tag]
        if not (basis_dir / f"{basis}.json").exists():
            skipped.add(basis)
            continue
        rel = f"../../{xyz.name}"
        if n == 1:
            # The relaxed monomer, twice, 100 A apart, at level 1: 2 E(1).
            atoms = read_xyz(xyz)
            two = GEOMS / f"w1x2_{tag}.xyz"
            body = [f"{s} {x:16.9f} {y:16.9f} {z:16.9f}" for s, x, y, z in atoms]
            body += [f"{s} {x + 100.0:16.9f} {y:16.9f} {z:16.9f}" for s, x, y, z in atoms]
            two.write_text(f"6\ntwo copies of {xyz.name}, 100 A apart; Angstrom\n" + "\n".join(body) + "\n")
            for mtag, method in METHODS.items():
                write(out / tag / f"w1ref__{mtag}.json",
                      deck(f"w1ref_{tag}_{mtag}", f"../../{two.name}", basis, method, args,
                           fragments=[[0, 1, 2], [3, 4, 5]], level=1))
            continue
        frags = water_fragments(read_xyz(xyz))
        assert len(frags) == n
        level = n if n <= FULL_ORDER_MAX else LARGE_CLUSTER_LEVEL
        for mtag, method in METHODS.items():
            for cp in ("none", "ssfc"):
                write(out / tag / f"{cluster}__{mtag}__{cp}.json",
                      deck(f"{cluster}_{tag}_{mtag}_{cp}", rel, basis, method, args,
                           fragments=frags, level=level, cp=cp))
            if level < n:
                write(out / tag / f"{cluster}__{mtag}__full.json",
                      deck(f"{cluster}_{tag}_{mtag}_full", rel, basis, method, args))
        nterms = sum(comb(n, k) for k in range(1, level + 1))
        ao, aux = NBF[tag]
        rows.append((tag, cluster, n, level, nterms, nterms + n, n * ao, n * aux))

    for tag, basis in BASES.items():
        if not (basis_dir / f"{basis}.json").exists():
            continue
        # The optimiser is on the CPU backend's RI-MP2 gradient.
        d = deck(f"w1_opt_{tag}", "../../w1_start.xyz", basis, "ri-mp2", args, backend="cpu")
        d["driver"] = "Optimize"
        write(out / tag / "w1__opt.json", d)
        if not (GEOMS / f"w1_{tag}.xyz").exists():
            print(f"note: w1_{tag}.xyz is missing, so no w1ref decks for {tag} "
                  f"(run {args.out}/{tag}/w1__opt.json, or take the monomer from the SI)")

    print(f"\n{'basis':5} {'cluster':10} {'N':>2} {'level':>5} {'terms':>6} {'ssfc':>5} "
          f"{'nbf(full)':>9} {'naux(full)':>10}")
    for r in sorted(rows, key=lambda r: ("DTQ5".index(r[0][2]), r[2], r[1])):
        print(f"{r[0]:5} {r[1]:10} {r[2]:2d} {r[3]:5d} {r[4]:6d} {r[5]:5d} {r[6]:9d} {r[7]:10d}")
    for basis in sorted(skipped):
        print(f"skipped: no {basis}.json in {basis_dir}, so nothing written for it")
    print(f"\n{sum(1 for _ in out.rglob('*.json'))} decks under {out}")


if __name__ == "__main__":
    main()
