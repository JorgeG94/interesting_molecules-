# more_water_clusters

Geometries of small water clusters from the Supporting Information of
L. H. Nguyen and S. S. Xantheas, "Many-Body Expansion of the MP2 Correlation
Contribution to the Binding Energies of Water Clusters: A Localized Excitation
Analysis", J. Chem. Theory Comput. (SI file `ct6c01287_si_001`, pages S7 to S22).

Each file is a standard xyz file named `w<n>_<isomer>_<basis>.xyz`. The basis
is the one the geometry was optimized with at the MP2 level: `aVDZ`, `aVTZ`,
`aVQZ` and `aV5Z` are aug-cc-pVDZ to aug-cc-pV5Z. The comment line holds the HF
and MP2 correlation energies the SI lists for that structure (Eh).

| Cluster | Isomers | Basis sets |
|---|---|---|
| (H2O)2 | | aVDZ, aVTZ, aVQZ, aV5Z |
| (H2O)3 | ring | aVDZ, aVTZ, aVQZ, aV5Z |
| (H2O)3 | aa, dd | aVDZ, aVTZ, aVQZ |
| (H2O)4 | ring | aVDZ, aVTZ, aVQZ |
| (H2O)5 | ring | aVDZ, aVTZ, aVQZ |
| (H2O)6 | ring, book, cage, prism | aVDZ, aVTZ, aVQZ |
| (H2O)20 | dodecahedron | aVDZ |

Notes:

- The SI heading says the coordinates are in a.u., but the O-H distances are
  about 0.96, so they are in Angstrom, and the files are written as such.
- The monomer geometries listed in the SI are not included.
- The coordinates were extracted from the PDF text and checked only for atom
  counts and that every O has two H atoms at O-H distances of 0.9 to 1.05.
