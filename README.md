# ancientDNASymmetry

Substitution symmetry analysis for ancient DNA (aDNA) and ancient environmental DNA
(aeDNA), from:

> Lemmon-Kishi M, Pipes L, De Sanctis B, Nielsen R. **Molecular Clock Dating of Ancient
> Environmental DNA Reveals Damage Beyond Deamination.** *bioRxiv* (2026).
> [doi:10.64898/2026.07.03.735781](https://doi.org/10.64898/2026.07.03.735781)

This is one of two code repositories for the paper. The molecular dating method
(ratePlacer), gfl, the simulations and helper scripts are in
[lemmonquiche/ratePlacer](https://github.com/lemmonquiche/ratePlacer). This repository
contains the damage and symmetry analysis.

> Work in progress. Will be updating as I transition to a postdoc.

## Motivation

Post-mortem damage in aDNA shows up as mismatches between reads and the reference. The
best-known signature is cytosine deamination: in double-stranded libraries, elevated
C → T at 5′ termini and complementary G → A at 3′ termini. Most damage-aware tools model
only deamination, so any other damage is left in the data. If it isn't accounted for, it
is mistaken for evolutionary divergence, lengthening ancient branches and biasing
molecular age estimates.

This repository looks for damage that isn't deamination, in two complementary ways.

**Substitution symmetry**. Under any
time-reversible substitution model, detailed balance (πᵢQᵢⱼ = πⱼQⱼᵢ) means i → j and
j → i substitutions are expected in equal numbers. The symmetry ratio Rᵢⱼ = Nᵢⱼ / Nⱼᵢ
should therefore be about 1, and an excess of one direction points to damage. This
needs no model fit, so it can be compared across many taxa and between ancient and
modern reads.

**Observed-to-expected ratios**.
Without damage or sequencing error, mismatches between reads and the reference reflect
only evolution. Their expected counts follow from a substitution model fitted to the
taxon's reference phylogeny: nπᵢQᵢⱼt for i → j. The elapsed time *t* is unknown, but it
scales all substitution types equally, so the observed-to-expected ratio
Oᵢⱼ = Nᵢⱼ / nπᵢQᵢⱼ should be roughly uniform across substitution types. Any type that
stands out points to post-mortem misincorporation. Unlike symmetry ratios, this puts all
12 substitution types on one scale, so they can be ranked. It can also be resolved by
read position, which shows where along the read damage is concentrated. Because it needs
a per-taxon model fit, it is run one taxon at a time. Expected counts are computed under
GTR+Γ, and recomputed under a non-reversible UNREST model to check that the
reversibility assumption doesn't drive the result from the symmetry analysis.

Together, these analyses revealed, beyond deamination, elevated G → T at 5′ termini and
C → A at 3′ termini, consistent with putative oxidative damage of guanine (8-oxoguanine). Both
deamination and oxidative signatures are strongest at the termini but persist across
the read interior.

## Repository layout

The two analyses are independent pipelines. Each has its own scripts, example data and
documentation.

| Directory | Analysis | Paper |
| --- | --- | --- |
| [`single_taxon/`](single_taxon/) | Aggregate and positional substitution counts and observed-to-expected ratios for one BAM (one taxon), under GTR or UNREST | Sections 4.2.2, 4.5.2; Table 1, Fig. 2, Supp. Figs S3–S5 |
| [`aggregate/`](aggregate/) | Ancient vs modern symmetry ratios across taxa, χ² / Cramér's V, label permutation, Watson-Crick complement pairs | Sections 4.2.1, 4.2.3; Tables 2–3, Fig. 3, Supp. Figs S6–S7 |

## Requirements

Python and R dependencies are listed in [`environment.yml`](environment.yml):

```bash
conda env create -f environment.yml
```

## Citation

If you use this code, please cite:

```bibtex
@article{lemmonkishi2026ratePlacer,
  title   = {Molecular Clock Dating of Ancient Environmental DNA Reveals Damage Beyond Deamination},
  author  = {Lemmon-Kishi, Maya and Pipes, Lenore and De Sanctis, Bianca and Nielsen, Rasmus},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.07.03.735781}
}
```

## Contact

Maya Lemmon-Kishi — maya_lemmon-kishi@berkeley.edu
