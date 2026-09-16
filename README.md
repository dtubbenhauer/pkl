# Code and Erratum for *Torsion proliferation*

This repository contains data and plotting material for the computations accompanying the paper

*Torsion proliferation*

by Joseph Baine and Daniel Tubbenhauer.

The computations compare the p-canonical basis element `pC(w)` with the ordinary Kazhdan--Lusztig basis element `C(w)` in finite Coxeter types.  An element is counted as **changed** if `pC(w) != C(w)`.  Equivalently, in the normalization used in the scripts, at least one local stalk polynomial below `w` differs from the characteristic-zero Kazhdan--Lusztig stalk polynomial.

## Contact

If you find any errors in the paper **please email me**:

[dtubbenhauer@gmail.com](mailto:dtubbenhauer@gmail.com?subject=[GitHub]%web-reps)

Same goes for any errors related to this page.

## Shoutout

All the computations were done using Gibson and Williamson's IHecke and ASLoc, see [IHecke Repository](https://github.com/joelgibson/IHecke) and [ASLoc Repository](https://github.com/joelgibson/ASLoc)

## What is in this repository?

```text
data/raw/          raw Magma output, grouped by summary type
data/processed/    CSV files parsed from the raw output
docs/              GitHub Pages page with plots
docs/figures/      PNG and PDF versions of all figures
scripts/           Magma summary scripts used to generate the raw output
```

The GitHub Pages entry point is `docs/index.html`.

## Main figures

The most useful paper/GitHub figures are:

- `docs/figures/01_changed_percentage_p2.png`: overview of all available finite-type computations at `p = 2`.
- `docs/figures/02_classical_rank_trends_p2.png`: rank trends in the classical families at `p = 2`.
- `docs/figures/03_prime_comparison_p2_vs_p3.png`: direct comparison of `p = 2` and `p = 3` on the same finite types.
- `docs/figures/04_correction_severity_stacked_p2.png`: mild/moderate/wild decomposition among changed elements.
- `docs/figures/05_genuinely_graded_fraction_p2.png`: fraction of changed elements whose correction is not concentrated in degree zero.
- `docs/figures/06_length_heatmap_p2.png`: changed percentage by Coxeter length.
- `docs/figures/09_full_vs_parabolic_support_p2.png`: full-support versus proper-parabolic support comparison.

## Processed data

The processed CSV files are:

- `summary_percentage.csv`: total elements, changed elements, percentage changed, and first changed length.
- `length_distribution.csv`: changed counts and percentages by Coxeter length.
- `summary_polynomial.csv`: global correction statistics in the ordinary KL basis.
- `polynomial_length_distribution.csv`: polynomial correction statistics by length.
- `full_support_summary.csv`: full-support versus proper-parabolic support.
- `support_size_profile.csv`: changed percentage by support size.
- `descent_profile.csv`: changed percentage by number of left and right descents.
- `bruhat_depth_distribution.csv`: distribution of correction summands by Bruhat depth.

## Correction conventions

For the polynomial summaries, the correction is measured in the ordinary KL basis:

```text
Delta_w = pC(w) - C(w) = sum_x f_{x,w}(v) C(x).
```

The summary scripts also record several crude size measures:

- `signed_mass_at_one`: sum of `f_{x,w}(1)`.
- `abs_mass_at_one`: sum of the absolute values of all monomial coefficients in all `f_{x,w}(v)`.
- `monomial_complexity`: total number of nonzero monomials in the correction.
- `genuinely_graded`: the correction has some nonzero-degree term.

For the mild/moderate/wild plots, the convention is:

```text
mild      signed_mass_at_one = 1
moderate  2 <= signed_mass_at_one <= 5
wild      signed_mass_at_one > 5
```

## Reproducing the summaries

The Magma scripts in `scripts/` expect the ASLoc setup and saved p-canonical bases.  Typical usage is:

```bash
magma -b type:=B5 prime:=2 saveDir:=saves summary-percentage.m > b5-2.txt
magma -b type:=B5 prime:=2 saveDir:=saves summary-polynomial.m > b5-2-polynomial.txt
magma -b type:=B5 prime:=2 saveDir:=saves summary-structure-KL.m > b5-2-structure-KL.txt
```

The raw outputs in `data/raw/` are enough to regenerate all plots.

## Current computed range

The uploaded data cover:

- `p = 2`: `A7`, `B2`--`B6`, `C2`--`C6`, `D4`--`D6`, `E6`, `F4`, `G2`.
- `p = 3`: `B6`, `C6`, `D6`, `E6`, `F4`, `G2`.

Some low-rank cases at `p = 3` are absent because the available run set focused on the first nontrivial/high-value comparisons.

## Notes for the paper

The computations are not used in the proof.  They are intended to show what already happens in small rank and to illustrate that the corrections are not merely present but can become large, spread out in Bruhat order, and genuinely graded.

## Erratum

Empty so far.
