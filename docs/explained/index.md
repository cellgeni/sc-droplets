---
title: Explained
nav_order: 4
has_children: true
---

# Explained
{: .no_toc }

The [Methods](../methods/index.md) pages say exactly what each tool computes: its rule, release
by release, and how far this package reproduces it. These pages say why the rule works, and
where it stops working. Each takes one idea at a time, draws it on real data, derives what
follows from it, and shows where it breaks. A detail that matters for reproducing a tool but not
for understanding it is named and linked, not dwelt on.

1. TOC
{:toc}

## How to read them

Start with [The barcode-rank curve](rank-curve.md). Every family looks at that curve, or at
something made from it. Then read a family's page, which builds the idea its members share. A
tool's page is that idea plus the choices one tool made: which curve, which barcodes, which of
two equal maxima.

The pages assume you are at ease with derivatives, logarithms and a little probability. Formulas
come after the picture they describe, never instead of it.

## The families

| Family | The idea | Pages |
|---|---|---|
| [Knee geometry](knee/index.md) | the point where the rank curve bends | the family, and nine tools |
| Rank and threshold | a count read off the top of the curve: a robust maximum divided by ten, the top N | to come |
| Density | the valley between the humps of log counts | to come |
| Mixture | a fitted mixture whose components are cells and background | to come |
| Test against ambient | is a droplet's content unlike the ambient RNA? | to come |
| Classifier | a model trained on droplets' features | to come |

The families are those of the [Overview](../overview/index.md#families).
[Same idea, different tools](../overview/families.md) runs each family's members side by side on
the same matrices.

## Notation

| Symbol | Meaning |
|---|---|
| $$t_b$$ | the total of barcode $$b$$: UMIs, unless a page says reads |
| $$c_0 \ge c_1 \ge \dots \ge c_{n-1}$$ | the totals above 0, sorted; the barcode at index $$i$$ has rank $$r = i + 1$$ |
| $$v_i = c_0 + \dots + c_i$$ | the cumulative count of the top $$i + 1$$ barcodes |
| $$\log$$ | base 10 in the figures; the base never changes where a knee falls |

## The data behind the figures

Every figure and number comes from two of 10x's raw matrices:

- **PBMC 1k** (whole cells, 3′ v3): 234,600 barcodes with a count, of which 10x called 1,222
  cells.
- **Mouse spinal cord nuclei** (3′ v4, Cell Ranger 10.0): 1,168,887 barcodes with a count, of
  which 10x called 15,257.

The first has a knee anyone can see. The second does not, and that is where the rules part.

`docs/_scripts/explained.py` draws the figures with the methods' own code. The identities the
pages derive are checked against that code twice: by `tests/test_explained.py` on simulated
data, and by the script on both matrices. The interactive figures compute from those identities
on the same two matrices. They illustrate; the numbers the pages quote come from the package.
