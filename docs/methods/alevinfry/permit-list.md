---
title: permit list
parent: alevin-fry
grand_parent: Methods
nav_order: 1
method_ids: [alevinfry-permit-list]
status: implemented
---

# alevin-fry generate-permit-list
{: .no_toc }

`alevin-fry generate-permit-list`: the barcodes kept as cells, before quantification. The rule
gives a threshold on each barcode's count, and every barcode at or above it is kept, so ties never
matter. The knee is a distance knee like [UMI-tools'](../umitools/distance.md), but repeated on
five times (not three times) as many barcodes, and the last of equally distant points wins.

Why the rule lands where it does, drawn on real data:
[alevin-fry's knee](../../explained/knee/alevinfry.md), in *Explained*.

1. TOC
{:toc}

## At a glance

| Id | `alevinfry-permit-list` |
|---|---|
| Tool | alevin-fry 0.1.0 and later (`generate-permit-list -k / -e / -f / -u`; simpleaf `--knee`, `--expect-cells`, `--forced-cells`, `--unfiltered-pl`) |
| Category | knee |
| Counts | reads per barcode, before barcode correction, in alevin-fry; here UMIs, or reads with `--unit reads` |
| Needs expected cells | no (optional) |
| Deterministic | yes |
| Original licence | BSD-3-Clause |
| Validated against | alevin-fry 0.1.0, 0.11.2, 0.17.1 and 0.19.0's own functions (exact) |

## How it works

The barcodes with a count are sorted by decreasing count: $$f_0 \ge f_1 \ge \dots \ge f_{n-1}$$
(0-based). Each rule gives a threshold $$t$$; the cells are the barcodes with count
$$\ge t$$.

**Knee** (`--knee-distance`, the default here). The cumulative counts are
$$v_i = f_0 + \dots + f_i$$.

1. The farthest point of the first $$m$$ cumulative counts: with $$x_i = i / m$$ and
   $$y_i = v_i / v_{m-1}$$, the line runs from $$(0, y_0)$$ to $$(1, 1)$$, and the distance of
   $$(x_i, y_i)$$ from it is
   $$|(1 - y_0)\,x_i - y_i + y_0| \,/\, \sqrt{(1 - y_0)^2 + 1}$$, in alevin-fry's order of
   operations. The farthest point is the last one at the largest distance.
2. Start with $$k$$ = the farthest point of all $$n$$ counts. While $$k$$ changes (at most 100
   times), $$k$$ becomes the farthest point of the first $$\min(n - 1, 5k)$$ counts. alevin-fry
   stops with an error when $$k$$ is 0 or fewer than two counts are left.
3. $$t = f_k$$.

**Expected cells** (`--expect-cells` $$E$$): $$r = \min(n - 1, \operatorname{round}(0.99\,E))$$
and $$t = \max(1, \operatorname{round}(f_r / 10))$$, rounding halves away from zero (Rust's
`round`). This is a tenth of the count at rank $$0.99E$$, not of the 99th percentile of the top
$$E$$ barcodes (Cell Ranger 2's rule), so it keeps many more barcodes.

**Forced cells** (`--force-cells` $$F$$): $$t = f_{\min(F, n) - 1}$$. Ties with the $$F$$-th
barcode are kept, so more than $$F$$ barcodes can pass.

**Unfiltered permit list** (`--unfiltered-pl` with `--min-reads`, 10 by default): $$t$$ is the
floor itself, applied to the barcodes on the given list. Here every barcode of the matrix is on
the list.

**Multi-barcode input** (10x Flex: a sample barcode and a cell barcode), 0.12.0 – 0.17.1 only,
`--flex`. Each sample's permit list uses rules of its own:

- knee: $$t = f_{k-1}$$, one rank above the single-barcode knee;
- expected cells: $$s = \lfloor f_0 / E \rfloor$$, $$p$$ is the first rank with $$f_p < s$$
  (else $$n$$), $$c = \min(\lfloor 10p \rfloor, n)$$, and $$t = f_{c-1}$$;
- forced cells: $$t = f_{\min(F, n) - 1}$$.

From 0.18.0 multi-barcode input uses the single-barcode rules.

## Versions

| alevin-fry | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.1.0 – 0.4.1 | `v0.1` | `0.1.0` | the knee, expected cells, forced cells; no multi-barcode input |
| 0.4.2 – 0.11.2 | `v0.4.2` | `0.4.2` | adds the unfiltered permit list |
| 0.12.0 – 0.17.1 | `v0.12` | `0.12.0` | multi-barcode input (`--flex`) has its own rules |
| 0.18.0 – 0.19.0 | `v0.18` | `0.18.0` | one rule set for both |


## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expect-cells` | none | `--expect-cells`, `-e` | the rule from expected cells |
| `--force-cells` | none | `--force-cells`, `-f` | the rule from forced cells |
| `--unfiltered` | off | `--unfiltered-pl`, `-u` | keep every barcode with at least `--min-reads` |
| `--min-reads` | 10 | `--min-reads` | the floor of the unfiltered permit list |
| `--flex` | off | multi-barcode RAD input | the 0.12 – 0.17 rules for 10x Flex |
| `--unit` | `umis` | reads per barcode | count UMIs (the matrix) or reads (`--reads`) |

With none of the first three, the knee is used (`--knee-distance`, `-k`).

## Run it

```bash
sc-droplets call alevinfry-permit-list raw_feature_bc_matrix.h5 --json
sc-droplets call alevinfry-permit-list raw_feature_bc_matrix.h5 --expect-cells 5000 --barcodes cells.tsv
sc-droplets call alevinfry-permit-list raw_feature_bc_matrix.h5 --tool-version 0.17.1 --flex --force-cells 3000
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "alevinfry-permit-list")
result.thresholds["umis"], result.thresholds["knee_rank"]
```

The thresholds hold the cut (in the unit counted) and, for the knee, its rank.
`sc-droplets plot diagnostics` draws the rank curve with the cut, and the cumulative curve of the
last pass with its line and each pass's knee.

## Differences from the original

- alevin-fry counts the reads of each barcode as sequenced, before barcode correction, so a
  barcode one error away from a cell's adds nothing to it. A matrix holds corrected barcodes, and
  sc-droplets counts their UMIs by default; reads need `--reads`.
- `--explicit-pl` (`-b`, a list given) is not offered: it is no rule.

## Validation

`oracles/alevinfry/run.sh` builds a small Rust program per release from alevin-fry's own source
(fetched at the tag): the knee functions as they are, and each rule's lines as they are in
`cellfilter.rs`. It runs them on the per-barcode UMI totals of 10x's pbmc_1k_v3 and pbmc4k
raw matrices, for 0.1.0, 0.11.2, 0.17.1 (single and multi-barcode rules) and 0.19.0: the
knee, expected cells and forced cells, and on 0.19.0 the unfiltered floor. Identical cells, and
identical thresholds where the rule computes one.

## References

- He, Zakeri, Sarkar, Soneson, Srivastava & Patro 2022, Alevin-fry unlocks rapid, accurate and
  memory-frugal quantification of single-cell RNA-seq data, *Nature Methods* 19:316–322,
  [doi:10.1038/s41592-022-01408-3](https://doi.org/10.1038/s41592-022-01408-3).
- Source: [`knee_finding.rs`](https://github.com/COMBINE-lab/alevin-fry/blob/v0.19.0/src/knee_finding.rs)
  and [`cellfilter.rs`](https://github.com/COMBINE-lab/alevin-fry/blob/v0.19.0/src/cellfilter.rs)
  (`select_retained_barcodes`, 0.19.0);
  [`cellfilter.rs` at 0.17.1](https://github.com/COMBINE-lab/alevin-fry/blob/v0.17.1/src/cellfilter.rs)
  (the multi-barcode rules) and
  [`libradicl/src/cellfilter.rs` at 0.1.0](https://github.com/COMBINE-lab/alevin-fry/blob/v0.1.0/libradicl/src/cellfilter.rs).
