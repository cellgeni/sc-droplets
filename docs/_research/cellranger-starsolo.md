# Cell Ranger, Cell Ranger ATAC/ARC and STARsolo: cell calling by version

Source code is the authority wherever it exists. Tags: **[SRC]** read at the tag named, **[DOC]**
10x documentation or release notes only, **UNVERIFIED** inferred.

## 0. What source exists

| Product | Public source | Where |
|---|---|---|
| Cell Ranger 2.0.2, 2.1.0, 2.1.1, 3.0.2, plus development commits from 2.1 to 3.0.2 | Yes (upstream history was replaced; forks keep it) | `github.com/lfperales/cellranger` (tags `2.0.2`, `2.1.0`, `3.0.2`) |
| Space Ranger 1.0 (≈ Cell Ranger 3.1 code) | Yes (third-party mirror) | `github.com/sbooeshaghi/spaceranger` |
| Cell Ranger 1.x, 2.2, 3.1, 4.0, 5.0, 6.0, 6.1, 7.0, 7.1 | No. Only in binary tarballs under the 10x EULA; this project does not use them | — |
| Cell Ranger 7.2.0, 8.0.0, 8.0.1, 9.0.0, 9.0.1, 10.0.0, 10.1.0 | Yes | `github.com/10XGenomics/cellranger`, tags `cellranger-X.Y.Z` |
| cellranger-atac 1.2 | Yes (single commit) | `github.com/10XGenomics/cellranger-atac` |
| cellranger-atac 1.0, 1.1, 2.0–2.2 | No | — |
| cellranger-arc | Only the CLI wrapper (`cellranger/lib/rust/cr_wrap/src/arc/`) | — |
| STAR 2.7.0a – 2.7.11b and master | Yes | `github.com/alexdobin/STAR` |

Key files:

- Cell Ranger 2.x: `lib/python/cellranger/stats.py`, `constants.py`,
  `mro/stages/counter/filter_barcodes/__init__.py`.
- Cell Ranger 3.0.2: adds `lib/python/cellranger/cell_calling.py`, `sgt.py`,
  `analysis/diffexp.py` (`adjust_pvalue_bh`), `feature/antibody/analysis.py`.
- Cell Ranger 7.2+: `lib/python/cellranger/cell_calling.py`, `cell_calling_helpers.py`,
  `stats.py` (renamed `cell_calling_stats.py` in 10.1), `sgt.py`, `chemistry.py` with
  `chemistry_defs.json`, `feature/antibody/analysis.py`,
  `mro/rna/stages/counter/filter_barcodes/__init__.py`.

URL pattern: `https://github.com/10XGenomics/cellranger/blob/cellranger-10.1.0/lib/python/cellranger/cell_calling.py`;
old versions `https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/cell_calling.py`.

## 1. Licences

| Code | Licence in the tree | For this project |
|---|---|---|
| Cell Ranger 2.0.2, 2.1.0, 2.1.1 | AGPLv2 | Not copied. Implemented from our own spec. |
| Cell Ranger commits 2018-08-03 ("Change license to MIT") to 2020-07-30, including tag 3.0.2 | MIT (Louvain excepted) | 3.0.2 code may be reused with its notice. |
| Cell Ranger after 2020-07-31, all 7.2–10.1 tags; cellranger-atac | 10x licence: internal use, only with 10x products, no redistribution, derivatives bound by the same terms | **Never copied.** Implemented only from the spec on our docs page. |
| STAR | MIT | May be reused with notice. `source/SimpleGoodTuring/sgt.h` is © David Elworthy 2004, permissive (keep the notice). |

## 2. Version timeline

| Rule | Versions | Step 1 | Step 2 / extras | Status |
|---|---|---|---|---|
| `cellranger_v1_ordmag` | CR 1.0–1.3 | 99th pct of top N, ÷10 | — | [DOC] Zheng 2017 methods |
| `cellranger_v2_ordmag` | CR 2.0–2.2 | OrdMag, 100 bootstraps, N=3000 default | per genome | [SRC] 2.0.2, 2.1.0 |
| `cellranger_v3_emptydrops` | CR 3.0.0–3.0.2 (3.1 ≈ same) | `v2_ordmag` | ambient ranks [45000, 90000), SGT, top 20000 candidates, min UMI max(500, 1 % median), 10000 sims, BH ≤ 0.01 | [SRC] 3.0.2; Space Ranger 1.0 ≈ 3.1 |
| (gap) | CR 4.0–6.1 | 4.0 targeted gradient; FB-only = OrdMag on antibody. 5.0 protein aggregate filter. 6.0 genome-assignment bug fixed. 6.1 per-species EmptyDrops, "changed certain parameters" | — | [DOC]; code UNVERIFIED |
| `cellranger_v7_ordmag_auto` + `cellranger_v7_emptydrops` | 7.0 (cap 2^18), 7.1/7.2 (cap = ambient lower bound) | auto-estimated expect-cells | no candidate cap; Flex/LT ranges; HOG filter (Flex) | [SRC] 7.2.0; 7.0/7.1 [DOC] |
| `cellranger_v8_emptydrops` | 8.0.0, 8.0.1 | as v7 | min UMI = max(500, 1 + max ambient UMI); v4/5′v3 range 80k–160k and FDR 0.001; 8.0.1 HT 80k–160k; hidden global-min-UMI and mito filters | [SRC] |
| `cellranger_v9_emptydrops` | 9.0.0, 9.0.1 | as v7 | 100000 sims, new simulator (PCG64 seed 42); FDR 0.001 for 3′v3/3′v4/5′/HT; Flex 90000×nBC; OCM 40000×nBC | [SRC] |
| `cellranger_v10_emptydrops` | 10.0.0, 10.1.0 | as v7 | Flex v2 160000×nBC; feature-ID order in simulation; HOG uses NB; 10.1 drops the LT special case | [SRC] |
| `cellranger_targeted_gradient` | 4.0–9.x (auto for hybrid-capture targeted); 10.1 hidden `gradient` mode | — | — | [SRC] 7.2+ |
| `cellranger_atac_v1_nbmix` | ATAC 1.0–1.2 | 2-component NB mixture EM, odds 100000 | — | [SRC] 1.2 |
| `cellranger_atac_v2` | 2.0–2.2 | same model per docs; 2.1 less stringent low-targeting; 2.2 drops the 20k cap | — | [DOC] |
| `cellranger_arc_joint_v1/v2` | ARC 1.x/2.x | 2-D OrdMag then k-means (k=2) | — | [DOC] |
| `starsolo_*` | STAR 2.7.3a+ | see §6 | — | [SRC] |

## 3. Cell Ranger, rule by rule

### 3.0 Shared primitives

- **UMI counts per barcode**: column sums of the feature × barcode UMI matrix, restricted to the
  chosen feature types (Gene Expression by default) and, where stated, to one genome's features.
- **BH adjustment** (`adjust_pvalue_bh`, identical 3.0.2 → 10.1) [SRC]:
  ```
  desc = argsort(p)[::-1]; scale = n / arange(n, 0, -1)
  q = minimum(1, minimum.accumulate(scale * p[desc])); return q[argsort(desc)]
  ```
- **SGT ambient profile** (`sgt.py`, identical 3.0.2 → 10.1; a port of Gale & Sampson's S code) [SRC]:
  ```
  sgt_proportions(f):              # f = positive ambient counts, one per feature
    nr = bincount(f); r = flatnonzero(nr)
    if len(r) < 10: raise SGTError -> EmptyDrops step aborted (no extra cells)
    rstar, p0 = simple_good_turing(r, nr[r])
    pstar_i = (1-p0) * rstar[f_i] / sum(nr[r]*rstar)
  simple_good_turing(r, n):
    N = sum(r*n); d = [1] + diff(r); dr = [0.5*(d[1:]+d[:-1]), d[-1]]; Z = n/dr
    b, a = linregress(log r, log Z); if b > -1: raise SGTError
    lgt_rel = (1+1/r)**(1+b)
    cons[i] = r[i+1] == r[i]+1 (last = False)
    tur_rel[i] = (r_i+1)/r_i * n_{i+1}/n_i if cons[i] else 0
    tursd[i] = (i+2)/n_i * sqrt(n_{i+1}*(1+n_{i+1}/n_i)) if cons[i] else 1   # i = 0-based POSITION
    use_turing = True
    for i: if use_turing and |lgt_rel_i - tur_rel_i|*(1+i)/tursd_i > 1.65: comb_i = tur_rel_i
           else: use_turing = False; comb_i = lgt_rel_i
    comb *= (1 - n[0]/N) / sum(comb*r*n/N); p0 = n[0]/N      # n[0] = count of the SMALLEST r, even if r[0] != 1
    return r*comb, p0
  ```
  Quirks to keep: positions instead of r values; z = 1.65; p0 uses `n[0]`; a non-consecutive r met
  while still in Turing mode gets Turing estimate 0.
  `estimate_profile_sgt`: p0 is split equally over features with zero ambient count but nonzero
  somewhere in the matrix. With no such features 3.0.2 divides by zero (crash); Space Ranger 1.0 /
  ≥ 3.1 and later renormalise `p_smoothed`.
- **Multinomial log-likelihood**: `lgamma(n+1) - Σ lgamma(x_g+1) + Σ x_g log p_g` over profile
  features (scipy `multinomial.logpmf` up to 8.x; sparse formula from 9.0).
- **Monte Carlo p-value**: `p = (1 + #{sims with loglik < obs}) / (1 + num_sims)`, simulated at
  the candidate's own total n.

### 3.1 `cellranger_v1_ordmag` (CR 1.x) [DOC only]

Zheng et al. 2017, Methods: all top barcodes within one order of magnitude (> 10 % of the top n-th
barcode, n = 1 % of the expected recovered cells) were called. Release note 1.2.0 mentions a
default cell count. Treat as `v2_ordmag` with N = 3000; the bootstrap is UNVERIFIED.

### 3.2 `cellranger_v2_ordmag` (CR 2.0–2.2) [SRC 2.0.2/2.1.0]

Functions: `stats.filter_cellular_barcodes_ordmag`, `find_within_ordmag`,
`summarize_bootstrapped_top_n`, `filter_cellular_barcodes_fixed_cutoff` (force-cells),
`filter_cellular_barcodes_manual`, called from `filter_barcodes()`.

Constants: `ORDMAG_RECOVERED_CELLS_QUANTILE = 0.99`, `ORDMAG_NUM_BOOTSTRAP_SAMPLES = 100`,
`DEFAULT_RECOVERED_CELLS_PER_GEM_GROUP = 3000`, `FILTER_BARCODES_MAX_RECOVERED_CELLS_MULTIPLE = 6`
(metric and assert only). Inputs: `--expect-cells` N (absent → 3000 per GEM group),
`--force-cells` (top N), `cell_barcodes` JSON (manual). Priority: manual > force > ordmag.

```
random.seed(0); np.random.seed(0)                # once per stage; stream shared across loops
for genome in genomes:                            # each genome's own features
  for gg in gem_groups:
    counts = UMIs per barcode (this genome, this gg)
    N = int(float(expect)/n_gg) if expect else 3000
    nz = counts[counts > 0]; if empty: no cells
    b = int(round(N*(1-0.99)))                    # Python-2 round; N*(1-0.99) is slightly > N/100
    b = min(b, len(nz)-1)
    top = [find_within_ordmag(np.random.choice(nz, len(nz)), b) for _ in range(100)]
       # find_within_ordmag: xs = sort(x); baseline = xs[-b]        <- 30th largest for N=3000; b==0 -> xs[0]!
       #                     cutoff = max(1, round(0.1*baseline)); return len(x) - searchsorted(xs, cutoff)
    n = round(mean(top))                          # 2.1+: int(round(...))
    cells = np.sort(np.argsort(counts)[::-1][:n]) # default quicksort: tie order arbitrary
force: top min(force//n_gg, #nonzero) by count.  manual: listed ∩ present.
```

Each genome's filtered matrix gets its own barcode set (barnyard). Metrics: lb/ub =
round(norm.ppf(0.025/0.975, mean, sd)), cv = sd/mean.

### 3.3 `cellranger_v3_emptydrops` (CR 3.0.0–3.0.2; 3.1 ≈ same) [SRC 3.0.2 + Space Ranger 1.0]

Step 1 = `v2_ordmag`, per genome and GEM group, GEX features only. Step 2 =
`cell_calling.find_nonambient_barcodes`, per GEM group, **genome-agnostic**, on the full GEX
matrix. Constants: `N_PARTITIONS=90000`, `N_CANDIDATE_BARCODES=20000`, `MIN_UMIS=500` (100 before
Sept 2018, pre-3.0.0), `MIN_UMI_FRAC_OF_MEDIAN=0.01`, `MAX_ADJ_PVALUE=0.01`, `num_sims=10000`,
`MAX_OCCUPIED_PARTITIONS_FRAC=0.5` (unused).

```
orig = ∪_genomes OrdMag calls for gg;  umis = colsums
desc = np.argsort(umis)[::-1]                      # quicksort, ties implementation-defined
empty = sort(desc[45000:90000]); ambient = empty ∩ {umis>0}
if ambient empty or SGT error: return (no extra cells)
F = features with rowsum>0 in this gg matrix;  p = SGT profile from ambient columns
min_umis = int(max(500, round(ceil(median(umis[orig]) * 0.01))))
cand = argsort(masked umis)[:n_unmasked][-20000:]  # not in orig, umis >= min_umis; top 20000 by UMI
obs = multinomial logpmf(M[F, cand]; umis, p)
distinct_n, sim = simulate_v1(p, umis[cand], 10000)
pval = (1+#sim<obs)/(1+10000); padj = BH; nonambient = padj <= 0.01
genome of each rescued bc = argmax per-genome UMI      # BUG until 6.0: counts and barcode strings misaligned
```

`simulate_v1` (3.0–8.x), called once per stage; RNG is the global `np.random` (3.x: seeded at
stage start; 7.2/8.x: `np.random.seed(0)` inside the function):

```
distinct_n = flatnonzero(bincount(umis_cand)); block = np.random.choice(G, 1_000_000, p=p); k = 0
for s in range(num_sims):
  x = multinomial.rvs(distinct_n[0], p); ll = logpmf(x); L[0,s] = ll
  for i>=1: step = n_i - n_{i-1}
    if step >= 1000: x += multinomial.rvs(step, p); ll = logpmf(x, n_i)
    else: for n in n_{i-1}+1..n_i: j = block[k]; k += 1 (refill at 1e6); x[j] += 1; ll += log p_j + log(n/x[j])
          # 7.2+: if 20 <= step < 1000, add `step` block draws at once:
          #       ll += Σc_f log p_f + lgamma(n_i+1) - lgamma(n_{i-1}+1) - Σ[lgamma(x_f+c_f+1) - lgamma(x_f+1)]
    L[i,s] = ll
```

Antibody aggregates (3.0), before step 1 if an Antibody library is present: remove a barcode if
`umi_corrected_reads/reads > 0.5` and `reads > 10000` (needs per-barcode read counts, not in the
matrix). 3.1 (Space Ranger 1.0 snapshot): identical, except the SGT `n0 == 0` case renormalises.

### 3.4 CR 4.0–6.1 [DOC; code UNVERIFIED]

- 4.0: targeted gradient (described from 7.2 code below); Feature-Barcode-only runs use step 1 on
  Antibody counts.
- 5.0: protein aggregation filter rebuilt (`cellranger_ab_aggregates_v5`, §3.9).
- 6.0: species annotation of EmptyDrops-rescued cells fixed (the misalignment above); a raw-matrix
  change slightly altered background.
- 6.1: "changed certain parameters in the cell calling algorithm"; EmptyDrops for multispecies
  experiments uses each species' UMI counts separately.

Present by 7.2 without a dated version: removal of the 20000 candidate cap, the LT range of 9000,
the changed OrdMag index. Exact versions unverifiable.

### 3.5 `cellranger_v7_ordmag_auto` (step 1; 7.0–10.1) [SRC 7.2.0, identical through 10.1 apart from typing]

`cell_calling_helpers.py`: `filter_cellular_barcodes_ordmag`, `estimate_recovered_cells_ordmag`,
`find_within_ordmag`, `summarize_bootstrapped_top_n`. Constants:
`ORDMAG_NUM_BOOTSTRAP_SAMPLES=100`, `ORDMAG_RECOVERED_CELLS_QUANTILE=0.99`,
`NP_SORT_KIND="stable"`, `MIN_RECOVERED_CELLS_PER_GEM_GROUP=50`,
`MAX_RECOVERED_CELLS_PER_GEM_GROUP=1<<18`.

```
rs = np.random.RandomState(0)                           # fresh per (genome, gg) call
nz = counts[counts>0]
find_within_ordmag(x, b): xs = sort(x); baseline = xs[-(b+1)]          # (b+1)-th largest; was xs[-b] in v2/v3
                          cutoff = max(1, np.round(0.1*baseline)).astype(int)   # banker's rounding
                          return len(x) - searchsorted(xs, cutoff)
if expect-cells absent:
   maxE = min(get_empty_drops_range(chem, nProbeBC)[0], 2**18) if chem else 2**18    # 7.0: always 2**18 [DOC]
   grid = unique(round(2**linspace(1, log2(maxE), 2000)).astype(int))              # 2 .. maxE
   per bootstrap s = rs.choice(nz, len(nz)) (100 of them):
       b = minimum(round(grid*(1-0.99)).astype(int), len(nz)-1); f = find_within_ordmag(s, b)
       loss = (f-grid)**2/grid; keep (grid[argmin], loss[argmin])
   N = max(int(np.round(mean of 100 grid picks)), 50)
else: N = max(expect//n_gg, 50)
b = min(int(np.round(N*(1-0.99))), len(nz)-1)
top = [find_within_ordmag(rs.choice(nz, len(nz)), b) for 100]     # continues the same rs stream
n = int(np.round(mean(top)))
# tie extension, as coded:
s = nz sorted desc (stable argsort reversed); cutoff = s[n-1]; idx = n-1
7.2:  while idx+1 < len(s) and s[idx] == cutoff: idx += 1; if idx+1-n > 0.2*n: return n
      n_final = idx+1
8.0+: same loop, but sets n_final = idx+1 inside the loop before testing again
      (on overflow it returns the partial value, not n)
cells = sort(argsort(counts, kind="stable")[::-1][:n_final])
```

As written, the tie-extension loop always runs once, so `n_final = n + 1` whenever
n < #nonzero, even with no ties (checked with a minimal replica: no ties, n = 5 → 6). The ≤ 7.2 and
≥ 8.0 variants differ only in the overflow case. To confirm against 10x public datasets.

7.1 release notes: grid cap became the ambient lower bound; "fixed a bug in OrdMag… all barcodes
called when very few cells" (the `xs[-0]` case). Which 7.x introduced the +1 index, tie extension
and 50 floor is unverifiable. Docs (2026): grid "2 to ~45,000 (NextGEM), ~80,000 (GEM-X)". Unknown
chemistry: cap 262144. Antibody-only runs: `ORDMAG` (no step 2) on Antibody features.

### 3.6 `cellranger_v7_emptydrops` (7.x; [SRC] 7.2.0)

Changes from v3:

- EmptyDrops **per genome and per GEM group**; features = that genome's GEX features;
  `orig` = union of all genomes' initial calls for the GEM group. A barcode can be added for
  several genomes.
- Ambient range `get_empty_drops_range` → (N/2, N): LT (`SC3Pv3LT`) 9000; multiplexed Flex
  (nProbeBC > 1) 45000 × nProbeBC; otherwise 90000.
- Candidates: every non-initial barcode with
  `umis ≥ int(max(emptydrops_minimum_umis, ceil(0.01·median(orig))))`; default 500, settable via
  multi-config `emptydrops-minimum-umis` (parsed in 7.2, documented from 9.0). **No 20000 cap.**
  Candidates sorted by index.
- `simulate_v1`, `NUM_SIMS=10000`, FDR 0.01.

Pipeline (`filter_barcodes/__init__.py`): method precedence `cell_barcodes` → MANUAL;
`force_cells` → TOP_N; override → that mode; antibody-only → ORDMAG; hybrid-capture targeted (not
RTL) → TARGETED; else ORDMAG_NONAMBIENT. Without expect-cells: None for OrdMag methods (auto),
3000 otherwise. Per-sample chunks for Flex (probe barcode) and OCM (overhang barcode): raw matrix
subset to barcodes whose substring `[offset:offset+len]` is in the sample's tags;
`num_probe_barcodes` = tags of that sample; expect/force/min-UMI per sample. The aggregate filter
(§3.9) runs first. In join: the high-occupancy GEM filter (§3.9); targeted non-spatial runs drop
cells with 0 target-gene UMIs.

### 3.7 `cellranger_v8_emptydrops` (8.0.0, 8.0.1) [SRC]

- Candidate minimum `max(emptydrops_minimum_umis, 1 + max(umis[empty_bcs]))`; the 1 %-of-median
  rule is gone; candidates lie above the ambient range.
- FDR 0.001 if chemistry in `SC3P_V4 + SC5P_V3`, else 0.01.
- Ranges: LT 9000; v4 / 5′v3 160000, or 80000 × nProbeBC if > 1; 8.0.1 adds HT (`SC3Pv3HT`,
  `SC5PHT`) 160000 (only when chemistry is set manually); else 90000, or 45000 × nProbeBC if > 1.
- New post-filters (not for MANUAL or TOP_N), no-ops by default:
  `apply_global_minimum_umis_threshold` (keep UMIs ≥ `global-minimum-umis`, default 0) and
  `apply_mitochondrial_threshold` (drop mito % > `max-mito-percent`, default 100.0; 13 human and
  13 mouse Ensembl IDs, matched after splitting the ID on "_"). Undocumented multi-config keys.
- Force-cells metric indexing fixed.

### 3.8 `cellranger_v9_emptydrops` (9.0.x) and `cellranger_v10_emptydrops` (10.0.0, 10.1.0) [SRC]

9.0:

- `NUM_SIMS=100000`.
- FDR 0.001 for `SC3P_V4 + SC3P_V3 + SC5P (incl. 5′ v2 and v3) + HT`; else 0.01 (Flex, LT, 3′v2,
  ARC, SC-FB…).
- Ranges: LT 9000; HT / 3′v4 / 5′v3 no tag 160000; HT / 3′v4 / 5′v3 with OCM tags 40000 × nTags;
  Flex v1 multiplexed 90000 × nProbeBC; else 90000.
- Sparse log-likelihood and new simulator:
  ```
  RNG = np.random.default_rng(seed=42)    # created at module import; stream continues across genome/gg calls
  max_n = max(cand umis); pcum = cumsum(p)
  per sim: draws = searchsorted(pcum, RNG.random(max_n))
           rank_t = occurrence number of draws[t] among draws[0..t] (1-based)
           ll = cumsum(log t - log rank_t + log p[draws_t]); compare ll[n_b - 1] < obs_b
  ```
- A Dirichlet-multinomial option exists in 9.0/10.0 (`method="dirichlet"`) but the pipeline never
  passes it; removed in 10.1.
- EmptyDrops output has a `sample` column; loop order gg then genome.

10.0 / 10.1:

- 10.0: Flex v2 N = 160000 × (nProbeBC or 1); features sorted by feature ID before profile
  estimation (changes the simulation stream); TARGETED no longer auto-selected; high-occupancy
  filter uses NB.
- 10.1: LT special case removed (LT now 90000); profile kept in matrix order but permuted to
  feature-ID order for simulation only (same effect as 10.0); hidden multi-config option
  `cell-calling-mode ∈ {ordmag, ordmag_nonambient, gradient}`; `cell_calling_stats.py` replaces
  `stats.py`.

### 3.9 Ancillary Cell Ranger rules

**`cellranger_targeted_gradient`** [SRC 7.2]. Constants: `N_CANDIDATE_BARCODES_GRADIENT=20000`
(2000 for LT, LT removed in 10.0), `TARGETED_CC_MIN_UMIS_ADDITIONAL_CELLS=10`,
`TARGETED_CC_MIN_UMIS_FROM_TARGET_GENES=1`.

1. N = expect or 3000, floored at 50.
2. On descending nonzero counts: `base = s[round(N*0.01)]`; `lower = #(s ≥ base/10) - 1`;
   `upper = clip(min(lower+maxAdd, #(s≥10)), lower, len-1)`.
3. Points: for each unique count u, x = log10(#(s≥u)), y = log10(u); append (log10(1+Σs), 0).
4. `UnivariateSpline(k=min(3,len-1), s=0)`; with > 50 points refit `LSQUnivariateSpline` on
   `get_spline_num_knots(n)` interior knots from the original knots by
   `linspace(1, len-2, n-2, int)`.
5. Knots: n<50 → n; <200 → 2^(log2 50+(log2 100−log2 50)(n−50)/150); <800 →
   2^(log2 100+(log2 140−log2 100)(n−200)/600); <3200 → 2^(log2 140+(log2 200−log2 140)(n−800)/2400);
   else 200+(n−3200)^0.2.
6. Derivatives at the points, zeroed outside rank ∈ [lower, upper]; `cut = round(10**y[argmin])`;
   `n = max(#(s > cut), lower+1)`.
7. TARGETED keeps cells with ≥ 1 target-gene UMI.

**`cellranger_ab_aggregates_v3`** (3.0–4.x?): the read-based "highly corrected" rule above.

**`cellranger_ab_aggregates_v5`** (5.0+; [SRC] 7.2–10.1), before cell calling unless disabled:

- Rule A: `frac_corrected_reads > 0.5` and `reads > 10000` (needs read data).
- Rule B (counts only): drop antibodies with total UMIs < 1000
  (`BACKGROUND_ANTIBODY_UMI_THRESHOLD`); skip if fewer than 5 remain; K = (nProbeBC or 1) × 25
  (`TOP_UMI_BCS`); candidates = top-K barcodes by total antibody UMI; aggregate if in the top-K for
  at least `round(nSignal × f)` antibodies, f = −0.02·nSignal + 1.1 (0.6 if nSignal > 26).
- Antigen (BEAM): among the top 100 barcodes by antigen UMI, remove those with UMI ≥ Q3 + 3·IQR
  (only if that threshold ≥ 1000).

**`cellranger_flex_hog_poisson`** (7.1?–9.x; [SRC] 7.2) and **`cellranger_flex_hog_nb`** (10.0+),
multiplexed Flex only. Constants `TOTAL_INSTRUMENT_PARTITIONS=115000`, `RECOVERY_FACTOR=1/1.65`.
GEM barcode = `[:16]` + probe barcode; `total_gems = max(int(115000/1.65), #GEMs with cells)`;
λ = cells/total_gems. Simulation (`np.random.seed(0)`): 1,000,000 GEMs, Poisson(λ) (7.2) or
NB(n=λ, p=0.5) (10.0+); keep occupied; assign each cell a probe barcode with observed frequencies;
count distinct probe barcodes per GEM. Threshold: 7.2 `ceil(quantile(counts, 0.999))`; 10.0
smallest k with cumulative fraction of cells (k·freq) ≥ 0.999. Remove every cell in a GEM above it.

**TOP_N / force-cells**: `min(force//n_gg, #nonzero)` top barcodes by stable argsort.
**MANUAL**: a list of barcodes.

## 4. Cell Ranger ATAC and ARC

### `cellranger_atac_v1_nbmix` [SRC 1.2]

Source: `mro/atac/stages/processing/cell_calling/{remove_low_targeting_barcodes,
remove_gel_bead_doublet_barcodes, remove_barcode_multiplets, detect_cell_barcodes}/__init__.py`,
`lib/python/analysis/em_fitting.py`, `lib/python/constants.py`.

Count used: fragments overlapping peaks, per barcode and species (not cut sites).

Exclusions (fragment-level, not reproducible from a matrix): low-targeting (1.1+;
frac_in_peaks < bp within peaks padded ±2000 / genome bp; 2.1+ docs say ±250), gel-bead doublets
(1.2; `MINIMUM_COUNTS=250`, `MAXIMUM_BARCODES=15000`), barcode multiplets (1.2; `MINIMUM_COUNT=100`,
`SELF_SIGNAL_THRESHOLD_MULTIPLIER=0.4`, `MAXIMUM_PILEUP=20`).

Constants: `CELL_CALLING_THRESHOLD=100000` (odds), `MAXIMUM_CELLS_PER_SPECIES=20000`,
`MINIMUM_COUNT=5`, `WHITELIST_CONTAM_RATE=0.02`.

```
shift = max(5, int(total_fragments*0.02/len(whitelist)))      # whitelist = 737K-cratac-v1
y = counts[counts >= shift] - shift;  cd = Counter(y);  if len(cd) < 10: 0 cells
init threshold: smallest-mass cut so that values <= t hold <= 1% of Σ value*mult (else values[1]);
     clamp to [values[1], values[-2]]
NB(mu, alpha) as scipy nbinom(n=1/alpha, p=1/(1+alpha*mu))
EM (maxiter 250, eps 1e-6 relative on all 5 params): responsibilities (columns where both pmfs are 0
   zeroed); mu = weighted mean; alpha = MLE dispersion (fixed point from MOM, floor 1e-6, maxiter
   10000, tol 1e-6; MOM if var <= mu); frac_noise; swap so mu_signal > mu_noise
threshold = first k in [max(int(mu_noise), nbinom.ppf(0.001, signal)), int(mu_signal)]
     with (1-f)·NB_s(k) / (f·NB_n(k)) >= 100000 (else int(mu_signal)); + shift
cells = counts >= threshold
if #cells > 20000, or --force-cells N (capped at 20000): threshold = N-th largest of y>0 + shift
```

Barnyard second pass is coded with an inverted condition, so it does not mask non-cells as
documented. Versions [DOC]: 1.0 no low-targeting filter; 1.1 adds it; 1.2 GB doublets and
multiplets; 2.0 new peak caller; 2.1 less stringent low-targeting, improved multi-species; 2.2
lifts the 20000 limit. No 2.x source.

### `cellranger_arc_joint_v1/_v2` [DOC only]

Per barcode: ATAC transposition events in peaks and GEX UMIs. Prefilters: ATAC low-targeting,
gel-bead doublets, > 1 count in each library. De-duplicate identical count pairs; OrdMag per axis
(cell if above both); refine with k-means K=2 from the two groups' centroids; map labels back;
cap < 20k cells. OrdMag N, log scaling and k-means metric undocumented. Override
`--min-atac-count X --min-gex-count Y` (AND in 2.0, OR in 1.x).

## 5. Exactness caveats

- RNG: v2/v3 numpy legacy `RandomState`, global, seeded 0 once, consumed in loop order (genomes →
  GEM groups → then the simulation). 7.2+: OrdMag `RandomState(0)` per call; simulation re-seeds
  the global with 0. 9.0+: PCG64(42) created once per process.
- Sort ties: ambient selection uses `np.argsort(umis)` (quicksort) whose tie order can differ
  between numpy versions and SIMD hardware. OrdMag selection uses stable sort from 7.2.
- Rounding: Python 2 `round` (half away from zero) up to 3.x; `np.round` (half to even) in 7.2+.
- Not in a count matrix: read-based antibody rule A, ATAC exclusions and fragments, the per-sample
  tag map (supply it), the chemistry name (supply it).

## 6. STARsolo

Sources (`github.com/alexdobin/STAR/blob/<tag>/`): `source/SoloFeature_cellFiltering.cpp`,
`source/SoloFeature_emptyDrops_CR.cpp`, `source/SimpleGoodTuring/sgt.h`,
`source/ParametersSolo.cpp` (`ParametersSolo::cellFiltering`), `source/parametersDefault`,
`source/SoloFeature_statsOutput.cpp`, `source/SoloFeature_loadRawMatrix.cpp`, `source/Solo.cpp`,
`CHANGES.md`, `RELEASEnotes.md`, `docs/STARsolo.md`, `extras/scripts/soloBasicCellFilter.awk`.

Common: `nUMIperCB` = UMIs per detected CB for the main (first) `--soloUMIdedup`, counting
unique-gene reads. Sorting descending; EmptyDrops breaks ties by ascending CB index. Applies to
Gene, GeneFull, GeneFull_Ex50pAS, GeneFull_ExonOverIntron and loaded matrices; Velocyto reuses the
Gene calls (2.7.8a+). Other features are not filtered.

| Rule | STAR versions | Definition |
|---|---|---|
| `starsolo_none` | 2.7.0a–2.7.2d (no option); `--soloCellFilter None` from 2.7.3a | — |
| `starsolo_basic_awk` | script added 2.7.2a | `nMax` = asc-sorted[len − int((1−p)·maxCells)] (1-based, the 31st largest for 3000/0.99); `nMin = nMax/maxMinRatio` (float); cells: total ≥ nMin; `exactCells` → nMin = sorted[len − exactCells]. Defaults 3000 / 0.99 / 10. |
| `starsolo_cellranger2_2_v1` | 2.7.3a–2.7.7a (default `CellRanger2.2 3000 0.99 10`; the numbers required) | `ind = (uint64)(N·(1−p)+0.5)` → 30; `nUMImax = sorted[min(nCB−1, ind)]` (31st largest; Cell Ranger 2.x used the 30th and bootstrapped); `nUMImin = int(nUMImax/ratio + 0.5)`; cells: nUMI ≥ nUMImin. No bootstrap, no floor. |
| `starsolo_cellranger2_2_v2` | 2.7.8a+ (numbers optional) | Same with `int(round(...))`, plus `nUMImin = max(nUMImin, 1)`. |
| `starsolo_topcells_buggy` | 2.7.3a–2.7.5c | `nUMImin = sorted[max(nCB−1, T)]`: calls all CBs if T < nCB−1; out of bounds otherwise. |
| `starsolo_topcells` | 2.7.6a+ (PR #1012) | `nUMImin = sorted[min(nCB−1, T)]`; cells: nUMI ≥ that, so ≥ T+1 cells plus ties (off by one as coded). |
| `starsolo_emptydrops_cr` | 2.7.8a – 2.7.11b / master; unchanged except a log message (2.7.10a) | Below. |

**`starsolo_emptydrops_cr`**: `EmptyDrops_CR nExpectedCells maxPercentile maxMinRatio indMin
indMax umiMin umiMinFracMedian candMaxN FDR simN` = `3000 0.99 10 45000 90000 500 0.01 20000 0.01
10000`; give 0 or all 10.

```
simple = starsolo_cellranger2_2_v2 (nCellsSimple)
if nCB <= indMin: stop
featDet = genes with >0 in any CB
amb[g] = Σ counts over CBs at desc ranks [indMin, min(nCB, indMax))
freq-of-freq over amb; remove undetected genes from the zero class
SGT (Elworthy/Sampson C): needs >= 5 distinct r; PZero = n_1/N; Z = 2n_r/(next − prev);
   OLS fit; switch to LGT at first gap or when |x−y| <= 1.96·sqrt((r+1)²·n_{r+1}/n_r²·(1+n_{r+1}/n_r));
   est_r = (1−PZero)·r*/Σ n_r r*
p_g = est[amb_g]; zero-class genes get PZero/#zeroGenes; normalise over featDet
candidates: desc ranks nCellsSimple .. nCellsSimple+candMaxN−1, stopping at count < minUMI,
   minUMI = max(umiMin, int(umiMinFracMedian · sorted[nCellsSimple/2]))   # truncation
obs: log multinomial pmf with a log-factorial table
sims: for s in 0..simN−1: std::mt19937 seeded 19760110*(s+1); std::discrete_distribution over p;
      walk 1..maxCount accumulating ll += log p_g + log c − log count_g
p = (1+#{sim[count] < obs})/(1+simN); BH (no min(1)); add if padj <= FDR
```

Differences from Cell Ranger 3.0: deterministic threshold at the 31st-largest count instead of the
100-bootstrap mean; median rule truncates the upper-median element instead of ceil of
`np.median`; SGT variant (z 1.96 vs 1.65, r values vs positions, ≥ 5 vs ≥ 10 distinct
frequencies, no slope check, PZero = n₁/N); independent per-sim MT19937 streams, single draws, up to
the maximum count (Cell Ranger: nested increments over a pre-drawn block); tie order; no genome
split, aggregate or HOG filters, no tie extension. The draw order of `discrete_distribution` is
library-specific (libstdc++ vs libc++).

Summary.csv "Estimated Number of Cells" is `filteredCells.nCells`, the filtered set. There is no
separate knee in STAR. 2.7.3a–2.7.7a print the cell block only for CellRanger2.2; from 2.7.8a for
any non-None filter. `UMIperCellSorted.txt` lists sorted UMIs; `Features.stats` holds read-level
stats; `CellReads.stats` (`--soloCellReadStats Standard`, 2.7.10a+) holds per-CB counts.

**`--runMode soloCellFiltering <raw_dir>/ <out_prefix>`** (first in 2.7.8a): needs no genome or
`--soloType`; only `--soloCellFilter` matters (default CellRanger2.2; pass `EmptyDrops_CR`
explicitly). Reads `<raw_dir>/` + `--soloOutFileNames` entries (`matrix.mtx`, `barcodes.tsv`,
`features.tsv`), **uncompressed**; `%` header lines skipped; coordinate MatrixMarket, rows =
features, columns = barcodes; values read as double and `std::round`ed; `nTot == 0` is an error.
`nCB` = columns with ≥ 1 entry. Writes filtered `barcodes.tsv`, `features.tsv`, `matrix.mtx` to
`<out_prefix>`. Bug #1230 segfaults this mode in 2.7.9a; use 2.7.10a or later.

## 7. Cannot be verified from source

Cell Ranger 1.x code; exact 3.1–7.1 code (removal of the 20k cap, LT 9000, the `-(b+1)` index,
the tie extension and the 50 floor, all present by 7.2; the 7.0 grid); whether 4.0–6.x
Feature-Barcode-only runs used N = 3000; cellranger-atac 2.x code; all cellranger-arc internals;
10.1's Flex v2 HOG partition constants (the code still uses 115000/1.65 for every chemistry).
