<!-- Clean-room specification, written by reading Cell Ranger 7.2.0-10.1.0 and describing
behaviour only (no code). The implementation in src/sc_droplets/methods/cellranger.py is written
from this text and the method pages, never from 10x's source. -->

# Cell Ranger cell calling, 7.2.0 to 10.1.0: behavioural specification (clean-room)

I read the source of every listed tag and describe behaviour only: no code, and no function bodies copied. 9.0.0 and 9.0.1 are byte-identical in every relevant file, so they appear below as "9.0.x". The only cell-calling difference between 8.0.0 and 8.0.1 is the HT ambient range (plus a metric fix that does not change calls). I downloaded the source only to my scratchpad; nothing was written under the project directory.

## 0. Conventions used throughout

- **round()** means round half to even, as numpy's round does, unless stated otherwise.
- **q** means `1 − 0.99` evaluated in float64, which is 0.010000000000000009 (slightly more than 0.01). It is used for every "baseline index". Because of it, products that should land exactly on .5 land slightly above and round **up**. For example, 50·q rounds to 1 and 250·q rounds to 3. Reproduce it with the float expression, not with the constant 0.01.
- **Stable-descending order** means an ascending *stable* argsort by count, then reversed. The result is descending by count, and **among equal counts the barcode with the larger column index comes first**.
- **Default argsort** means numpy's default `quicksort`, which is not stable. Its tie order depends on the implementation: numpy version, SIMD dispatch and dtype (see §12).
- **Barcode order** is the order of the matrix columns. Raw matrices are sorted by (GEM group, valid flag, sequence bytes); with one GEM well this is lexicographic order of the barcode string.
- **npb** is the number of probe-barcode or OCM tags assigned to the sample being processed. It is None when there is no per-sample split.

---

## 1. Chemistry groups → EmptyDrops ambient range and FDR

### How the lookup works (all versions)

- Membership is decided by an **exact match on the chemistry's `description` string**, not on its name. Any two chemistries with the same description therefore behave identically. Examples: in 7.2, SFRP and MFRP are both "Fixed RNA Profiling"; in 9.0+, SFRP and MFRP-RNA are both "Flex Gene Expression".
- A number of partitions N is chosen, and the ambient range is `(low, high) = (N // 2, N)` using integer floor division.
- `low` is also the cap on the OrdMag recovered-cells grid (§3).
- The rules are tested in the order listed for each version. The first rule that matches wins.

### 7.2.0

| Rule | N | FDR |
|---|---|---|
| Description is "Single Cell 3' v3 LT" (SC3Pv3LT) | 9000 | 0.01 |
| Otherwise, if npb > 1 | 45000·npb | 0.01 |
| Otherwise | 90000 | 0.01 |

The FDR is the constant 0.01 in 7.2.

### 8.0.0 and 8.0.1

| Rule (in order) | N |
|---|---|
| LT (SC3Pv3LT) | 9000 |
| **8.0.1 only:** HT group = SC3Pv3HT ("Single Cell 3' v3 HT") and SC5PHT ("Single Cell 5' HT") | 160000, whatever npb is |
| "V4" group = SC3Pv4, SC3Pv4-OH, SC5P-PE-v3, SC5P-R2-v3, SC5P-R2-OH-v3, SC5P-R1-v3 | 80000·npb if npb > 1, else 160000 |
| Anything else | 45000·npb if npb > 1, else 90000 |

- FDR is 0.001 if the description is in the V4 group, otherwise 0.01. HT therefore stays at 0.01 in 8.0.1.
- SC3Pv4HT and SC5PHT-v3 exist only in 8.0.0, and they take the default rule.

### 9.0.x, 10.0.0 and 10.1.0

| Rule (in order) | N |
|---|---|
| LT (9.0.x and 10.0.0 only; the branch is removed in 10.1) | 9000 |
| Range group **R** = SC3Pv3HT-polyA, SC3Pv3HT-CS1, SC5PHT, SC3Pv4-polyA, SC3Pv4-CS1, SC3Pv4-polyA-OCM, SC3Pv4-CS1-OCM, SC5P-PE-v3, SC5P-R2-v3, SC5P-R2-OCM-v3, SC5P-R1-v3 | 160000 if npb is None, else **40000·npb** |
| 10.0+ only: Flex v2 group = every chemistry whose name starts with "Flex-v2" (96-R1, 96-RNA-R2, R1, RNA-R2, Ab-R2:45, Ab-R2:64, CRISPR-R2:1, singleplex) | 160000 × (npb, or 1 if npb is None or 0) |
| npb is not None ("Flex" branch) | 90000·npb |
| Otherwise | 90000 |

FDR group **F** gets FDR 0.001; everything else gets 0.01. F contains:
- the four 3'v4 chemistries;
- SC3Pv3-polyA and SC3Pv3-CS1;
- SC5P-PE, SC5P-R2 and SC5P-R1;
- SC5P-PE-v3, SC5P-R2-v3, SC5P-R2-OCM-v3 and SC5P-R1-v3;
- SC3Pv3HT-polyA, SC3Pv3HT-CS1 and SC5PHT.

Chemistries in neither R nor F, which take the default range and FDR 0.01: SC3Pv1, SC3Pv2, SC3Pv3-polyA-OCM, SC3Pv3-CS1-OCM, SC5P-R2-OCM (the non-v3 one), SC5P-PE-OCM-v3, SC5P-R1-OCM-v3, SFRP, all MFRP-\* (including MFRP-96-\* in 10.x), SC-FB, ARC, spatial and custom.

### Summary matrix

Notation:
- D7 = (45000, 90000) when npb is None or ≤ 1; (22500·m, 45000·m) when npb = m ≥ 2.
- V8 = (80000, 160000) when npb ≤ 1; (40000·m, 80000·m) when m ≥ 2.
- R9 = (80000, 160000) when npb is None; (20000·m, 40000·m) when npb = m (including m = 1).
- D9 = (45000, 90000) when npb is None; (45000·m, 90000·m) otherwise.
- F2 = (80000·m, 160000·m), with m = 1 when npb is None.
- LT = (4500, 9000).

Each cell gives range / FDR.

| Kit | 7.2.0 | 8.0.0 | 8.0.1 | 9.0.x | 10.0.0 | 10.1.0 |
|---|---|---|---|---|---|---|
| 3' v1/v2 | D7/.01 | D7/.01 | D7/.01 | D9/.01 | D9/.01 | D9/.01 |
| 3' v3 (SC3Pv3 → SC3Pv3-polyA/CS1) | D7/.01 | D7/.01 | D7/.01 | D9/**.001** | D9/.001 | D9/.001 |
| 3' v3 OH/OCM | D7/.01 | D7/.01 | D7/.01 | D9/.01 | D9/.01 | D9/.01 |
| 3' v3 HT | D7/.01 | D7/.01 | (80000,160000)/.01 | R9/.001 | R9/.001 | R9/.001 |
| 3' v3 LT | LT/.01 | LT/.01 † | LT/.01 † | LT/.01 † | LT/.01 † | n/a (gone) |
| 3' v4 (SC3Pv4 / -polyA / -CS1) | n/a | V8/.001 | V8/.001 | R9/.001 | R9/.001 | R9/.001 |
| 3' v4 OH/OCM | n/a | V8/.001 | V8/.001 | R9/.001 | R9/.001 | R9/.001 |
| 3' v4 HT (SC3Pv4HT) | n/a | D7/.01 | n/a | n/a | n/a | n/a |
| 5' v1/v2 (SC5P-PE/R2/R1) | D7/.01 | D7/.01 | D7/.01 | D9/**.001** | D9/.001 | D9/.001 |
| 5' R2 OH/OCM (non-v3) | D7/.01 | D7/.01 | D7/.01 | D9/.01 | D9/.01 | D9/.01 |
| 5' HT (SC5PHT) | D7/.01 | D7/.01 | (80000,160000)/.01 | R9/.001 | R9/.001 | R9/.001 |
| 5' v3 PE/R2/R1 | n/a | V8/.001 | V8/.001 | R9/.001 | R9/.001 | R9/.001 |
| 5' v3 R2 OH/OCM | n/a | V8/.001 | V8/.001 | R9/.001 | R9/.001 | R9/.001 |
| 5' v3 PE-OCM, R1-OCM | n/a | n/a | n/a | D9/.01 | D9/.01 | D9/.01 |
| Flex v1 (SFRP, MFRP\*, MFRP-96\*) | D7/.01 | D7/.01 | D7/.01 | D9/.01 | D9/.01 | D9/.01 |
| Flex v2 (all Flex-v2-\*) | n/a | n/a | n/a | n/a | F2/.01 | F2/.01 |
| SC-FB, ARC, spatial, custom, any unknown | D7/.01 | D7/.01 | D7/.01 | D9/.01 | D9/.01 | D9/.01 |

† From 8.0 onward the pipeline refuses LT ("no longer supported, use 7.2"), so the LT branch cannot be reached in a real run.

### Unknown or custom chemistry

- A custom chemistry gets a description of the form "custom: <index-scheme name>". That matches no group, so it gets the default rule and FDR 0.01.
- If the description were None (direct API use only), OrdMag would use the 262144 cap. `get_empty_drops_range` itself would still take the default branch.
- In the pipeline the description is never None.

### The 8.0.1 "HT only when set manually" condition

- **What the code tests:** only that the primary chemistry's description is exactly "Single Cell 3' v3 HT" or "Single Cell 5' HT". npb plays no part, and this rule is checked before the V4 rule.
- **Why it is "manual only":** chemistry auto-detection resolves the whitelist ambiguity {SC3Pv3, SC3Pv3HT} to SC3Pv3, and {SC5P-R2, SC5PHT} to SC5P-R2.
- So HT is in effect only when the user explicitly sets chemistry SC3Pv3HT or SC5PHT (`count --chemistry`, or the multi `chemistry` setting).
- 9.0.x to 10.1 keep the same down-selection to the regular kit, so the R9 range for HT there also needs the chemistry to be set explicitly.

---

## 2. Where the chemistry and npb come from

### Which chemistry description is used

- **7.2:** the GEM well has a single chemistry definition, either detected or set by the user. Its description string goes into cell calling.
- **8.0+:** the stage receives a map from library type to chemistry definition.
  - The "primary" definition is the only entry if there is one; otherwise it is the "Gene Expression" entry.
  - If there are several entries and none is Gene Expression, the stage fails.
  - The primary definition's description is used for range, FDR and OrdMag cap. Its barcode layout is used for the probe and overhang offsets.

### Per-sample split (the stage's split phase)

- A sample split happens only when a multi graph exists and is either RTL-multiplexed (Flex probe barcodes) or overhang/OCM-multiplexed.
- **RTL (Flex):**
  - offset = gel-bead barcode length (16);
  - length = probe barcode length (8 for Flex v1 "right_probe", 10 for Flex v2);
  - if there are both a left and a right probe, length runs to the end of whichever finishes later.
- **OCM/OH:** offset and length come straight from the chemistry's "overhang" entry. For the 3'v4 and 5'v3 OCM chemistries this is offset 7, length 2, i.e. positions 7–8 of the 16-nt cell barcode.
- **Per sample:**
  - The sample's tag IDs (probe_barcode_ids or OCM IDs from the multi config) are translated into sequences through the whitelist's ID→sequence map.
  - **npb = the number of tags assigned to that sample** — not the total number in the run.
  - The raw matrix, and the barcode-correction table, are restricted to barcodes whose `barcode[offset : offset+length]` is one of those sequences.
  - Each sample runs as its own chunk, in its own process.
- Every other case (no multi graph, CellPlex/CMO, singleplex Flex) runs as one chunk with npb = None.

### Values that reach cell calling

| Case | Description string | npb | Ambient range / FDR |
|---|---|---|---|
| (a) Plain 3' v3 | 7.2–8.0.x: "Single Cell 3' v3"; 9.0+: "Single Cell 3' v3 (polyA)" (the GEX library's definition) | None | (45000, 90000) every version. FDR 0.01 up to 8.0.x, 0.001 from 9.0 |
| (b) Multiplexed Flex, per sample | GEX library's definition. 7.2/8.x: "Fixed RNA Profiling" or "Fixed RNA Profiling (probe barcode on R1)". 9.0+: "Flex Gene Expression" or "Flex Gene Expression (probe barcode on R1)". Flex v2: e.g. "GEM-X Flex v2 Gene Expression" | Number of probe barcodes for that sample, usually 1 | v1, 1 BC: (45000, 90000) in every version. v1, 2 BCs: 7.2/8.x (45000, 90000); 9.0+ (90000, 180000). Flex v2, 1 BC: (80000, 160000). FDR 0.01 |
| (c) OCM (9.0+) | e.g. "Single Cell 3' v4 (polyA) OCM" | Number of OCM tags for that sample, usually 1 | 1 tag: (20000, 40000), FDR 0.001 |

- In 8.0.x the overhang analogue ("Single Cell 3' v4 OH") with 1 tag gives (80000, 160000).

### Other inputs

- **GEM groups:** always one per chunk (`[gem_well]`; written as `[1]` in 10.1).
- **Recovered/force/ED-min/global-min/mito settings:** a library-level ("per GEM well") value wins and is applied unchanged to every sample chunk. Otherwise the per-sample value is used.
- **Expect-cells and force-cells with no sample** (a single chunk): the per-sample expect values are summed, and if any is missing the result is None. Force-cells is not summed.
- **Per-sample value of 0 treated as unset (quirk):** the fallback tests truthiness, so a per-sample 0 becomes the default — ED-min becomes 500 and max-mito becomes 100, and expect/force become "not provided". A library-level 0 is honoured.
- **Genome list:**
  - 7.2–9.0.x: the reference's metadata list, in file order. In 9.0 this falls back to the probe set's `reference_genome`, then to ["NONE"].
  - 10.0+: the **sorted** set of non-empty genome tags in the matrix's feature reference, or ["NONE"].

---

## 3. OrdMag (initial calls)

### Inputs

- Counts c for all barcodes of the GEM group, in column order.
- Features: those of the selected feature types, with Gene Expression features limited to the current genome. Non-GEX features of a selected type are kept whatever their genome.
- nz = the counts that are > 0, in column order. If nz is empty, no cells are called.

### Random generator

- Each OrdMag call (one per genome × GEM group) creates a **fresh, private legacy numpy `RandomState(0)`** (MT19937).
- Every bootstrap draw is that object's `choice(nz, size=len(nz))`: uniform, with replacement, no weights.
- The draws come in order: 100 for estimation (only when no expected count is supplied), then 100 for the final call, continuing the same stream.

### Expected cells supplied

1. The value is an int (8.0+ truncates a float with `int()`).
2. The stage computes `int(float(E) / G)`, where G is the number of GEM groups (always 1).
3. The helper then does `E // G` again if E is truthy; 0 counts as "not supplied".
4. E' = max(E, 50).
5. No estimation bootstrap is run.

### Auto-estimation (no expected cells)

1. Cap: M = min(low end of the ambient range from §1 for this chemistry and npb, 262144).
2. Grid of exponents: 2000 evenly spaced float values from 1 to log2(M) inclusive (numpy linspace).
3. Candidate cell counts: each is round(2^exponent) as an int. Drop duplicates and sort ascending; call this list R.
4. For each of the 100 bootstrap samples S:
   - For each r in R:
     - i = min(int(round(r·q)), len(nz) − 1);
     - base = the (i+1)-th largest value of S;
     - cut = max(1, round(0.1·base));
     - f(r) = number of values in S that are ≥ cut (searchsorted, left side, on the ascending copy);
     - loss(r) = (f − r)² / r. The square is taken in integers; the division is float.
   - Pick the **first** index with the minimum loss. Since R ascends, ties go to the smallest r.
   - Keep that r.
5. E_est = mean of the 100 kept r values (float), rounded, converted to int, then max(·, 50).

### Final call

1. i = min(int(round(float(E′ or E_est)·q)), len(nz) − 1).
2. Draw 100 more bootstrap samples from the same generator. For each, f_b = the number of values ≥ max(1, round(0.1 × (i+1)-th largest)).
3. n0 = int(round(mean(f_b))). The population variance is computed too, but only for metrics.

### Tie extension (applied when n0 > 0)

Let s = nz sorted descending (the method is stable-descending, but only the values matter here). Let c = s[n0 − 1].

Shared loop:
- Start with an index j = n0 − 1.
- While j + 1 < len(s) **and s[j] = c**: increase j by 1, then test whether (j + 1 − n0) > 0.2·n0, with 0.2·n0 computed in float64.

The test is made on s[j] *before* the increment. So the barcode at the new j is included without being checked. The loop therefore always takes **one barcode past the end of the tie run**, i.e. the first barcode with a count below c. With no ties and n0 ≥ 5, the result is n0 + 1. For n0 ≤ 4 the very first step exceeds 20%.

| Version | When the 20% test fires | Normal exit |
|---|---|---|
| 7.2.0 | Revert: return **n0** | Result j + 1. A guard "c > 0" is always true |
| 8.0+ | Return the **last accepted value**: the extension is truncated at the largest m with m − n0 ≤ 0.2·n0, or n0 if the first step already fails | Result j + 1, updated inside the loop after each passed test. The "c > 0" guard is removed; no behavioural change |

- If the run reaches the end of the array, the result is len(s).
- 8.0+ also records c as a metric only.

### Selection

- Take the first top_n barcodes of the **stable-descending order over all of c** (zero counts included). Among tied counts, higher column index wins.
- Return them sorted by index.
- In practice the single below-c barcode added by the off-by-one is the one with the highest index among barcodes having that next-lower count.

### Sort kinds used

| Where | Sort |
|---|---|
| Estimation and bootstrap | Default `np.sort` of values (ties irrelevant) plus left searchsorted |
| Tie extension | Stable argsort reversed (values only) |
| Final selection | Stable argsort reversed, then ascending sort of indices |

---

## 4. EmptyDrops in 7.2 and 8.x

### Inputs

Per (genome, GEM group): the genome's GEX features (§6), and the per-barcode UMI totals u over those features.

### Ambient set (unchanged up to 10.1)

1. Order all barcodes by the **default (non-stable) argsort ascending of u**, then reverse.
2. Take ranks [low, high), 0-based, of that descending list. This is the "empty" set.
   - Ties at the two boundaries follow the unstable sort and are implementation-defined (§12).
3. Ambient = empty set ∩ {u > 0}.
4. If there are fewer than `low` barcodes, the ambient set is empty and EmptyDrops ends with no calls.

### Profile (unchanged up to 10.1)

1. **Evaluated features** = features with a non-zero total over **all** barcodes of this genome and GEM group (not only the ambient ones).
2. Profile counts = the sum over ambient barcodes of each evaluated feature.
3. Simple Good–Turing (SGT) is applied to the positive counts. This is a direct port of Gale & Sampson's published S code:
   - **Abort** (EmptyDrops returns nothing for this genome/GEM group) if there are fewer than 10 distinct frequency values.
   - Let r = the distinct frequencies (ascending) and n_r = their counts.
   - Averaging transform: Z_r = n_r / dr, where dr is half the sum of the gaps to the previous and the next r. The **gap before the first r is fixed at 1**, even if the smallest r is not 1. For the last r, dr is its own preceding gap.
   - The slope b is the ordinary least-squares slope of log Z on log r.
   - **7.2/8.x:** if b > −1, abort.
   - **9.0+:** if b ≥ −1, set b = −1 and continue.
   - Smoothed (linear Good–Turing) relative estimate: (1 + 1/r)^(1 + b).
   - Turing relative estimate: (r+1)/r · n_{r+1}/n_r, but only where the next observed frequency is exactly r+1; otherwise 0.
   - Switching rule, walking the frequencies by 0-based index j:
     - while "use Turing" is on, compare |smoothed − Turing| · (j+1) / sd_j with 1.65;
     - sd_j = (j+2)/n_j · sqrt(n_{j+1}·(1 + n_{j+1}/n_j)) where the Turing estimate exists, and 1 otherwise;
     - if the comparison is > 1.65, use the Turing estimate; otherwise switch to the smoothed estimate permanently.
     - These index-based factors are faithful to the S code.
   - p0 = n_first / N, where n_first is the count of the *smallest* observed frequency and N = Σ r·n_r.
   - Each feature with frequency f gets (1 − p0) · r*(f) / Σ_r n_r·r*(r), where r* = r × the chosen relative estimate. The intermediate renormalisation cancels out.
4. Zero-count evaluated features: if there are n0 > 0 of them, each gets p0/n0. If there are none, the observed-feature probabilities are divided by their own sum and p0 is discarded.

### Candidates

**7.2:**
- minU = int(max(E_min, ceil(median_u × 0.01))), where E_min is the configured value (default 500).
- median_u is the numpy median (mean of the two middle values when even) of u over the "original" cells. Compute median × 0.01 in float64 before the ceiling.
- Candidates are barcodes that are not original cells and have u ≥ minU.
- The ambient range is **not** excluded.

**8.0+:**
- minU = max(E_min, 1 + max(u over the *empty* set), where the max of an empty set is 0).
- Candidates are barcodes that are not original cells and have u ≥ minU. They are therefore disjoint from the ambient range.
- The 0.01×median rule is removed.
- minU is reported as the metric `<genome>_gem_group_<g>_ordmag_nonambient_threshold`.

**Both:**
- "Original cells" = the union of all current call lists for this GEM group across genomes (see §6).
- Candidates are ordered by **ascending barcode index**.
- There is no cap on the number of candidates.
- EmptyDrops returns nothing when:
  - there are no candidates;
  - the profile is empty;
  - the SGT aborted;
  - the sum of the original-cell *indices* is 0. This quirk is also hit when the only original cell is column 0.

### Observed log-likelihood (7.2/8.x)

- scipy's multinomial log-pmf over the dense candidate rows, with n = each row's sum and p = the profile.
- scipy replaces the last probability with 1 − Σ(the others) before evaluating.

### Simulation (7.2/8.x)

- 10000 simulations; jump = 1000; block size B = 1,000,000.
- **Seeding:** the global numpy RandomState is seeded with 0 **at the start of every simulation call**, i.e. once per genome × GEM group. The stage also seeds 0 on import, but that has no effect here.

Steps, in order:

1. D = distinct candidate UMI totals, ascending.
2. **Block** = a legacy weighted `choice` over feature indices, size B:
   - cumulative distribution = cumsum(p) divided by its last value, so it ends at exactly 1;
   - B uniform values are drawn;
   - each is mapped by searchsorted, **right** side.
   - Position k = 0.
3. For each simulation, in order:
   1. Draw multinomial counts with total D₀ via scipy (global RandomState; scipy's adjusted p; legacy sequential-binomial algorithm). Record the full log-pmf.
   2. For each next total D_i, with step = D_i − D_{i−1}:
      - **step ≥ 1000:** draw a fresh multinomial with total `step` from the global RandomState, add it, and **recompute** the full log-pmf at D_i.
      - **step < 20, one feature at a time:** for n = D_{i−1}+1 … D_i:
        - take feature x = block[k] and advance k;
        - **if k now equals the block length, refill at once with a full new block of B draws** and set k = 0;
        - increment that feature's count;
        - add log p_x + log(n / new count of x). This uses the unadjusted log p.
      - **20 ≤ step < 1000, batched:** if k + step ≥ block length (note the **≥**: refill even when the batch would end exactly at the block's end), **refill before using the batch**:
        - the new block = the unused tail of the old block (block length − k features, in order) followed by **k** freshly drawn features, from one weighted-choice call of size k;
        - the length stays B (B would only grow if step ≥ B, which cannot happen);
        - set k = 0.
        - Then take features [k, k + step) and advance k.
        - Group those features (sorted unique values with their multiplicities) and add:
          Σ mult·log p + lgamma(D_i + 1) − lgamma(D_{i−1} + 1) − Σ [lgamma(old + mult + 1) − lgamma(old + 1)].
      - Store the running log-likelihood for D_i.

**Is the feature sequence the same as one-at-a-time?**
- The *consumed feature stream* is continuous: nothing is skipped or reused, so a batch consumes exactly the next `step` features.
- But the *schedule of RNG calls* differs. A batch refill happens early, before consumption, and draws k values rather than B. The multinomial calls (at the start of each simulation, and for steps ≥ 1000) share the same global stream.
- So to reproduce results exactly, a reimplementation must copy this refill schedule precisely.

### p-values and BH (all versions)

- p = (1 + number of simulations at that barcode's N with simulated LL **<** observed LL) / (1 + number of simulations).
- BH adjustment:
  1. sort p descending;
  2. multiply by N_tests / rank;
  3. take the running minimum;
  4. cap at 1;
  5. map back to the original order.
  - The result does not depend on how ties in p are ordered.
- A barcode is non-ambient (called) when **adjusted p ≤ FDR**.

---

## 5. EmptyDrops simulator, 9.0+ (9.0.x, 10.0, 10.1)

### Random generator

- numpy `Generator(PCG64)` created as `default_rng(seed=42)` **when the module is imported**, so once per process. Each sample chunk is a separate process and starts fresh.
- The stream **continues** across all EmptyDrops calls in that process, in loop order (§6). Nothing else draws from it.
- A call that returns early (no candidates, SGT abort, and so on) draws nothing.

### Number of simulations

100000.

### Mapping uniforms to features

- The cumulative distribution is the sequential float64 cumsum of the profile vector in its current feature order (see "Feature order" below). It is **not normalised and not forced to end at 1**.
- Each simulation makes **one** call: `random(size = max_n)`, where max_n is the largest candidate total. Each value is a 53-bit double.
- Mapping is searchsorted, **left** side: the first index whose cumulative value is ≥ u.
- A uniform larger than the last cumulative value would give an out-of-range index and crash. There is no clamping.

### Log-likelihood per simulation

1. R_t = the number of draws among positions 1..t equal to draw t: its occurrence rank, 1-based. 10.1 computes it vectorised with identical results.
2. Increment_t = (log t − log R_t) + log p[x_t], evaluated in that order.
3. Running sum by sequential cumsum. The value for total N is the cumsum at position N − 1.
4. The simulated value for a barcode with total N is compared with its observed value using a strict "<". 10.1 accumulates the counts as it goes; 9.0/10.0 store the full matrix first. The result is the same.

### Observed log-likelihood (9.0+)

- For each candidate: lgamma(n + 1) − Σ lgamma(x + 1) + Σ x·log p, summed over the barcode's non-zero entries.
- No scipy, and no adjustment of the last probability.

### Feature order

| Version | Order used |
|---|---|
| 9.0.x | The matrix's feature-index order (evaluated features ascending) for the profile, the observed log-likelihood and the simulation cumsum |
| 10.0.0 | Immediately after the evaluated features are found and **before profile estimation**, they are permuted so their feature IDs ascend. The profile, observed-LL rows, log p and cumsum are all in that order |
| 10.1.0 | The profile and observed log-likelihood stay in index order. **Only for the simulation** is the profile vector permuted into ascending-ID order (and log p and cumsum taken from the permuted vector) |

- The ID sort is numpy's default argsort of the **byte-string IDs**: bytewise lexicographic, where a prefix sorts before its extensions. IDs are unique within a feature reference, so stability does not matter.
- 10.0 and 10.1 give the same simulation, up to last-bit floating-point differences. These can appear in the renormalisation sum when no evaluated feature is zero in the ambient set, and in the order of the observed-LL sums.

### Shared with 8.x

Ambient set, profile (with the SGT slope handling given in §4), candidates (8.0+ rule), p-value formula, BH, and the ≤ FDR test.

---

## 6. Multiple genomes and GEM groups

### Loop order

| Step | 7.2–8.0.x | 9.0+ |
|---|---|---|
| OrdMag (all versions) | for genome: for GEM group | same |
| EmptyDrops | for genome: for GEM group | **for GEM group: for genome** |
| Global minimum UMI (8.0+) | for genome: for GEM group | same |
| Mito filter (8.0+) | for genome: for GEM group | same |

### Which barcodes count as "original cells" in EmptyDrops

- **7.2/8.x:** the set is rebuilt at each genome × GEM group step from the *current* lists. So genome k treats the EmptyDrops calls already made for genomes 1…k−1 as original cells: they are excluded from candidates and, in 7.2, they enter the median.
- **9.0+:** the set is fixed once per GEM group, before the genome loop, as the union of the initial OrdMag calls across genomes.

### Features each step sees

- **OrdMag:** that genome's GEX features, plus any non-GEX features of the selected types. For example, in an antibody-only run every antibody is used for every genome.
- **EmptyDrops and global minimum UMI:** that genome's GEX features only (GEX ∩ selected types).
- **Mito:** that genome's GEX features.
- **UMI totals** are always genome-specific.
- **Expected cells** for OrdMag are passed unsplit to every genome.

### Combining and reporting

- Calls are stored per (GEM group, genome).
- The cell set = the union over genomes and samples. 7.2 sorts by (GEM group, sequence); 8.0+ uses a sorted set ordered by (GEM group, valid, sequence).
- Per-genome sets are kept. The filtered-barcodes CSV has (genome, barcode) rows, so a barcode called in two genomes appears twice. Metrics are prefixed per genome.
- The non-ambient CSV has a genome column (and a sample column from 9.0).
- In 9.0+, that CSV's "umis" column holds the **all-genome** total; in 7.2/8.x it held the genome-specific total. This affects reporting only.

---

## 7. Filters applied after calling

### Inside each sample (8.0+ only)

These run in the order shown, for every method **except MANUAL and TOP_N**:

1. **Global minimum UMI:**
   - for each genome and GEM group, keep the called barcodes whose genome-GEX UMI total is **≥ min**;
   - the default min is 0, which removes nothing;
   - set by multi `global-minimum-umis` (library-level or per sample) or the count equivalent.
2. **Mitochondrial filter:**
   - set by `max-mito-percent`, a percentage; the default 100 removes nothing;
   - a feature is mitochondrial if the part of its ID after the last "_" (to allow barnyard prefixes) is one of 26 Ensembl IDs: the 13 protein-coding MT genes of human and of mouse:
     - human ENSG00000198888, …198763, …198804, …198712, ENSG00000228253, …198899, …198938, …198840, ENSG00000212907, …198886, …198786, …198695, …198727;
     - mouse ENSMUSG00000064341, …064345, …064351, …064354, …064356, …064357, …064358, …064360, ENSMUSG00000065947, …064363, …064367, …064368, …064370;
   - mito% = 100 × mito UMIs / genome-GEX UMIs for each called barcode;
   - remove the barcode if mito% **>** threshold;
   - a 0/0 result (NaN) is never removed.

### After all samples are merged

1. **High-occupancy GEM filter**, for multiplexed Flex (§9).
2. **Zero-target-UMI removal:**
   - applies to non-spatial runs whose feature reference has target features (hybrid-capture targeted GEX *and* Flex probe sets), in all versions;
   - removes every called barcode with 0 UMIs summed over all target features;
   - the filtered matrix is then limited to target features.

7.2 has no global-minimum or mito filter. Nothing else changes which barcodes are called; spatial antibody normalisation is not a filter.

**Quirk:** in an antibody-only run with global minimum > 0, the counts are computed over "selected types ∩ genome GEX", which is empty. Every total is therefore 0 and every cell would be removed.

---

## 8. Antibody and antigen aggregate filters (same rules 7.2–10.1)

### When they run

- In each sample chunk, **before OrdMag**, on the sample-restricted raw matrix.
- Only if all of these hold:
  - the method is not MANUAL;
  - a barcode-correction table exists and is non-empty;
  - that table has Antibody or Antigen rows.
- The antibody rule runs if "Antibody Capture" is among the library types (taken from the multi graph, else the matrix). The antigen rule runs if "Antigen Capture" is.
- The barcodes found are **removed from the matrix** used by OrdMag, EmptyDrops and everything after.
- If aggregate filtering is disabled (`filter-aggregates=false`, or the count flag), detection and metrics still run but nothing is removed.

### Antibody: highly-corrected rule

A barcode is an aggregate if, in the Antibody rows of the table, umi_corrected_reads / reads **> 0.5** and reads **> 10000**.

### Antibody: count-only rule (combined with the rule above by union)

1. Build a dense barcode × antibody count matrix over **all** barcodes.
2. Drop antibodies whose total across barcodes is **< 1000**. Let n be the number of remaining ("signal") antibodies. If n < 5, this rule returns nothing.
3. K = (npb, or 1 if None) × 25.
4. Candidates = the last K barcodes of a **default argsort ascending** of each barcode's total over the signal antibodies. Boundary ties are implementation-defined.
5. For each signal antibody, its top-K set = the last K of an ascending **pandas `sort_values` (default quicksort, not stable)** of that antibody's column. Membership is checked by barcode label.
6. hits(b) = the number of signal antibodies whose top-K set contains candidate b.
7. f = 0.6 if n > 26; otherwise (−0.02·n) + 1.1, computed in float64. At n = 26 this gives 0.58.
8. b is an aggregate if hits(b) ≥ int(round(n·f)), half to even.

### Antigen (BEAM) rule

1. Per-barcode antigen totals.
2. Take the top 100 barcodes: the first 100 of a default argsort of the negated totals.
3. Q1 and Q3 = numpy's default (linear) quantiles at 0.25 and 0.75 of those 100 totals.
4. Threshold T = Q3 + 3·(Q3 − Q1).
5. If T < 1000, nothing is flagged. Otherwise flag the top-100 barcodes with total ≥ T.
6. There is no npb scaling.

---

## 9. High-occupancy GEM filter (merge step)

### When it runs

- After all samples are merged, if the primary chemistry has a probe-barcode offset that is non-null and non-zero, and the filter is not disabled (`filter-high-occupancy-gems`, default on).
- In practice: multiplexed Flex v1 and v2. It does not run for OCM, singleplex Flex or 3'/5' kits.

### Inputs

- F = all called barcodes after every per-sample filter. 7.2 sorts them by (GEM group, sequence); 8.0+ uses (GEM group, valid, sequence).
- Barcode strings have the form GEM(16) + probe(8 for v1 / 10 for v2) + "-" + GEM group. Here off = 16, the gel-bead length.
- GEM = barcode[:off]. Probe label = barcode[off:], i.e. the probe sequence plus the suffix.

### Common steps

1. G = number of distinct GEMs in F. per_gem[g] = the number of called barcodes in GEM g.
2. P = 115000 × (1/1.65) = 69696.97…
3. λ = |F| / max(int(P), G), with int(P) = 69696.
   - 7.2–9.0.x compute this as a weighted mean of a histogram whose zero bin is max(0, int(P − G)); this is equivalent.
4. Probe-label probabilities = how often each label occurs in F, in **first-appearance order** while scanning F. That order changes the mapping of uniforms to labels and so the result.
5. RNG: `np.random.seed(0)` on the global legacy RandomState, called inside the threshold routine.
6. Simulate 1,000,000 GEM occupancies (distribution in the table below). Keep those > 0 and let m = their maximum.
7. Draw a legacy weighted choice of labels with shape (number kept, m). The cumsum is normalised and mapped with the right-side search.
8. For each kept GEM: its distinct-label count = the number of distinct labels among the first (occupancy) entries of its row.

### Threshold T

| Version | Occupancy distribution | Threshold rule |
|---|---|---|
| 7.2–9.0.x | Poisson(λ) | T = int(ceil(0.999 quantile of the distinct-label counts)), using numpy's default linear quantile |
| 10.0+ | Negative binomial with n = λ, p = 0.5 (legacy gamma–Poisson) | w[k] = k × (number of simulated GEMs with k distinct labels), for k = 0..max. T = the smallest k whose cumulative Σw / total Σw is **≥ 0.999** |

- If λ = 0, T = 0.

### What is removed

Every called barcode whose GEM has per_gem **> T**, from both the overall list and every per-genome list. Per-barcode read metrics are used only for reporting.

---

## 10. Method selection and special modes

### Precedence (per sample chunk)

1. A cell-barcodes JSON is configured (multi `barcode-sample-assignment`) → **MANUAL**.
2. Force-cells is set (library-level, or the sample's value) → **TOP_N**.
3. An override mode is set → that method. Accepted strings: manual, topn, ordmag, ordmag_nonambient, gradient, targeted.
   - 7.2/8.x: always null in the shipped pipelines.
   - 9.0.x/10.0: comes from an internal pipeline parameter, not a user CSV option.
   - **10.1:** the hidden `[gene-expression] cell-calling-mode` option, values `ordmag | ordmag_nonambient | gradient`.
4. Antibody-only run → **ORDMAG**.
5. 7.2–9.0.x only: hybrid-capture targeted (§11) → **TARGETED**.
6. Otherwise → **ORDMAG_NONAMBIENT** (OrdMag followed by EmptyDrops).

### Feature types

- Gene Expression, or Antibody Capture for antibody-only runs.
- The library-type override is always null in the shipped pipelines.

### What each method does

| Method | What it does |
|---|---|
| MANUAL | For every genome, the listed barcodes present anywhere in the chunk's matrix. No aggregate filter, EmptyDrops, global minimum or mito filter. The merge-step filters still apply |
| TOP_N | Runs after the aggregate filter. Takes min(force // G, number of non-zero barcodes) in stable-descending order. No EmptyDrops, global minimum or mito filter |
| ORDMAG | OrdMag only, then the global minimum and mito filters. This is what antibody-only runs get, with the OrdMag grid capped by the chemistry's `low` |
| GRADIENT | §11, with expected cells or 3000, max-add 20000 and minimum 10, then the global minimum and mito filters |
| TARGETED | GRADIENT followed by the target-UMI ≥ 1 test (§11) |

10.1's `cell-calling-mode` deliberately offers no "targeted" option; the Rust documentation calls it dead code.

---

## 11. Gradient and targeted method: corrections to your summary

1. **N.** Expected cells (passed through `int(float/G)` and then `// G`, as in §3) or 3000 if not supplied. The floor of 50 applies only to a supplied value; 3000 is already above it. ✔
2. **Ordering.** s = non-zero counts in **descending** order (a stable sort of values). base = s[min(round(N·q), len − 1)], where q is the float 0.010000000000000009, so exact .5 cases round up. ✔ with that caveat.
3. **lower.** lower = (number of s ≥ base / 10.0) − 1, then min(lower, len − 1). This uses **unrounded float division**, unlike OrdMag's rounded cutoff. ✔
4. **upper.** upper = min(lower + maxAdd, number of s ≥ 10), then max(upper, lower), then min(upper, len − 1). maxAdd = 20000. The value 2000 applies **only in TARGETED**, and only when the description is LT (7.2–10.0); GRADIENT always uses 20000. ✔
5. **Points.** For each unique u in descending order: (log10(number of s ≥ u), log10 u), plus a final point (log10(1 + Σs), 0). The x values strictly increase. ✔
6. **Spline.**
   - Degree k = min(3, P − 1), where P = the number of points.
   - First fit an interpolating spline with smoothing 0 (scipy `UnivariateSpline`, s = 0).
   - If P > 50: work out K from P:

     | P | K |
     |---|---|
     | < 200 | int(2^(log2 50 + (log2 100 − log2 50)·(P − 50)/150)) |
     | < 800 | int(2^(log2 100 + (log2 140 − log2 100)·(P − 200)/600)) |
     | < 3200 | int(2^(log2 140 + (log2 200 − log2 140)·(P − 800)/2400)) |
     | otherwise | int(200 + (P − 3200)^0.2) |

     (P < 50 would give K = P, but this is never used.)
   - If K is smaller than the number of knots the interpolating spline returns, take K − 2 of those knots, at the integer positions floor(linspace(1, knot count − 2, K − 2)), and refit a least-squares spline (`LSQUnivariateSpline`) of the same degree with those as interior knots.
7. **Derivative.** Evaluate the first derivative at every point **except the appended end-point**. Set it to 0 where the **raw** x (count of s ≥ u, not its log) falls outside [lower, upper], inclusive at both ends. ✔
8. **Cut.** cut = round(10^y), where y is taken at the first minimum of the derivative. ✔ If every value is 0, the first point is used, so cut = the maximum count.
9. **Number called.** n = max(number of s **>** cut (strict), lower + 1), then top_n = min(n, len(s)). ✔
10. **Selection.** Stable-descending order over all barcodes (ties go to the higher index), then sort the indices. Gradient counts are over the genome's GEX features, including non-targeted genes.
11. **TARGETED only.** Then keep the barcodes with ≥ 1 UMI in the genome's target features.
12. **When TARGETED is chosen automatically.** 7.2.0–9.0.1, when all of these hold:
    - the feature reference has target features;
    - the target-set CSV's targeting method is **hybrid capture** (RTL and Flex probe sets are excluded);
    - none of manual, force, override or antibody-only applies.

    It is **never chosen automatically from 10.0** (only by override).

---

## 12. Not determined, or needing your attention

1. **Bundled numpy and scipy versions.** These are not pinned anywhere in the repository. They decide the tie order of numpy's default unstable argsort, and hence:
   - the barcodes at the two boundaries of the ambient range;
   - antibody top-K candidates and per-antibody top-K sets (pandas quicksort);
   - the antigen top 100.

   They also decide scipy's multinomial parameter handling in 7.2/8.x. You can only match exactly by using the same numpy build and CPU dispatch path, or by checking the numpy shipped in a Cell Ranger installation.
2. **Raw-matrix column order.** I inferred (GEM group, valid, sequence) from the Rust barcode ordering; I did not trace the matrix writer end to end.
3. **Library-level expect/force with multiplexed samples.** I did not check whether the multi-config validator forbids library-level expect-cells or force-cells for multiplexed Flex or OCM. Only the precedence rule is described.
4. **Uniform beyond the last cumulative value (9.0+).** The draw would index out of range. I inferred a crash from the code; I did not run it.
5. **8.0.x overhang chemistries.** I could not tell whether SC3Pv4-OH and its relatives were customer-facing.
