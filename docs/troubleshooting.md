---
title: Troubleshooting
nav_order: 7
---

# Troubleshooting
{: .no_toc }

1. TOC
{:toc}

## "this looks like a filtered matrix"

The input has few barcodes and almost none with fewer than 10 counts: it holds cells some tool
already called. Knee methods need the empty droplets too. Give the raw matrix
(`raw_feature_bc_matrix.h5`, `Solo.out/Gene/raw/`, `counts_unfiltered/`).

## "no raw count matrix found in ..."

The folder holds none of the layouts in [Inputs](reference/inputs.md). If it says "only a
filtered matrix", the pipeline's raw output is elsewhere or was not kept. For STARsolo, check
`--solo-feature` names a folder that exists under `Solo.out/`.

## "no features of type 'Gene Expression'"

The matrix has a feature-type column without Gene Expression in it (an antibody-only library,
for example). Name the type to count with `--feature-type`, or count everything with
`--feature-type all`.

## "... cannot be reproduced"

The release asked for falls in a range whose code is not public (Cell Ranger 4.0 to 7.1). The
method's page says what is known about it from the release notes.

## "... needs --n-cells" (or another option)

The method has a parameter without a default. `sc-droplets methods <id>` lists them; in
`compare` and `plot`, give it in the spec: `-m top-n:n_cells=3000`.
