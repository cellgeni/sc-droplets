---
title: Planned
parent: Methods
nav_order: 15
has_children: true
---

# Planned methods

Methods the package does not run yet, scheduled for a release. Each page describes the method
as its source does it, release by release, and what reproducing it will need. When a method is
implemented, its page moves to its tool and gains *Run it*, *Differences from the original* and
*Validation*.

## 0.2: tests and classifiers on RNA

| Method | Tool | Idea |
|---|---|---|
| [CB2](cb2.md) | scCB2 | cluster the candidates, test clusters and single barcodes against the ambient profile |
| [DIEM](diem.md) | DIEM | a mixture of debris and cell profiles, fitted to the barcodes above a floor |
| [DropletQC](dropletqc.md) | DropletQC | the nuclear fraction separates empty droplets, damaged cells and cells |
| [DropSift](dropsift.md) | DropSift | a classifier trained on exemplar cells and empty droplets |
| [valiDrops](validrops.md) | valiDrops | a curve threshold, then quality measures and clustering |
| [SiftCell](siftcell.md) | SiftCell | barcodes shuffled into a null, then a boosted classifier |
| [miQC](miqc.md) | miQC | a mixture of regressions of the mitochondrial fraction on genes detected |
| [dropEst classifier](../dropest/classifier.md) | dropEst | a classifier on quality measures, trained on the bounds' high and low classes |
| [alevin naive Bayes whitelist](../alevin/naive-bayes.md) | alevin | a naive Bayes classifier on the initial whitelist's cells and low-confidence band |
