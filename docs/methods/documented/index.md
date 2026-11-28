---
title: Documented only
parent: Methods
nav_order: 16
has_children: true
---

# Methods documented only

Methods described here but not run by the package: they need a trained model, or their code is
closed or carries no licence. Each page says what the method does, release by release, and how
its own output can still be set beside the package's calls.

| Method | Tool | Why documented only |
|---|---|---|
| [CellBender model](../cellbender/model.md) | CellBender | a trained variational model; needs a GPU to be practical |
| [EmptyNN](emptynn.md) | EmptyNN | a trained neural network; no licence |
| [dropkick](dropkick.md) | dropkick | a regularised logistic regression trained on its own labels |
| [BD Rhapsody](bd-rhapsody.md) | BD Rhapsody pipeline | closed source |
| [ScaleBio](scalebio.md) | ScaleBio ScaleRna | its own rules, read from its pipeline |
| [CeleScope](celescope.md) | Singleron CeleScope | runs STARsolo's EmptyDrops_CR since 2.0 |
| [Ambimux](ambimux.md) | Ambimux | genotypes and a trained model, on multiome |
| [AmbientMapper](ambientmapper.md) | AmbientMapper | genotypes, empty / singlet / doublet by model fits |
