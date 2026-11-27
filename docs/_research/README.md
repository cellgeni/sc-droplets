# Research notes

Working notes from the source reading done before the first release (November 2026). They are
the starting point for the method pages, not part of the site (`_config.yml` excludes this
folder).

| File | Covers |
|---|---|
| [cellranger-starsolo.md](cellranger-starsolo.md) | Cell Ranger 1.x–10.1, Cell Ranger ATAC and ARC, STARsolo 2.7.0a–2.7.11b |
| [dropletutils-cellbender-cellsweep.md](dropletutils-cellbender-cellsweep.md) | DropletUtils 1.0–1.34, CellBender 0.1–0.4, cellsweep 0.1–1.0.2, existing Python ports |
| [other-methods-and-reviews.md](other-methods-and-reviews.md) | Reviews and benchmarks, every other knee heuristic, statistical and ML methods, mixture models |
| [cellranger-4.0-7.1-outputs.md](cellranger-4.0-7.1-outputs.md) | Cell Ranger 4.0–7.1 from 10x's published outputs (no source), with `cellranger-4.0-7.1-survey.py` |
| [planned-rna.md](planned-rna.md) | The 0.2 methods: CB2, DIEM, DropletQC, DropSift, valiDrops, SiftCell, miQC, dropEst's classifier, alevin's naive Bayes whitelist |
| [planned-atac.md](planned-atac.md) | The 0.3 methods: Cell Ranger ATAC and ARC, EmptyDropsMultiome, BAP, SADE |
| [docs-only.md](docs-only.md) | Methods documented, not implemented: CellBender's model, EmptyNN, dropkick, BD Rhapsody, ScaleBio, CeleScope, Ambimux, AmbientMapper |

Tags used in the notes:

- **[SRC]**: read in the tool's source at the tag named.
- **[DOC]**: documentation or release notes only; no source.
- **[OUT]**: inferred from the tool's published outputs (no source).
- **UNVERIFIED**: inferred; nothing to check it against.

When a method page is written, the page becomes the specification and these notes are only
history. If a note and a page disagree, the page wins once its validation against the original
tool has passed.
